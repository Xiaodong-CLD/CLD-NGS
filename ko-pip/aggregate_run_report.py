"""Per-run aggregation of per-sample outputs into final_report/.

Reads each sample's persistent Results/<sample>/ artifacts and concatenates
them into run-level tables:

    final_report/all_snv.tsv            VarScan SNVs (var_freq > 0.01), clone-annotated
    final_report/all_indel.tsv          merged GATK+VarScan indels, clone-annotated
    final_report/all_uninformative.tsv  per-sample uninformative-read summaries
    final_report/all_breakpoints.tsv    per-locus aggregated breakpoint clusters
                                        (summed support > 0.01), clone-annotated
    final_report/all_consensus.fa       per-sample consensus haplotype FASTAs,
                                        concatenated with clone-annotated headers

`sample_records` is a list of (folder_name, sample_id, clone_id) tuples; the
folder_name is the Results/<folder_name>/ directory for that sample.
"""

import os

import pandas as pd

MIN_SNV_VAR_FREQ = 0.01  # keep SNVs with VarScan var_freq strictly above 1%

# Uninformative-site keep-gate for all_uninformative.tsv. A site's uninformative
# reads (sum of ALT+REF+uncountable) must make up at least this fraction of the
# sample's total reads (sample_read_frac); low-support sites are dropped.
MIN_UNINFORMATIVE_SAMPLE_FRAC = 0.01  # keep site if sample_read_frac >= 0.01

# Per-locus aggregation of soft-clip breakpoint clusters for all_breakpoints.tsv.
# Sub-clusters (split by sequence-level edit distance) that sit at the same
# locus are merged: same side (3p/5p), top_locus position within +/- POS_TOL bp,
# and soft-clip length within +/- LEN_TOL bp. The most frequent (highest-reads)
# sub-cluster represents the merged locus; support_frac is summed across members.
# A merged locus is kept only if its summed support_frac exceeds MIN_SUMMED_FRAC.
BREAKPOINT_POS_TOL_BP = 3        # +/- bp on top_locus position to call same locus
BREAKPOINT_LEN_TOL_BP = 3        # +/- bp on soft-clip length to call same locus
BREAKPOINT_MIN_SUMMED_FRAC = 0.01  # keep merged locus if summed support_frac > 0.01

# Per-caller attributed-depth keep-gate for all_indel.tsv. A row is kept if its
# caller attributed at least this fraction of the raw locus depth to an allele:
#   GATK    : gatk_sum_ad_frac    = sum_AD / DP_INFO  (INFO/DP)
#   VarScan : varscan2_sum_ad_frac = (AD + RD) / SDP
# A 'both' row is kept if either caller passes.
GATK_GATE_SUM_AD_FRAC = 0.01      # GATK keep if sum_AD/DP_INFO >= 0.01
VARSCAN_GATE_SUM_AD_FRAC = 0.01   # VarScan keep if (AD+RD)/SDP >= 0.01

# Consensus corroboration: a haplotype's internal indel (from its sig) is
# confirmed when a gate-passing indel of the same type and length sits within
# this many bp of it, or a gated breakpoint locus sits within this many bp.
# The tolerance absorbs left-alignment / anchor-convention shifts in
# low-complexity tracts (e.g. sig 72-75:DEL:TTTG vs indel pos 71 DEL 4).
CONSENSUS_CONFIRM_POS_TOL_BP = 5


def _read_tsv(path):
    """Read a TSV into a DataFrame, or return None if missing/empty."""
    if not os.path.exists(path):
        return None
    try:
        return pd.read_csv(path, sep='\t')
    except pd.errors.EmptyDataError:
        return None


def aggregate_snv(results_dir, sample_records, output_path):
    """Write per-sample VarScan SNV blocks (var_freq > 0.01) in banner-block form.

    Mirrors all_indel.tsv: a global note, a single column header (with leading
    sample/clone_id), then one '# sample' banner per sample followed by that
    sample's kept rows (var_freq > 0.01, sorted by freq).
    """
    header_cols = None
    blocks = []  # (folder_name, clone_id, kept_df_or_None)
    fallback_cols = [
        'sample_id', 'chrom', 'position', 'ref', 'alt', 'consensus',
        'ref_support', 'alt_support', 'total_informative', 'var_freq',
        'p_value', 'alt_strands_plus', 'alt_strands_minus',
    ]
    for folder_name, _sample_id, clone_id in sample_records:
        df = _read_tsv(os.path.join(results_dir, folder_name, 'snv', 'varscan_snp_report.tsv'))
        if df is None or df.empty:
            blocks.append((folder_name, clone_id, None))
            continue
        df = df.drop(columns=['sample', 'clone_id'], errors='ignore')
        if 'var_freq' in df.columns:
            df = df[pd.to_numeric(df['var_freq'], errors='coerce').fillna(0.0) > MIN_SNV_VAR_FREQ]
            df = df.sort_values('var_freq', ascending=False)
        if header_cols is None and not df.empty:
            header_cols = list(df.columns)
        blocks.append((folder_name, clone_id, df))

    if header_cols is None:
        header_cols = fallback_cols

    n_rows = 0
    with open(output_path, 'w') as fh:
        fh.write(f'# all_snv: per-sample VarScan SNVs with var_freq > {MIN_SNV_VAR_FREQ:g}\n')
        fh.write('sample\tclone_id\t' + '\t'.join(header_cols) + '\n')
        for folder_name, clone_id, df in blocks:
            fh.write(f'# sample: {folder_name}\n')
            if not (df is None or df.empty):
                for _, r in df.iterrows():
                    vals = []
                    for c in header_cols:
                        v = r.get(c, '')
                        vals.append('' if pd.isna(v) else str(v))
                    fh.write(f'{folder_name}\t{clone_id}\t' + '\t'.join(vals) + '\n')
                    n_rows += 1
            else:
                fh.write(f'No SNVs with var_freq > {MIN_SNV_VAR_FREQ:g} detected\n')
            fh.write('\n')
    print(f'[Takeda Pipeline Message] all_snv wrote {n_rows} rows -> {output_path}')


def _parse_total_mapped_reads(path):
    """Read the '# total_mapped_reads: N' banner from a merged_indel.tsv, or None."""
    try:
        with open(path) as fh:
            for line in fh:
                if not line.startswith('#'):
                    break
                if 'total_mapped_reads' in line:
                    return line.split(':', 1)[1].strip()
    except OSError:
        return None
    return None


def _parse_comment_field(path, key):
    """Read a '# <key>: value' comment banner from a TSV's header, or None."""
    try:
        with open(path) as fh:
            for line in fh:
                if not line.startswith('#'):
                    break
                if key in line:
                    return line.split(':', 1)[1].strip()
    except OSError:
        return None
    return None


def _indel_gate_pass(row):
    """True if a merged-indel row clears the per-caller attributed-depth gate.

    GATK rows judged on gatk_sum_ad_frac, VarScan rows on varscan2_sum_ad_frac;
    a 'both' row passes if either caller clears its threshold. Rows flagged
    edge_artifact='yes' (indel start within a few bp of an amplicon terminus)
    are dropped outright -- they are excluded from all_indel.tsv and never
    corroborate a consensus haplotype.
    """
    if str(row.get('edge_artifact', 'no')) == 'yes':
        return False
    source = str(row.get('source', ''))

    def _f(col):
        try:
            return float(row.get(col, ''))
        except (TypeError, ValueError):
            return 0.0

    gatk_pass = source in ('gatk_only', 'both') and _f('gatk_sum_ad_frac') >= GATK_GATE_SUM_AD_FRAC
    varscan_pass = source in ('varscan2_only', 'both') and _f('varscan2_sum_ad_frac') >= VARSCAN_GATE_SUM_AD_FRAC
    return gatk_pass or varscan_pass


def aggregate_indel(results_dir, sample_records, output_path):
    """Write per-sample merged-indel blocks, keeping only gate-passing rows.

    Output is a concatenation of per-sample blocks. Each block is preceded by a
    '# sample' / '# total_mapped_reads' banner; the rows are the sample's
    merged_indel.tsv rows that pass the per-caller attributed-depth gate. A
    single column header (with leading sample/clone_id) is written once.
    """
    header_cols = None
    blocks = []  # (folder_name, clone_id, total_mapped_reads, kept_df_or_None)
    for folder_name, _sample_id, clone_id in sample_records:
        path = os.path.join(results_dir, folder_name, 'indel', 'merged_indel.tsv')
        total = _parse_total_mapped_reads(path) if os.path.exists(path) else None
        df = None
        if os.path.exists(path):
            try:
                df = pd.read_csv(path, sep='\t', comment='#')
            except pd.errors.EmptyDataError:
                df = None
        if df is None or df.empty:
            blocks.append((folder_name, clone_id, total, None))
            continue
        keep = df[df.apply(_indel_gate_pass, axis=1)].copy()
        if header_cols is None and not keep.empty:
            header_cols = list(keep.columns)
        blocks.append((folder_name, clone_id, total, keep))

    n_rows = 0
    with open(output_path, 'w') as fh:
        fh.write('# all_indel: per-sample merged GATK+VarScan indels passing the attributed-depth gate\n')
        fh.write(f'# gate: keep row if GATK sum_AD/DP_INFO >= {GATK_GATE_SUM_AD_FRAC:g} '
                 f'OR VarScan (AD+RD)/SDP >= {VARSCAN_GATE_SUM_AD_FRAC:g}; '
                 'edge_artifact rows dropped\n')
        if header_cols is not None:
            fh.write('sample\tclone_id\t' + '\t'.join(header_cols) + '\n')
        for folder_name, clone_id, total, keep in blocks:
            fh.write(f'# sample: {folder_name}\n')
            fh.write(f'# total_mapped_reads: {total if total is not None else "NA"}\n')
            if not (keep is None or keep.empty or header_cols is None):
                for _, r in keep.iterrows():
                    vals = [str(r.get(c, '')) for c in header_cols]
                    fh.write(f'{folder_name}\t{clone_id}\t' + '\t'.join(vals) + '\n')
                    n_rows += 1
            else:
                fh.write(f'No merged indels passing the gate detected\n')
            fh.write('\n')
    print(f'[Takeda Pipeline Message] all_indel wrote {n_rows} rows -> {output_path}')


def aggregate_uninformative(results_dir, sample_records, output_path):
    """Write per-sample uninformative-read summaries in banner-block form.

    Mirrors all_indel.tsv: a global note, a single column header (with leading
    sample/clone_id), then one '# sample' banner per sample followed by that
    sample's uninformative_summary.tsv rows.
    """
    header_cols = None
    blocks = []  # (folder_name, clone_id, df_or_None)
    for folder_name, _sample_id, clone_id in sample_records:
        df = _read_tsv(os.path.join(results_dir, folder_name,
                                    'uninformative_analysis', 'uninformative_summary.tsv'))
        if df is None or df.empty:
            blocks.append((folder_name, clone_id, None))
            continue
        df = df.drop(columns=['sample', 'clone_id'], errors='ignore')
        if header_cols is None:
            header_cols = list(df.columns)
        blocks.append((folder_name, clone_id, df))

    n_rows = 0
    with open(output_path, 'w') as fh:
        fh.write('# all_uninformative: per-sample composition of reads GATK ignored at its called sites\n')
        fh.write(f'# gate: keep site if sample_read_frac >= {MIN_UNINFORMATIVE_SAMPLE_FRAC:g} '
                 '(uninformative reads make up >= 1% of the sample); low-support sites dropped\n')
        if header_cols is not None:
            fh.write('sample\tclone_id\t' + '\t'.join(header_cols) + '\n')
        for folder_name, clone_id, df in blocks:
            fh.write(f'# sample: {folder_name}\n')
            if not (df is None or df.empty or header_cols is None):
                for _, r in df.iterrows():
                    try:
                        sample_read_frac = float(r.get('sample_read_frac', 0) or 0)
                    except (TypeError, ValueError):
                        sample_read_frac = 0.0
                    if sample_read_frac < MIN_UNINFORMATIVE_SAMPLE_FRAC:
                        continue
                    vals = [str(r.get(c, '')) for c in header_cols]
                    fh.write(f'{folder_name}\t{clone_id}\t' + '\t'.join(vals) + '\n')
                    n_rows += 1
            else:
                fh.write(f'No uninformative sites with sample_read_frac >= {MIN_UNINFORMATIVE_SAMPLE_FRAC:g} detected\n')
            fh.write('\n')
    print(f'[Takeda Pipeline Message] all_uninformative wrote {n_rows} rows -> {output_path}')


BREAKPOINT_AGG_COLUMNS = [
    'side', 'anchor', 'length', 'summed_support_frac', 'rep_reads',
    'rep_sequence', 'blast_result', 'flag', 'variant',
]


def _summarize_blast(construct_blast, reference_blast):
    """Summarize construct/reference BLAST hits into one field.

    Returns '.' when neither matched, 'construct' or 'reference' for a single
    match, or 'construct,reference' when a locus matched both.
    """
    hits = []
    if str(construct_blast) == 'construct_match':
        hits.append('construct')
    if str(reference_blast) == 'reference_match':
        hits.append('reference')
    return ','.join(hits) if hits else '.'


def _aggregate_breakpoint_clusters(df):
    """Merge per-locus soft-clip sub-clusters into representative rows.

    Sub-clusters are merged when they share the same side and chrom and sit
    within BREAKPOINT_POS_TOL_BP / BREAKPOINT_LEN_TOL_BP of an existing group's
    representative. Rows are processed in descending read support so the first
    (most frequent) sub-cluster seeds and represents each group. A group is
    returned only if its summed support_frac exceeds BREAKPOINT_MIN_SUMMED_FRAC
    AND its representative is construct_blast == 'construct_match' OR
    reference_blast == 'reference_match' OR flag contains large_insert.
    """
    rows = df.to_dict('records')
    rows.sort(key=lambda r: float(r.get('cluster_reads', 0) or 0), reverse=True)
    groups = []  # each: dict with rep fields + summed_support_frac
    for r in rows:
        try:
            pos = int(float(r.get('anchor', 0) or 0))
        except (TypeError, ValueError):
            pos = 0
        side = str(r.get('side'))
        length = int(float(r.get('length', 0) or 0))
        support_frac = float(r.get('support_frac', 0) or 0)
        placed = False
        for g in groups:
            if (g['side'] == side
                    and abs(pos - g['pos']) <= BREAKPOINT_POS_TOL_BP
                    and abs(length - g['length']) <= BREAKPOINT_LEN_TOL_BP):
                g['summed_support_frac'] += support_frac
                placed = True
                break
        if not placed:
            groups.append({
                'pos': pos,
                'side': side,
                'length': length,
                'rep_reads': int(float(r.get('cluster_reads', 0) or 0)),
                'rep_sequence': r.get('rep_sequence'),
                'construct_blast': r.get('construct_blast'),
                'reference_blast': r.get('reference_blast', '.'),
                'flag': str(r.get('flag', '.')),
                'variant': str(r.get('variant', '.')),
                'summed_support_frac': support_frac,
            })
    kept = []
    for g in groups:
        if g['summed_support_frac'] <= BREAKPOINT_MIN_SUMMED_FRAC:
            continue
        # Keep locus if it matches the known construct, matches the reference
        # (suggesting a large deletion), or looks like a large insertion.
        flag = g.get('flag', '.')
        if not (
            str(g['construct_blast']) == 'construct_match'
            or str(g['reference_blast']) == 'reference_match'
            or 'large_insert' in str(flag)
        ):
            continue
        kept.append({
            'side': g['side'],
            'anchor': g['pos'],
            'length': g['length'],
            'summed_support_frac': round(g['summed_support_frac'], 6),
            'rep_reads': g['rep_reads'],
            'rep_sequence': g['rep_sequence'],
            'blast_result': _summarize_blast(g['construct_blast'], g['reference_blast']),
            'flag': flag,
            'variant': g['variant'],
        })
    kept.sort(key=lambda g: g['summed_support_frac'], reverse=True)
    return kept


def aggregate_breakpoints(results_dir, sample_records, output_path):
    """Write per-sample aggregated breakpoint loci (summed support > 0.01).

    Reads breakpoint_analysis/top_clusters.tsv per sample (which carries '#'
    comment metadata lines above the header), merges sub-clusters that map to
    the same locus, and keeps merged loci that clear the support gate and are a
    construct match or a large insertion. Output mirrors all_indel.tsv: a
    global gate note, a single column header, then one '# sample' /
    '# total_analyzed_reads' banner block per sample followed by that sample's
    kept rows.
    """
    blocks = []  # (folder_name, clone_id, total_analyzed_reads, kept_rows)
    for folder_name, _sample_id, clone_id in sample_records:
        path = os.path.join(results_dir, folder_name,
                            'breakpoint_analysis', 'top_clusters.tsv')
        total = _parse_comment_field(path, 'total_analyzed_reads') if os.path.exists(path) else None
        kept = []
        if os.path.exists(path):
            try:
                df = pd.read_csv(path, sep='\t', comment='#')
            except pd.errors.EmptyDataError:
                df = None
            if df is not None and not df.empty and 'anchor' in df.columns:
                kept = _aggregate_breakpoint_clusters(df)
        blocks.append((folder_name, clone_id, total, kept))

    n_rows = 0
    with open(output_path, 'w') as fh:
        fh.write('# all_breakpoints: per-sample soft-clip loci aggregated within '
                 f'+/-{BREAKPOINT_POS_TOL_BP} bp position and +/-{BREAKPOINT_LEN_TOL_BP} bp length\n')
        fh.write(f'# gate: keep locus if summed support_frac > {BREAKPOINT_MIN_SUMMED_FRAC:g} '
                 'AND (construct_blast == match OR reference_blast == match OR flag contains large_insert); '
                 'most frequent sub-cluster represents the locus\n')
        fh.write('sample\tclone_id\t' + '\t'.join(BREAKPOINT_AGG_COLUMNS) + '\n')
        for folder_name, clone_id, total, kept in blocks:
            fh.write(f'# sample: {folder_name}\n')
            fh.write(f'# total_analyzed_reads: {total if total is not None else "NA"}\n')
            if kept:
                for g in kept:
                    vals = [str(g.get(c, '')) for c in BREAKPOINT_AGG_COLUMNS]
                    fh.write(f'{folder_name}\t{clone_id}\t' + '\t'.join(vals) + '\n')
                    n_rows += 1
            else:
                fh.write('No large indel-like breakpoint detected\n')
            fh.write('\n')
    print(f'[Takeda Pipeline Message] all_breakpoints wrote {n_rows} rows -> {output_path}')


def _gated_indel_calls(results_dir, folder_name):
    """Return [(pos, type, length)] for merged_indel rows passing the indel gate."""
    path = os.path.join(results_dir, folder_name, 'indel', 'merged_indel.tsv')
    if not os.path.exists(path):
        return []
    try:
        df = pd.read_csv(path, sep='\t', comment='#')
    except pd.errors.EmptyDataError:
        return []
    if df is None or df.empty:
        return []
    calls = []
    for _, r in df[df.apply(_indel_gate_pass, axis=1)].iterrows():
        try:
            pos = int(float(r.get('start_pos')))
            length = int(float(r.get('length')))
        except (TypeError, ValueError):
            continue
        calls.append((pos, str(r.get('type')), length))
    return calls


def _gated_breakpoint_positions(results_dir, folder_name):
    """Return [pos] for gate-passing breakpoint loci (soft-clip junctions)."""
    path = os.path.join(results_dir, folder_name,
                        'breakpoint_analysis', 'top_clusters.tsv')
    if not os.path.exists(path):
        return []
    try:
        df = pd.read_csv(path, sep='\t', comment='#')
    except pd.errors.EmptyDataError:
        return []
    if df is None or df.empty or 'anchor' not in df.columns:
        return []
    positions = []
    for g in _aggregate_breakpoint_clusters(df):
        positions.append(g['anchor'])
    return positions


def _parse_sig_indels(sig):
    """Yield (start, type, length) for the DEL/INS variants in a haplotype sig.

    A sig is '|'-joined 'start-end:TYPE:allele' parts; 'WT' has none. SNV/MNV
    parts are ignored (only indels are cross-checked against indel/breakpoint
    gates).
    """
    variants = []
    if not sig or sig == 'WT':
        return variants
    for part in sig.split('|'):
        fields = part.split(':')
        if len(fields) < 3:
            continue
        span, vtype, allele = fields[0], fields[1], fields[2]
        if vtype not in ('DEL', 'INS'):
            continue
        try:
            start = int(span.split('-')[0])
        except ValueError:
            continue
        variants.append((start, vtype, len(allele)))
    return variants


def _variant_confirmed(variant, gated_indels, gated_bp_positions):
    """True if a haplotype indel is backed by a gated indel or breakpoint locus."""
    start, vtype, length = variant
    for pos, gtype, glength in gated_indels:
        if gtype == vtype and glength == length and abs(pos - start) <= CONSENSUS_CONFIRM_POS_TOL_BP:
            return True
    for pos in gated_bp_positions:
        if abs(pos - start) <= CONSENSUS_CONFIRM_POS_TOL_BP:
            return True
    return False


def _apply_unconfirmed_flag(header_content, gated_indels, gated_bp_positions):
    """Append 'unconfirmed' to a haplotype header's flags= field when its indel
    variant(s) are not corroborated by a gated indel/breakpoint call.

    header_content is the header text without the leading '>'. Reference lines
    (no sig=/flags=) are returned unchanged.
    """
    if 'sig=' not in header_content or 'flags=' not in header_content:
        return header_content
    parts = [p.strip() for p in header_content.split(';')]
    sig = ''
    for p in parts:
        if p.startswith('sig='):
            sig = p[len('sig='):]
    variants = _parse_sig_indels(sig)
    unconfirmed = bool(variants) and any(
        not _variant_confirmed(v, gated_indels, gated_bp_positions) for v in variants
    )
    if not unconfirmed:
        return header_content
    new_parts = []
    for p in parts:
        if p.startswith('flags='):
            existing = [x for x in p[len('flags='):].split(',') if x and x != 'none']
            if 'unconfirmed' not in existing:
                existing.append('unconfirmed')
            p = 'flags=' + ','.join(existing)
        new_parts.append(p)
    return '; '.join(new_parts)


def aggregate_consensus(results_dir, sample_records, output_path):
    """Concatenate per-sample consensus haplotype FASTAs into one run FASTA.

    Reads each sample's consensus/<folder>.consensus.fa and re-emits every
    record with the sample folder (and clone id) prefixed onto the header, so
    the run-level file stays a valid FASTA where each haplotype is traceable to
    its sample. A haplotype whose indel (from its sig) is not backed by a
    gate-passing indel or breakpoint locus for that sample gets 'unconfirmed'
    appended to its flags= field.
    """
    n_records = 0
    n_samples = 0
    with open(output_path, 'w') as out:
        for folder_name, _sample_id, clone_id in sample_records:
            path = os.path.join(results_dir, folder_name, 'consensus',
                                f'{folder_name}.consensus.fa')
            if not os.path.exists(path):
                continue
            n_samples += 1
            gated_indels = _gated_indel_calls(results_dir, folder_name)
            gated_bp_positions = _gated_breakpoint_positions(results_dir, folder_name)
            clone_tag = f' clone={clone_id}' if clone_id else ''
            with open(path) as fh:
                for line in fh:
                    if line.startswith('>'):
                        content = _apply_unconfirmed_flag(
                            line[1:].rstrip('\n'), gated_indels, gated_bp_positions
                        )
                        out.write(f'>{folder_name}{clone_tag} {content}\n')
                        n_records += 1
                    else:
                        out.write(line)
    print(f'[Takeda Pipeline Message] all_consensus wrote {n_records} record(s) '
          f'from {n_samples} sample(s) -> {output_path}')


def aggregate_run_report(results_dir, final_report_dir, sample_records):
    """Write the run-level tables under final_report_dir."""
    import shutil
    os.makedirs(final_report_dir, exist_ok=True)
    aggregate_snv(results_dir, sample_records, os.path.join(final_report_dir, 'all_snv.tsv'))
    aggregate_indel(results_dir, sample_records, os.path.join(final_report_dir, 'all_indel.tsv'))
    aggregate_uninformative(results_dir, sample_records,
                            os.path.join(final_report_dir, 'all_uninformative.tsv'))
    aggregate_breakpoints(results_dir, sample_records,
                          os.path.join(final_report_dir, 'all_breakpoints.tsv'))
    aggregate_consensus(results_dir, sample_records,
                        os.path.join(final_report_dir, 'all_consensus.fa'))

    # Collect per-sample coverage plots into final_report/coverage_plots/
    coverage_plots_dir = os.path.join(final_report_dir, 'coverage_plots')
    os.makedirs(coverage_plots_dir, exist_ok=True)
    for folder_name, _sid, _clone in sample_records:
        src = os.path.join(results_dir, folder_name, 'plots', f'{folder_name}_coverage_plot.jpeg')
        if os.path.exists(src):
            shutil.copy2(src, coverage_plots_dir)


if __name__ == '__main__':
    pass
