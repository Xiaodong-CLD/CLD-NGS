"""Per-run aggregation of per-sample TLA-seq outputs into final_report/.

Reads each sample's Results/<sample_id>/ artifacts and concatenates them into
run-level tables:

    final_report/vector_snp_summary.tsv
        per-sample VarScan SNV reports (all contigs concatenated)
    final_report/vector_indel_summary.tsv
        per-sample VarScan indel reports (all contigs concatenated)
    final_report/vector_sv_summary.tsv
        per-sample vector structural breakpoints (GRIDSS vector-vector)
    final_report/resolved_integration_sites_summary.tsv
        per-sample resolved integration sites
    final_report/coverage_plots/
        per-sample coverage PNGs copied and prefixed with sample ID

Each TSV has one global column header followed by `# Sample: <id>` banner
blocks.  Comment lines (starting with `#`) can be skipped in downstream
parsers with ``comment='#'``.
"""

import glob
import os
import shutil

import pandas as pd


def _read_tsv(path):
    """Return a DataFrame for path, or None if missing/empty."""
    if not os.path.exists(path):
        return None
    try:
        return pd.read_csv(path, sep='\t', comment='#')
    except pd.errors.EmptyDataError:
        return None


def _write_block_tsv(output_path, file_label, blocks, no_data_msg='# (no data)'):
    """Write a banner-block TSV from (sample_id, df_or_None) pairs.

    One global column header is written at the top; each sample gets a
    ``# Sample: <id>`` line followed by its rows (or no_data_msg when empty).
    """
    header_cols = None
    for _sid, df in blocks:
        if df is not None and not df.empty:
            header_cols = list(df.columns)
            break

    n_rows = 0
    with open(output_path, 'w') as fh:
        fh.write(f'# {file_label}\n')
        if header_cols is not None:
            fh.write('sample_id\t' + '\t'.join(header_cols) + '\n')
        for sample_id, df in blocks:
            fh.write(f'# Sample: {sample_id}\n')
            if df is not None and not df.empty and header_cols is not None:
                for _, r in df.iterrows():
                    vals = ['' if pd.isna(v) else str(v) for v in (r.get(c, '') for c in header_cols)]
                    fh.write(sample_id + '\t' + '\t'.join(vals) + '\n')
                    n_rows += 1
            else:
                fh.write(no_data_msg + '\n')
            fh.write('\n')
    print(f'[TLA Pipeline] {os.path.basename(output_path)}: {n_rows} rows -> {output_path}')


def aggregate_snv(results_dir, sample_ids, output_path):
    """Combine per-sample, per-contig vector SNV reports."""
    blocks = []
    for sid in sample_ids:
        pattern = os.path.join(
            results_dir, sid, 'vector_characterization', f'{sid}.*.snv_report.tsv')
        hits = sorted(glob.glob(pattern))
        if not hits:
            print(f'[TLA Pipeline] Warning: no SNV report for {sid}; skipping.')
            blocks.append((sid, None))
            continue
        dfs = [d for d in (_read_tsv(p) for p in hits) if d is not None and not d.empty]
        blocks.append((sid, pd.concat(dfs, ignore_index=True) if dfs else None))
    _write_block_tsv(output_path, 'vector_snp_summary: per-sample VarScan SNVs', blocks)


def aggregate_indel(results_dir, sample_ids, output_path):
    """Combine per-sample, per-contig vector indel reports."""
    blocks = []
    for sid in sample_ids:
        pattern = os.path.join(
            results_dir, sid, 'vector_characterization', f'{sid}.*.indel_report.tsv')
        hits = sorted(glob.glob(pattern))
        if not hits:
            print(f'[TLA Pipeline] Warning: no indel report for {sid}; skipping.')
            blocks.append((sid, None))
            continue
        dfs = [d for d in (_read_tsv(p) for p in hits) if d is not None and not d.empty]
        blocks.append((sid, pd.concat(dfs, ignore_index=True) if dfs else None))
    _write_block_tsv(output_path, 'vector_indel_summary: per-sample VarScan indels', blocks)


def aggregate_sv(results_dir, sample_ids, output_path):
    """Combine per-sample vector structural breakpoints."""
    blocks = []
    for sid in sample_ids:
        path = os.path.join(
            results_dir, sid, 'vector_characterization',
            f'{sid}.vector_structural_breakpoints.tsv')
        df = _read_tsv(path)
        if df is None:
            print(f'[TLA Pipeline] Warning: no vector SV file for {sid}; skipping.')
        blocks.append((sid, df))
    _write_block_tsv(
        output_path,
        'vector_sv_summary: per-sample GRIDSS vector-vector structural breakpoints',
        blocks,
        no_data_msg='No vector-vector SV detected.',
    )


def aggregate_resolved_sites(results_dir, sample_ids, output_path):
    """Combine per-sample resolved integration sites."""
    blocks = []
    for sid in sample_ids:
        path = os.path.join(results_dir, sid, 'integration', f'{sid}.resolved_sites.tsv')
        df = _read_tsv(path)
        if df is None:
            print(f'[TLA Pipeline] Warning: no resolved integration sites for {sid}; skipping.')
        blocks.append((sid, df))
    _write_block_tsv(
        output_path,
        'resolved_integration_sites_summary: per-sample resolved integration sites',
        blocks,
    )


def aggregate_integration_summary(results_dir, sample_ids, output_path):
    """Concatenate per-sample integration/<sid>.summary.tsv into one flat table."""
    dfs = []
    for sid in sample_ids:
        path = os.path.join(results_dir, sid, 'integration', f'{sid}.summary.tsv')
        df = _read_tsv(path)
        if df is None:
            print(f'[TLA Pipeline] Warning: no integration summary for {sid}; skipping.')
        else:
            dfs.append(df)
    if dfs:
        combined = pd.concat(dfs, ignore_index=True)
        tmp = output_path + '.tmp'
        combined.to_csv(tmp, sep='\t', index=False)
        os.replace(tmp, output_path)
        print(f'[TLA Pipeline] {os.path.basename(output_path)}: {len(combined)} row(s) -> {output_path}')
    else:
        print(f'[TLA Pipeline] Warning: no integration summaries found; {output_path} not written.')


def copy_coverage_plots(results_dir, sample_ids, plots_dir):
    """Copy per-sample PNG plots into plots_dir, prefixed with sample ID."""
    os.makedirs(plots_dir, exist_ok=True)
    n_copied = 0
    for sid in sample_ids:
        src_dirs = [
            os.path.join(results_dir, sid, 'coverage'),
            os.path.join(results_dir, sid, 'vector_characterization'),
        ]
        for src_dir in src_dirs:
            for png in sorted(glob.glob(os.path.join(src_dir, '*.png'))):
                dest = os.path.join(plots_dir, f'{sid}_{os.path.basename(png)}')
                if not os.path.exists(dest):
                    shutil.copy2(png, dest)
                n_copied += 1
    print(f'[TLA Pipeline] coverage_plots: {n_copied} plot(s) -> {plots_dir}')


def aggregate_run_report(results_dir, final_report_dir, sample_ids):
    """Write run-level aggregate reports under final_report_dir."""
    os.makedirs(final_report_dir, exist_ok=True)
    aggregate_snv(results_dir, sample_ids,
                  os.path.join(final_report_dir, 'vector_snp_summary.tsv'))
    aggregate_indel(results_dir, sample_ids,
                    os.path.join(final_report_dir, 'vector_indel_summary.tsv'))
    aggregate_sv(results_dir, sample_ids,
                 os.path.join(final_report_dir, 'vector_sv_summary.tsv'))
    aggregate_resolved_sites(results_dir, sample_ids,
                             os.path.join(final_report_dir,
                                          'resolved_integration_sites_summary.tsv'))
    aggregate_integration_summary(results_dir, sample_ids,
                                  os.path.join(final_report_dir,
                                               'integration_summary.tsv'))
    copy_coverage_plots(results_dir, sample_ids,
                        os.path.join(final_report_dir, 'coverage_plots'))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Aggregate TLA-seq run-level report from completed sample Results/.')
    parser.add_argument('-r', '--results_dir', required=True,
                        help='Path to Results/ (i.e. output_dir/Results)')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Destination directory for aggregate report files')
    args = parser.parse_args()

    discovered = sorted(
        d for d in os.listdir(args.results_dir)
        if os.path.isdir(os.path.join(args.results_dir, d))
    )
    if not discovered:
        print(f'No sample directories found in {args.results_dir}')
    else:
        print(f'[TLA Pipeline] Aggregating {len(discovered)} sample(s): {discovered}')
        aggregate_run_report(args.results_dir, args.output_dir, discovered)
