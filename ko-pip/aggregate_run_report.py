"""Per-run aggregation of per-sample outputs into final_report/.

Reads each sample's persistent Results/<sample>/ artifacts and concatenates
them into run-level tables:

    final_report/all_snv.tsv            VarScan SNVs (var_freq > 1%), clone-annotated
    final_report/all_indel.tsv          merged GATK+VarScan indels, clone-annotated
    final_report/all_uninformative.tsv  per-sample uninformative-read summaries
    final_report/all_breakpoints.tsv    flagged breakpoint clusters (construct
                                        match / insert like), clone-annotated

`sample_records` is a list of (folder_name, sample_id, clone_id) tuples; the
folder_name is the Results/<folder_name>/ directory for that sample.
"""

import os

import pandas as pd

MIN_SNV_VAR_FREQ = 0.01  # keep SNVs with VarScan var_freq strictly above 1%
BREAKPOINT_FLAGS = ('construct match', 'insert like')  # flags to aggregate


def _read_tsv(path):
    """Read a TSV into a DataFrame, or return None if missing/empty."""
    if not os.path.exists(path):
        return None
    try:
        return pd.read_csv(path, sep='\t')
    except pd.errors.EmptyDataError:
        return None


def _collect(results_dir, sample_records, rel_path):
    """Read rel_path under each sample folder, tag with clone_id, concat."""
    frames = []
    for folder_name, _sample_id, clone_id in sample_records:
        df = _read_tsv(os.path.join(results_dir, folder_name, rel_path))
        if df is None or df.empty:
            continue
        df = df.copy()
        df.insert(0, 'clone_id', clone_id)
        df.insert(0, 'sample', folder_name)
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def aggregate_snv(results_dir, sample_records, output_path):
    """Concatenate per-sample SNV reports, keep var_freq > 1%, sort by freq."""
    df = _collect(results_dir, sample_records, os.path.join('snv', 'varscan_snp_report.tsv'))
    if not df.empty and 'var_freq' in df.columns:
        df = df[df['var_freq'].astype(float) > MIN_SNV_VAR_FREQ]
        df = df.sort_values('var_freq', ascending=False)
    df.to_csv(output_path, sep='\t', index=False)
    print(f'[Takeda Pipeline Message] all_snv wrote {len(df)} rows -> {output_path}')


def aggregate_indel(results_dir, sample_records, output_path):
    """Concatenate per-sample merged-indel tables."""
    df = _collect(results_dir, sample_records, os.path.join('indel', 'merged_indel.tsv'))
    df.to_csv(output_path, sep='\t', index=False)
    print(f'[Takeda Pipeline Message] all_indel wrote {len(df)} rows -> {output_path}')


def aggregate_uninformative(results_dir, sample_records, output_path):
    """Concatenate per-sample uninformative-read summaries."""
    df = _collect(results_dir, sample_records,
                  os.path.join('uninformative_analysis', 'uninformative_summary.tsv'))
    df.to_csv(output_path, sep='\t', index=False)
    print(f'[Takeda Pipeline Message] all_uninformative wrote {len(df)} rows -> {output_path}')


def aggregate_breakpoints(results_dir, sample_records, output_path):
    """Concatenate flagged breakpoint clusters (construct match / insert like).

    Reads breakpoint_analysis/top_clusters.tsv per sample (which carries '#'
    comment metadata lines above the header) and keeps only rows whose ``flag``
    column is a real breakpoint flag -- unflagged background clusters are
    dropped so the run-level table holds only candidate events.
    """
    frames = []
    for folder_name, _sample_id, clone_id in sample_records:
        path = os.path.join(results_dir, folder_name,
                            'breakpoint_analysis', 'top_clusters.tsv')
        if not os.path.exists(path):
            continue
        try:
            df = pd.read_csv(path, sep='\t', comment='#')
        except pd.errors.EmptyDataError:
            continue
        if df.empty or 'flag' not in df.columns:
            continue
        df = df[df['flag'].isin(BREAKPOINT_FLAGS)].copy()
        if df.empty:
            continue
        df.insert(0, 'clone_id', clone_id)
        df.insert(0, 'sample', folder_name)
        frames.append(df)
    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    out.to_csv(output_path, sep='\t', index=False)
    print(f'[Takeda Pipeline Message] all_breakpoints wrote {len(out)} rows -> {output_path}')


def aggregate_run_report(results_dir, final_report_dir, sample_records):
    """Write the run-level tables under final_report_dir."""
    os.makedirs(final_report_dir, exist_ok=True)
    aggregate_snv(results_dir, sample_records, os.path.join(final_report_dir, 'all_snv.tsv'))
    aggregate_indel(results_dir, sample_records, os.path.join(final_report_dir, 'all_indel.tsv'))
    aggregate_uninformative(results_dir, sample_records,
                            os.path.join(final_report_dir, 'all_uninformative.tsv'))
    aggregate_breakpoints(results_dir, sample_records,
                          os.path.join(final_report_dir, 'all_breakpoints.tsv'))


if __name__ == '__main__':
    pass
