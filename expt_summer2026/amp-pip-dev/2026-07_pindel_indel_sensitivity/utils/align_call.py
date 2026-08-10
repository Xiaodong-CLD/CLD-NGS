#!/usr/bin/env python3
"""Per-condition pipeline tail: align -> call -> plot.

Given already-simulated reads for one condition, reuse the production alignment
pipeline UNCHANGED (bwa_runner with trimming disabled -> samtools_runner), call
indels with the sim-specific Pindel runner (sim_pindel, -e 1 base VCFs, no
deletion filter), and render a raw-depth coverage plot. The experiment scripts
build the condition grids and call process_condition / process_control here.
"""

import concurrent.futures
import glob
import os
import subprocess
import sys
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))            # utils/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config          # noqa: E402
import simulate        # noqa: E402
import sim_pindel      # noqa: E402
import coverage_plot   # noqa: E402

# Reuse the production alignment wrappers unchanged (from the repo root).
sys.path.insert(0, config.REPO_ROOT)
from bwa_wrapper import bwa_runner            # noqa: E402
from samtools_wrapper import samtools_runner  # noqa: E402


def condition_status(cond_dir):
    """Return 'complete', 'incomplete', or 'pending' for a condition directory.

    'complete'   – at least one base per-type VCF (size > 0) exists under
                   pindel-output/ AND a *.depth.txt file exists.
    'incomplete' – directory exists but those criteria are not met.
    'pending'    – directory does not exist yet.
    """
    if not os.path.isdir(cond_dir):
        return "pending"
    depth_files = glob.glob(os.path.join(cond_dir, "*.depth.txt"))
    has_depth = bool(depth_files)
    pindel_dir = os.path.join(cond_dir, "pindel-output")
    has_vcfs = (
        os.path.isdir(pindel_dir)
        and any(
            (n.endswith("_SI.vcf") or n.endswith("_D.vcf") or n.endswith("_LI.vcf"))
            and os.path.getsize(os.path.join(pindel_dir, n)) > 0
            for n in os.listdir(pindel_dir)
        )
    )
    return "complete" if (has_vcfs and has_depth) else "incomplete"


def align_call_plot(cond_id, cond_dir, r1, r2, wt_ref, indel_pos, read_len,
                    insert_size, outdir=None):
    """bwa -> samtools -> sim_pindel -> coverage plot for one condition.

    indel_pos may be None (WT-only controls) to omit the indel marker.
    insert_size is passed to Pindel (derived from config.dataset_pindel_insert).
    outdir: root benchmark output directory; defaults to config.OUTDIR.
    """
    effective_outdir = outdir if outdir is not None else config.OUTDIR
    bam = bwa_runner(wt_ref, r1, r2, disable_trimming=True)
    if not bam:
        raise RuntimeError(f"{cond_id}: bwa alignment failed")
    _, sorted_bam, _ = samtools_runner(bam, wt_ref)
    if not sorted_bam:
        raise RuntimeError(f"{cond_id}: samtools processing failed")
    if insert_size < read_len:
        print(f"[align_call] SKIP Pindel for {cond_id}: "
              f"insert_size={insert_size} < read_len={read_len}; "
              f"Pindel requires insert_size >= read_len")
    else:
        sim_pindel.run_sim_pindel(sorted_bam, wt_ref, insert_size)

    depth_file = bam.replace(".bam", ".depth.txt")
    if os.path.exists(depth_file):
        cov_dir = os.path.join(effective_outdir, "results", "coverage")
        os.makedirs(cov_dir, exist_ok=True)
        coverage_plot.plot_depth(
            depth_file,
            os.path.join(cov_dir, f"{cond_id}.coverage.png"),
            title=cond_id,
            indel_pos=indel_pos,
            read_len=read_len,
        )
    else:
        print(f"[align_call] WARNING: depth file not found, skipping plot: {depth_file}")


def process_condition(cond_id, mutant_ref, wt_ref, af, seed, depth, read_len,
                      indel_pos, insert_size, force=False, outdir=None,
                      fragment_mean=None, fragment_sd=None):
    effective_outdir = outdir if outdir is not None else config.OUTDIR
    cond_dir = os.path.join(effective_outdir, "conditions", cond_id)
    os.makedirs(cond_dir, exist_ok=True)
    if condition_status(cond_dir) == "complete" and not force:
        print(f"[align_call] skip {cond_id}: already complete")
        return cond_dir
    r1, r2 = simulate.simulate_condition(cond_dir, mutant_ref, wt_ref, af, seed,
                                         depth, read_len, fragment_mean, fragment_sd)
    align_call_plot(cond_id, cond_dir, r1, r2, wt_ref, indel_pos, read_len,
                    insert_size, effective_outdir)
    return cond_dir


def process_control(cond_id, wt_ref, seed, depth, read_len, insert_size,
                    force=False, outdir=None, fragment_mean=None, fragment_sd=None):
    effective_outdir = outdir if outdir is not None else config.OUTDIR
    cond_dir = os.path.join(effective_outdir, "conditions", cond_id)
    os.makedirs(cond_dir, exist_ok=True)
    if condition_status(cond_dir) == "complete" and not force:
        print(f"[align_call] skip {cond_id}: already complete")
        return cond_dir
    r1, r2 = simulate.simulate_control(cond_dir, wt_ref, seed, depth, read_len,
                                       fragment_mean, fragment_sd)
    align_call_plot(cond_id, cond_dir, r1, r2, wt_ref, None, read_len,
                    insert_size, effective_outdir)
    return cond_dir


# ---------------------------------------------------------------------------
# Parallel execution infrastructure
# ---------------------------------------------------------------------------

def row_from_task(task):
    """Assemble a manifest row from a task dict without running anything."""
    if task["is_control"]:
        return simulate.control_row(
            task["cond_id"], task["experiment"], task["dataset"],
            task["read_len"], task["depth"], task["seed"],
        )
    return simulate.condition_row(
        task["cond_id"], task["experiment"], task["truth"],
        task["af"], task["read_len"], task["depth"], task["seed"],
    )


def _run_task(task):
    """Worker: run one condition and return its manifest row.

    Top-level so it is picklable by ProcessPoolExecutor.
    Sets OMP_NUM_THREADS=1 to pin the aligner's OpenMP thread pool;
    per-condition parallelism is at the process level only.
    """
    os.environ["OMP_NUM_THREADS"] = "1"
    cond_id = task["cond_id"]
    if task["is_control"]:
        process_control(
            cond_id, task["wt_ref"], task["seed"], task["depth"],
            task["read_len"], task["insert_size"],
            outdir=task["outdir"],
            fragment_mean=task["fragment_mean"],
            fragment_sd=task["fragment_sd"],
        )
        return simulate.control_row(
            cond_id, task["experiment"], task["dataset"],
            task["read_len"], task["depth"], task["seed"],
        )
    process_condition(
        cond_id, task["mut_path"], task["wt_ref"],
        task["af"], task["seed"], task["depth"], task["read_len"],
        task["indel_pos"], task["insert_size"],
        outdir=task["outdir"],
        fragment_mean=task["fragment_mean"],
        fragment_sd=task["fragment_sd"],
    )
    return simulate.condition_row(
        cond_id, task["experiment"], task["truth"],
        task["af"], task["read_len"], task["depth"], task["seed"],
    )


def pre_index_references(refs):
    """Index a set of reference FASTAs sequentially before parallel dispatch.

    bwa_runner always calls 'bwa index' unconditionally, creating a write-write
    race when multiple workers share the same reference file.  Indexing here
    first ensures the index is complete before any worker starts; each worker's
    redundant re-index then only overwrites identical data.  For a ~400 bp
    amplicon the re-index window is <1 ms.
    """
    for ref in sorted(refs):
        if not os.path.exists(ref + ".bwt"):
            print(f"[align_call] bwa index {ref}")
            subprocess.run(["bwa", "index", ref], check=True)


def run_tasks_parallel(tasks, jobs=None):
    """Submit tasks to a ProcessPoolExecutor; return (rows, failed_cond_ids).

    jobs=1 uses the same ProcessPoolExecutor code path as jobs=N, so there
    is only one execution path to debug.
    Rows are returned unsorted; callers sort by condition_id for a stable,
    run-order-independent manifest.
    Failures are logged but do not abort the remaining conditions.
    """
    if jobs is None:
        jobs = config.PARALLEL_JOBS
    rows, failures = [], []
    with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as pool:
        fut_to_id = {pool.submit(_run_task, t): t["cond_id"] for t in tasks}
        for fut in concurrent.futures.as_completed(fut_to_id):
            cond_id = fut_to_id[fut]
            try:
                rows.append(fut.result())
            except Exception as exc:
                print(f"[align_call] ERROR {cond_id}: {exc}")
                traceback.print_exc()
                failures.append(cond_id)
    return rows, failures


def startup_banner(jobs):
    """Print jobs / cpu_count / total RAM at experiment startup."""
    ram = "unknown"
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemTotal:"):
                    ram = f"{int(line.split()[1]) / 1e6:.1f} GB"
                    break
    except OSError:
        pass
    print(f"[parallel] jobs={jobs}  cpu_count={os.cpu_count()}  total_RAM={ram}")
    print("[parallel] NOTE: Pindel is memory-hungry; "
          "practical ceiling is often RAM, not cores.")
