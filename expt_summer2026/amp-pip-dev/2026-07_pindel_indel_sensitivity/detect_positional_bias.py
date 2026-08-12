#!/usr/bin/env python3
"""Experiment B - positional bias of indel detection along the amplicon.

QUESTIONS ADDRESSED
  1. Does Pindel's ability to call an indel depend on WHERE in the amplicon the
     breakpoint sits - near the primers (well covered by read ends) vs the
     amplicon center (the insert-gap blind spot between paired reads)?
  2. Is that positional pattern different for deletions vs insertions?
  3. Does a longer read length (which shrinks the center gap) remove the
     positional penalty?

A single moderately-detectable event (length = 10 bp, allele fraction = 0.10) is
walked across breakpoint positions from primer-proximal to center to the far
primer, holding everything else fixed, so any change in detection is attributable
to position alone. Scoring is at the raw-detectability ceiling
(min_support = config.CEILING_E = 1); WT-only controls give the FP baseline.
"""

import argparse
import itertools
import os
import sys

BENCH_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH_DIR)
sys.path.insert(0, os.path.join(BENCH_DIR, "utils"))

import config       # noqa: E402
import simulate     # noqa: E402
import align_call   # noqa: E402
import evaluate     # noqa: E402

EXPERIMENT = "B"

# --- Experiment B grid (owned here, not in config) ---------------------------
# Positions are computed relative to each dataset's amplicon length at runtime
# (see make_positions), so this experiment works for all amplicon sizes.
AF = 0.10
LENGTH = 10
DEPTH = config.DEFAULT_DEPTH


def make_positions(amplicon_len, read_len=75):
    """7 representative breakpoints spanning the amplicon for a given read length.

    Covers: mid-5'-read-zone, just inside 5' gap boundary, just outside it,
    dead center, just outside 3' gap boundary, just inside it, mid-3'-read-zone.
    For short amplicons where the read fully spans (read_len >= amplicon_len//2),
    positions collapse to primer-proximal and center only.
    """
    gap_start = read_len + 5          # first position a 75 bp R1 cannot cover
    gap_end   = amplicon_len - read_len - 5   # last position R2 cannot cover
    center    = amplicon_len // 2
    pts = [
        max(2, read_len // 2),                 # mid 5' read zone
        max(2, gap_start - 10),                # just inside 5' coverage
        min(amplicon_len - 2, gap_start + 10), # just outside (blind)
        center,
        max(2, gap_end - 10),                  # just outside (blind, 3' side)
        min(amplicon_len - 2, gap_end + 10),   # just inside 3' coverage
        min(amplicon_len - 2, amplicon_len - read_len // 2),  # mid 3' zone
    ]
    # Deduplicate and sort; clamp to valid 1-based coords.
    pts = sorted(set(max(2, min(amplicon_len - 2, p)) for p in pts))
    return pts


def main(datasets=None, jobs=None):
    if jobs is None:
        jobs = config.PARALLEL_JOBS
    align_call.startup_banner(jobs)
    if datasets is None:
        datasets = list(config.DATASETS)
    for dataset in datasets:
        # Per-dataset output tree: runs/sim_benchmark/exptB/<dataset>/
        dataset_outdir = os.path.join(config.OUTDIR, f"expt{EXPERIMENT}", dataset)

        wt_amplicon, wt_ref, refs_dir = simulate.ensure_wt_ref(
            dataset, dataset_outdir)
        amplicon_len = simulate.dataset_amplicon_len(dataset)
        insert_size  = config.DATASETS[dataset]["pindel_insert"]

        ds_cfg = config.DATASETS[dataset]
        fragment_mean = ds_cfg["fragment_mean"]
        fragment_sd   = ds_cfg["fragment_sd"]

        # Positions depend on the shortest read length so they stay meaningful
        # across all swept read lengths.
        # A per-dataset "read_lengths" key overrides the global READ_LENGTHS
        # (e.g. LAGA_LC omits r250 because pindel_insert < 250).
        _candidate_rls = ds_cfg.get("read_lengths", config.READ_LENGTHS)
        positions = make_positions(amplicon_len, read_len=min(_candidate_rls))
        # Skip read lengths longer than the amplicon or >= fragment_mean:
        # Pindel hard-rejects insert_size < read_len; ART requires
        # read_len < fragment_mean (strictly less than).
        active_read_lengths = [rl for rl in _candidate_rls
                               if rl <= amplicon_len and rl < fragment_mean]
        print(f"[B] {dataset}: amplicon_len={amplicon_len}, "
              f"insert_size={insert_size}, positions={positions}, "
              f"read_lengths={active_read_lengths}")
        simulate.validate_positions(wt_amplicon, positions, f"{dataset} Exp B")

        tasks = []
        for indel_type, pos in itertools.product(config.INDEL_TYPES, positions):
            mut_path, truth = simulate.build_mutant(
                wt_amplicon, dataset, indel_type, LENGTH, pos, refs_dir,
                dataset_outdir)
            for read_len, seed in itertools.product(active_read_lengths,
                                                    config.SEEDS):
                cond_id = config.condition_id(dataset, indel_type, LENGTH, AF,
                                              pos, read_len, seed, DEPTH)
                tasks.append({
                    "is_control": False,
                    "cond_id": cond_id,
                    "mut_path": mut_path,
                    "wt_ref": wt_ref,
                    "af": AF,
                    "seed": seed,
                    "depth": DEPTH,
                    "read_len": read_len,
                    "indel_pos": pos,
                    "indel_type": indel_type,
                    "length": LENGTH,
                    "dataset": dataset,
                    "insert_size": insert_size,
                    "outdir": dataset_outdir,
                    "fragment_mean": fragment_mean,
                    "fragment_sd": fragment_sd,
                    "experiment": EXPERIMENT,
                    "truth": truth,
                })

        # WT-only controls: one per read_len x seed at the fixed depth.
        for read_len, seed in itertools.product(active_read_lengths, config.SEEDS):
            cond_id = config.control_id(dataset, read_len, seed, DEPTH)
            tasks.append({
                "is_control": True,
                "cond_id": cond_id,
                "wt_ref": wt_ref,
                "seed": seed,
                "depth": DEPTH,
                "read_len": read_len,
                "insert_size": insert_size,
                "outdir": dataset_outdir,
                "fragment_mean": fragment_mean,
                "fragment_sd": fragment_sd,
                "experiment": EXPERIMENT,
                "dataset": dataset,
            })

        # --- pre-flight: categorise complete / incomplete / pending ---
        complete_tasks, incomplete_tasks, pending_tasks = [], [], []
        for t in tasks:
            cond_dir = os.path.join(dataset_outdir, "conditions", t["cond_id"])
            status = align_call.condition_status(cond_dir)
            if status == "complete":
                complete_tasks.append(t)
            elif status == "incomplete":
                incomplete_tasks.append(t)
            else:
                pending_tasks.append(t)

        if incomplete_tasks:
            print(f"[B] WARNING: {len(incomplete_tasks)} incomplete condition(s) "
                  f"-- skipping, not repairing:")
            for t in incomplete_tasks:
                print(f"  {t['cond_id']}")

        print(f"[B] {dataset}: {len(complete_tasks)} complete (skip), "
              f"{len(incomplete_tasks)} incomplete (warn+skip), "
              f"{len(pending_tasks)} pending (run)")

        print(f"[B] {dataset}: submitting {len(pending_tasks)} conditions (jobs={jobs})")
        unique_refs = {t["wt_ref"] for t in tasks} | {
            t["mut_path"] for t in tasks if not t["is_control"]}
        align_call.pre_index_references(unique_refs)
        rows, failures = align_call.run_tasks_parallel(pending_tasks, jobs)
        rows += [align_call.row_from_task(t)
                 for t in complete_tasks + incomplete_tasks]
        if failures:
            print(f"\n[B] {len(failures)} condition(s) FAILED: {failures}")
        rows.sort(key=lambda r: r["condition_id"])

        simulate.write_run_table(rows, dataset_outdir, EXPERIMENT)
        evaluate.evaluate(dataset_outdir, EXPERIMENT, min_support=config.CEILING_E)
        print(f"\n[B] {dataset} done. "
              f"See runs/sim_benchmark/exptB/{dataset}/results/B/")

    print("\n[B] all datasets done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Experiment B: positional bias of indel detection.")
    parser.add_argument("-d", "--dataset", choices=list(config.DATASETS),
                        help="Run one dataset only (omit to run all).")
    args = parser.parse_args()
    main(
        datasets=[args.dataset] if args.dataset else None,
    )
