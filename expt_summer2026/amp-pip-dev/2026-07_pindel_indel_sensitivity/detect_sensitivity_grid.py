#!/usr/bin/env python3
"""Experiment A - indel detection sensitivity grid.

QUESTIONS ADDRESSED
  1. What is Pindel's detection limit as a function of indel SIZE and ALLELE
     FRACTION? i.e. for each indel length, what is the lowest allele fraction
     at which the event is still called exactly?
  2. How does that limit differ between DELETIONS and INSERTIONS?
  3. How does read length (75 / 150 / 250 bp) shift the limit, especially for
     large indels that must fit inside a single read?
  4. How does sequencing depth (5000 / 10000 / 100000 pairs) shift the limit?
  5. Is the limit reproducible across sequence contexts (tla vs LAGA_LC)?

Scoring is at the raw-detectability CEILING (min_support = config.CEILING_E = 1):
"can Pindel see the event at all", independent of any calling threshold. WT-only
controls (AF=0) at every depth/read-length give the false-positive baseline.

Reads simulated once here are also reused by Experiment C, which re-scores the
same base VCFs across a support-threshold sweep.
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

EXPERIMENT = "A"

# --- Experiment A grid (owned here, not in config) ---------------------------
# Allele fractions: 0.005 is intentionally excluded - not reliably detectable
# with short-read amplicon data.
ALLELE_FREQS = [0.01, 0.02, 0.05, 0.10, 0.20]
# Depth tiers. 100000 is kept for now; drop it at a later stage if runtime is a
# problem.
DEPTHS = [5000, 10000, 100000]


def main(datasets=None, jobs=None):
    if jobs is None:
        jobs = config.PARALLEL_JOBS
    align_call.startup_banner(jobs)
    if datasets is None:
        datasets = list(config.DATASETS)
    for dataset in datasets:
        # Per-dataset output tree: runs/sim_benchmark/exptA/<dataset>/
        dataset_outdir = os.path.join(config.OUTDIR, f"expt{EXPERIMENT}", dataset)

        wt_amplicon, wt_ref, refs_dir = simulate.ensure_wt_ref(
            dataset, dataset_outdir)
        amplicon_len = simulate.dataset_amplicon_len(dataset)
        insert_size = config.DATASETS[dataset]["pindel_insert"]

        # Fragment stats for standard paired-end simulation.
        ds_cfg = config.DATASETS[dataset]
        fragment_mean = ds_cfg["fragment_mean"]
        fragment_sd   = ds_cfg["fragment_sd"]

        # Positions: relative fractions mapped to 1-based amplicon coordinates.
        positions = config.relative_positions(amplicon_len)
        simulate.validate_positions(wt_amplicon, positions, f"{dataset} Exp A")

        # Skip indel lengths that are >= half the amplicon: injecting a 100 bp
        # indel into a 130 bp amplicon is not meaningful.
        active_lengths = [l for l in config.INDEL_LENGTHS
                          if l < amplicon_len // 2]
        # Skip read lengths longer than the amplicon or >= fragment_mean:
        # Pindel hard-rejects insert_size < read_len; ART requires
        # read_len < fragment_mean (strictly less than).
        # A per-dataset "read_lengths" key overrides the global READ_LENGTHS
        # (e.g. LAGA_LC omits r250 because pindel_insert < 250).
        _candidate_rls = ds_cfg.get("read_lengths", config.READ_LENGTHS)
        active_read_lengths = [rl for rl in _candidate_rls
                               if rl <= amplicon_len and rl < fragment_mean]
        print(f"[A] {dataset}: amplicon_len={amplicon_len}, "
              f"insert_size={insert_size}, positions={positions}, "
              f"indel_lengths={active_lengths}, "
              f"read_lengths={active_read_lengths}")

        # Build mutant references once per (indel_type, length, pos) tuple;
        # reuse across AF / depth / read_len / seed.
        mutant_cache = {}
        for indel_type, length, pos in itertools.product(
                config.INDEL_TYPES, active_lengths, positions):
            # Guard: a deletion must not extend beyond the amplicon end.
            if indel_type == "del" and pos + length - 1 > amplicon_len:
                continue
            mut_path, truth = simulate.build_mutant(
                wt_amplicon, dataset, indel_type, length, pos, refs_dir,
                dataset_outdir)
            mutant_cache[(indel_type, length, pos)] = (mut_path, truth)

        tasks = []
        for indel_type, length, pos, af, read_len, depth, seed in itertools.product(
                config.INDEL_TYPES, active_lengths, positions, ALLELE_FREQS,
                active_read_lengths, DEPTHS, config.SEEDS):
            if indel_type == "del" and pos + length - 1 > amplicon_len:
                continue
            mut_path, truth = mutant_cache[(indel_type, length, pos)]
            cond_id = config.condition_id(dataset, indel_type, length, af, pos,
                                          read_len, seed, depth)
            tasks.append({
                "is_control": False,
                "cond_id": cond_id,
                "mut_path": mut_path,
                "wt_ref": wt_ref,
                "af": af,
                "seed": seed,
                "depth": depth,
                "read_len": read_len,
                "indel_pos": pos,
                "indel_type": indel_type,
                "length": length,
                "dataset": dataset,
                "insert_size": insert_size,
                "outdir": dataset_outdir,
                "fragment_mean": fragment_mean,
                "fragment_sd": fragment_sd,
                "experiment": EXPERIMENT,
                "truth": truth,
            })

        # WT-only controls: one per depth x read_len x seed.
        for read_len, depth, seed in itertools.product(
                active_read_lengths, DEPTHS, config.SEEDS):
            cond_id = config.control_id(dataset, read_len, seed, depth)
            tasks.append({
                "is_control": True,
                "cond_id": cond_id,
                "wt_ref": wt_ref,
                "seed": seed,
                "depth": depth,
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
            print(f"[A] WARNING: {len(incomplete_tasks)} incomplete condition(s) "
                  f"-- skipping, not repairing:")
            for t in incomplete_tasks:
                print(f"  {t['cond_id']}")

        print(f"[A] {dataset}: {len(complete_tasks)} complete (skip), "
              f"{len(incomplete_tasks)} incomplete (warn+skip), "
              f"{len(pending_tasks)} pending (run)")

        print(f"[A] {dataset}: submitting {len(tasks)} conditions (jobs={jobs})")
        unique_refs = {t["wt_ref"] for t in tasks} | {
            t["mut_path"] for t in tasks if not t["is_control"]}
        align_call.pre_index_references(unique_refs)
        rows, failures = align_call.run_tasks_parallel(pending_tasks, jobs)
        # Emit rows for complete and incomplete (skipped) conditions too.
        rows += [align_call.row_from_task(t)
                 for t in complete_tasks + incomplete_tasks]
        if failures:
            print(f"\n[A] {len(failures)} condition(s) FAILED: {failures}")
        rows.sort(key=lambda r: r["condition_id"])

        simulate.write_run_table(rows, dataset_outdir, EXPERIMENT)
        evaluate.evaluate(dataset_outdir, EXPERIMENT, min_support=config.CEILING_E)
        print(f"\n[A] {dataset} done. "
              f"See runs/sim_benchmark/exptA/{dataset}/results/A/")

    print("\n[A] all datasets done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Experiment A: indel detection sensitivity grid.")
    parser.add_argument("-d", "--dataset", choices=list(config.DATASETS),
                        help="Run one dataset only (omit to run all).")
    args = parser.parse_args()
    main(
        datasets=[args.dataset] if args.dataset else None,
    )
