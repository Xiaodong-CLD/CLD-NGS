#!/usr/bin/env python3
"""Guard tests for the two easy-to-break invariants of the benchmark:

  1. simulate.make_mutant_sequence coordinate/anchor convention (the VCF-style
     anchor is the unchanged base 5' of the event, and the mutant sequence is
     exactly WT with the event applied).
  2. simulate.mix_counts allele-fraction mixing math (mut = round(af*depth),
     wt + mut == depth).

Run: python test_benchmark.py   (also works under pytest)
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "utils"))

import config
import simulate


WT = "ACGTACGTGGCCAATTACGTACGTTTGGCCAAACGT"  # arbitrary, homopolymer-free-ish

# Local test grids (experiment grids no longer live in config).
TEST_ALLELE_FREQS = [0.01, 0.02, 0.05, 0.10, 0.20]
TEST_DEPTHS = [5000, 10000, 100000]


def test_deletion_anchor_convention():
    pos, length = 10, 3
    mutant, ref_allele, alt_allele, anchor = simulate.make_mutant_sequence(
        WT, "del", length, pos
    )
    idx = pos - 1
    assert anchor == pos - 1                      # anchor is base before event
    assert mutant == WT[:idx] + WT[idx + length:]  # exactly the deletion
    assert alt_allele == WT[anchor - 1]           # anchor base only
    assert ref_allele == WT[anchor - 1] + WT[idx:idx + length]
    assert len(WT) - len(mutant) == length


def test_insertion_anchor_convention():
    pos, length = 12, 4
    mutant, ref_allele, alt_allele, anchor = simulate.make_mutant_sequence(
        WT, "ins", length, pos
    )
    idx = pos - 1
    inserted = config.INSERTION_TEMPLATE[:length]
    assert anchor == pos - 1
    assert mutant == WT[:idx] + inserted + WT[idx:]
    assert ref_allele == WT[anchor - 1]
    assert alt_allele == WT[anchor - 1] + inserted
    assert len(mutant) - len(WT) == length


def test_mix_counts_sums_to_depth():
    for af in TEST_ALLELE_FREQS:
        for depth in TEST_DEPTHS:
            wt_pairs, mut_pairs = simulate.mix_counts(af, depth)
            assert mut_pairs == round(af * depth)
            assert wt_pairs + mut_pairs == depth
            assert mut_pairs > 0


def test_mix_counts_rejects_zero_mutant():
    try:
        simulate.mix_counts(0.0, config.DEFAULT_DEPTH)
    except ValueError:
        return
    raise AssertionError("mix_counts(0.0, depth) should raise ValueError")


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"[ok] {t.__name__}")
    print(f"[done] {len(tests)} tests passed")
