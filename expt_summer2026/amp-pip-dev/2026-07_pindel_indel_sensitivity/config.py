#!/usr/bin/env python3
"""Shared configuration for the Pindel indel-detection sensitivity benchmark.

Holds only constants shared across all experiments (paths, amplicon/dataset
geometry, sequencing primitives, indel building blocks, matching tolerances,
the insertion template, and the condition/control id helpers). Experiment-
specific grids (allele fractions, depths, breakpoint positions, support-
threshold sweeps) live in the individual experiment scripts, not here.

Used by utils/simulate.py, utils/sim_pindel.py, utils/align_call.py,
utils/evaluate.py and the three experiment scripts.
"""

import os
import random

# --- repository layout -------------------------------------------------------
_EXPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(_EXPT_DIR))), "amp-pip-dev"
)

# Where all simulation artefacts go (inside this experiment folder).
OUTDIR = os.path.join(_EXPT_DIR, "runs", "sim_benchmark")

# --- BAM for LAGA-LC fragment statistics -------------------------------------
# Set to an aligned BAM for the LAGA-LC amplicon to auto-derive fragment_mean/sd
# from the TLEN distribution. Leave None to use the defaults.
LAGA_LC_BAM = None


def fragment_stats_from_bam(bam_path):
    """Estimate (fragment_mean, fragment_sd) from positive TLEN values in a BAM.

    Uses samtools view. Filters TLEN to (0, 5000] to exclude unmapped/chimeric
    reads. Raises ValueError if fewer than 100 usable TLEN values are found.
    """
    import statistics
    import subprocess
    result = subprocess.run(
        ["samtools", "view", bam_path],
        capture_output=True, text=True, check=True,
    )
    tlens = []
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) > 8:
            try:
                t = int(fields[8])
                if 0 < t <= 5000:
                    tlens.append(t)
            except ValueError:
                pass
    if len(tlens) < 100:
        raise ValueError(
            f"Too few positive TLEN values in {bam_path}: {len(tlens)}"
        )
    return round(statistics.mean(tlens)), round(statistics.stdev(tlens))


_laga_lc_frag_mean, _laga_lc_frag_sd = 350, 50
if LAGA_LC_BAM and os.path.exists(LAGA_LC_BAM):
    try:
        _laga_lc_frag_mean, _laga_lc_frag_sd = fragment_stats_from_bam(LAGA_LC_BAM)
        print(f"[config] LAGA_LC fragment stats from BAM: "
              f"mean={_laga_lc_frag_mean}, sd={_laga_lc_frag_sd}")
    except Exception as _e:
        print(f"[config] WARNING: could not read LAGA_LC fragment stats: {_e}")

# --- amplicon datasets -------------------------------------------------------
# Each entry has the same four keys:
#   fasta         : absolute path to the reference FASTA.
#   region        : (offset, length) to carve from the FASTA (0-based), or None
#                   to use the full sequence.
#   fragment_mean : mean fragment length (bp) for ART -m.
#   fragment_sd   : fragment-length standard deviation for ART -s.
#   pindel_insert : insert size passed to Pindel -i (config column 2).
# NOTE: pindel_insert is deliberately DECOUPLED from fragment_mean. Setting
# Pindel -i to the true fragment mean suppresses split-read detection on these
# amplicons (each amplicon has an empirical workable -i window); LAGA_LC in
# particular yields zero calls at -i=fragment_mean (350). See benchmark notes.
DATASETS = {
    # Synthetic 400 bp window from the TLA D08 reference (benchmark baseline).
    "tla": {
        "fasta":         os.path.join(_EXPT_DIR, "references", "D08_HI_HS_AvB1.fa"),
        "region":        (3000, 400),
        "fragment_mean": 250,
        "fragment_sd":   30,
        "pindel_insert": 250,
    },
    # Real LC amplicon (SAD0309 run, 711 bp); primary benchmark dataset.
    "LAGA_LC": {
        "fasta":         os.path.join(_EXPT_DIR, "references", "SAD0309_P01_G04_LAGA-LC.fa"),
        "region":        None,
        "fragment_mean": _laga_lc_frag_mean,
        "fragment_sd":   _laga_lc_frag_sd,
        # Empirically the true variant is recovered only for -i ~155-175 (peak
        # ~165); fragment_mean (350) yields zero calls. Note: 165 < read_len 250,
        # so r250 is excluded entirely (Pindel requires insert_size >= read_len).
        "pindel_insert": 165,
        # r250 is structurally unevaluable: pindel_insert (165) < read_len (250).
        # Exclude it so no conditions, results rows, or output directories are
        # generated for this tier. r75 and r150 are unaffected.
        "read_lengths":  [75, 150],
    },
}


def amplicon_name(dataset):
    """FASTA record / file stem for a dataset's amplicon."""
    return f"amplicon_{dataset}"


# --- indel placement ---------------------------------------------------------
POSITION_FRACS = [0.10, 0.25, 0.50, 0.75, 0.90]


def relative_positions(amplicon_len, fracs=None, margin=2):
    """Convert relative fractions to 1-based amplicon coordinates.

    Positions are rounded to integers and clamped to [margin, amplicon_len - margin]
    to ensure a valid anchor base exists on both sides. Duplicates are removed.
    """
    if fracs is None:
        fracs = POSITION_FRACS
    positions = []
    for frac in fracs:
        pos = max(margin, min(amplicon_len - margin, round(frac * amplicon_len)))
        if pos not in positions:
            positions.append(pos)
    return sorted(positions)


# --- sequencing simulation ---------------------------------------------------
READ_LENGTHS = [75, 150, 250]
READ_LEN = 75
SEQ_SYS = "MSv3"       # ART built-in profile (MiSeq v3); truncated to read len
DEFAULT_DEPTH = 5000
PARALLEL_JOBS = 4   # ProcessPoolExecutor workers; BWA uses -t 16 per job

# --- indel building blocks (shared across experiments) -----------------------
INDEL_TYPES = ["del", "ins"]
INDEL_LENGTHS = [1, 5, 10, 20, 50, 100]   # bp
SEEDS = [5, 12, 42, 0, 1, 7, 13, 21, 99, 123]

# --- evaluation tolerances ---------------------------------------------------
POS_TOLERANCE = 0         # bp; exact breakpoint required
SIZE_TOLERANCE = 0        # bp; exact size required
SIZE_TOLERANCE_FRAC = 0.0
NEARMISS_POS_WINDOW = 15  # diagnostic only; near-miss window for manual inspection

# Supporting-read floor for the raw-detectability ceiling.
CEILING_E = 1


def _make_insertion_template(length=200, seed=1234):
    """Deterministic, non-repetitive sequence used as inserted bases."""
    rng = random.Random(seed)
    bases = "ACGT"
    seq = []
    while len(seq) < length:
        b = rng.choice(bases)
        if len(seq) >= 2 and seq[-1] == b and seq[-2] == b:
            continue
        seq.append(b)
    return "".join(seq)


INSERTION_TEMPLATE = _make_insertion_template(max(INDEL_LENGTHS))


def condition_id(dataset, indel_type, length, af, pos, read_len, seed, depth):
    """Stable, filesystem-safe identifier for one benchmark condition."""
    af_tag = f"{af * 100:g}".replace(".", "p")
    return (f"{dataset}_{indel_type}_L{length}_AF{af_tag}_p{pos}"
            f"_r{read_len}_d{depth}_s{seed}")


def control_id(dataset, read_len, seed, depth):
    """Identifier for an AF=0 WT-only control condition."""
    return f"{dataset}_control_r{read_len}_d{depth}_s{seed}"
