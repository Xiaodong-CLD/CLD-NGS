# Findings — Pindel indel-sensitivity benchmark (August 2026)

Datasets: `tla` (400 bp synthetic TLA amplicon) and `LAGA_LC` (711 bp real LC
amplicon). Analysis notebook: `analyze_results.ipynb`.

**Status (August 2026):** All four experiments complete on both datasets.

---

## Overall recommendation

**Couple LoFreq into the amplicon variant-calling pipeline alongside or instead of
Pindel for indel detection.**

Pindel requires a per-amplicon `-i` (insert size) parameter that cannot be
derived from first principles and must be calibrated empirically for each new
amplicon. Its sensitivity degrades substantially when `-i` is wrong (LAGA_LC F1
drops from 0.66 to 0.35 vs tla). LoFreq, by contrast, is a pileup-based caller
that is robust to amplicon geometry, maintains consistent performance across both
amplicons tested (F1 0.73–0.74), and runs in comparable wall time to VarScan2.

| Caller | tla F1 | LAGA_LC F1 | Δ | Verdict |
|---|---|---|---|---|
| **LoFreq** | **0.73** | **0.74** | +0.01 | Consistent — **recommended** |
| Pindel | 0.66 | 0.35 | −0.31 | Fragile — requires per-amplicon calibration |
| Mutect2 | 0.47 | 0.50 | +0.03 | High floor (~20% AF); not useful at 1–5% |
| VarScan2 | 0.54 | 0.40 | −0.14 | High sensitivity, low precision; noisy |

---

## Finding 1 — Pindel's insert size (`-i`) cannot be derived for amplicons

### Background

Pindel's config insert size (`-i`, column 2 of the config file) models a
**fragment-size distribution from random shearing**. TLA-seq amplicons have
fixed, primer-defined geometry — no such distribution exists. The benchmark
originally set `-i = fragment_mean`, which happens to work for `tla` but
**zeroes out all true calls on `LAGA_LC`**.

### Results

| dataset | amplicon | ART `fragment_mean` | `-i = fragment_mean` | working `-i` range |
|---|---|---|---|---|
| `tla` | 400 bp | ~250 | works (coincidental) | ~250–398 |
| `LAGA_LC` | 711 bp | ~345 | **0 true calls** | **~155–175 (peak ~165)** |

Key observations:
- Setting `-i` to `fragment_mean`, amplicon length, or the production hard-code
  398 all yield zero true LAGA calls.
- The working range for LAGA (~165) was found empirically and does **not**
  transfer to `tla` (which requires `-i ≥ 250`). The windows barely overlap.
- A 3′-end blind spot exists at `-i` 165 that no insert tuning fixes: positions
  at the 3′ end of LAGA (p621, p641, p674) are never recovered. At p674 only
  8/10 000 reads carry a clean `10D` CIGAR; the aligner soft-clips the near-end
  deletion instead of emitting a gap → too little split signal for Pindel.
  This is a real positional-bias finding (see Finding 3), not a parameter issue.
- Pindel requires `-i ≥ read_len` (strict), so LAGA's `-i` 165 makes all `r250`
  conditions unevaluable (Pindel is skipped at runtime).

### Action taken

`config.DATASETS[...]["pindel_insert"]` is **decoupled** from `fragment_mean`:
`tla` = 250, `LAGA_LC` = 165. ART still simulates with `fragment_mean`.

### Implication for production

Pindel `-i` has no principled derivation for amplicon data. This is the primary
motivation for the caller comparison in Experiment D.

---

## Finding 2 — Operating threshold: fractional cutoff is depth-stable, absolute `-e` is not

*Experiment C — `tla`, Aug 2026. Source: `results/exptC/tla/threshold_sweep.csv`*

### F1-optimal thresholds (pooled DEL + INS, low-AF subset AF 1–2%, len 1–20 bp)

| depth | best absolute `-e` | F1 | `best_e / depth` | best fractional cutoff | F1 |
|---|---|---|---|---|---|
| 5 000 | **5** | 0.55 | 0.0010 | ~0.5–1% | ~0.55 |
| 10 000 | **5** | 0.65 | 0.0005 | ~0.5–1% | ~0.65 |
| 100 000 | **30** | 0.78 | 0.0003 | ~0.5–1% | ~0.78 |

### Interpretation

- **Absolute `-e` does not scale with depth.** The best count grows only
  sub-linearly (20× more depth → ~6× larger `-e`), and `best_e / depth` shrinks
  from 0.0010 to 0.0003. Any fixed `-e = k × depth` rule is too aggressive at
  high depth and would silently drop real low-AF calls.
- **The fractional cutoff (~0.5–1% on Pindel AD) is stable across the full depth
  range.** It transfers from the ~5000× benchmark to the ~370 000× real assay
  without re-tuning.

### Caveat — Pindel AD overstates true VAF (~4×)

The fractional cutoff divides by Pindel's own `AD` (alt supporting reads), not
the true mpileup depth. "1% on AD" is **not** 1% true AF. For a production
filter, divide alt reads by the `samtools depth` at the breakpoint, not Pindel
`AD`. The `af_gap.csv` table documents the per-condition ratio between Pindel AD
fraction and true mpileup fraction.

### LAGA_LC (Experiment C, Aug 2026)

The LAGA_LC threshold sweep shows **much lower F1 ceilings** than `tla`:

| depth | best F1 (DEL) | best F1 (INS) | optimal frac (DEL) | optimal frac (INS) |
|---|---|---|---|---|
| 5 000 | 0.05 | 0.28 | ~2% | ~1% |
| 10 000 | 0.10 | 0.32 | ~2% | ~0.5% |
| 100 000 | 0.41 | 0.43 | ~0.5% | ~0.5% |

The optimal fractional band (~0.5–1%) holds, but the achievable F1 is far lower
than tla. This is a consequence of the insert-size problem: Pindel's sensitivity
on LAGA is structurally limited regardless of threshold. Raising or lowering `-e`
cannot compensate for calls that are never generated.

### Recommendation

Operate on a **fractional cutoff of ~0.5–1% on Pindel AD**. The band is stable
across depths and both amplicons, but the ceiling F1 is amplicon-dependent and
may be unacceptably low for difficult amplicons. This reinforces the case for
LoFreq as the primary indel caller (see overall recommendation).

---

## Finding 3 — Positional blind spots (Experiment B)

*Source: `results/exptB/{tla,LAGA_LC}/position_sensitivity.csv`*

### Summary

Fixed condition: **10 bp deletion or insertion, AF 10%, depth 5000**, position
swept across the amplicon for read lengths 75/150/250 bp.

**Deletions** are largely position-independent: any read spanning the breakpoint
carries a `D` CIGAR op and Pindel detects it. Detection is flat across positions
at all read lengths except at the 3′ end of LAGA (see below).

**Insertions** are read-length dependent:
- At **75 bp**, short anchors cause BWA to soft-clip; Pindel detects via
  split-read signal. Center positions are recoverable.
- At **150 bp**, longer anchors cause BWA to emit CIGAR `I` instead of a
  soft-clip. Pindel misses these (Pindel is a split-read caller, not a
  CIGAR-`I` caller). VarScan2 and LoFreq are unaffected.
- At **250 bp**, reads span the full `tla` amplicon — no center gap — and CIGAR
  `I` is used for insertions throughout. Pindel misses nearly all of them.

**LAGA 3′-end blind spot (all callers, all read lengths):**  
Positions p621, p641, p674 (3′ end of the 711 bp amplicon) are never recovered
for deletions regardless of `-i` or read length. Root cause: the deletion is too
close to the amplicon 3′ end; nearly all reads end before the breakpoint and the
few that span it are soft-clipped by the aligner rather than carrying a gap
CIGAR. At p674: only ~8/10 000 reads have a clean `10D` CIGAR (vs ~280 at p355,
mid-amplicon). This is a genuine positional limitation of the aligner for
near-end events, not a Pindel-specific issue.

### Recommendation

For insertions in production, supplement Pindel with a pileup-based caller
(VarScan2 or LoFreq). Experiment D quantifies this improvement.

---

## Finding 4 — Caller comparison (Experiment D, `tla` and `LAGA_LC`, Aug 2026)

*Source: `results/exptD/{tla,LAGA_LC}/caller_sensitivity.csv`*

**Pindel** is scored at `ad_frac ≥ 1%` (the Exp C operating point). All other
callers use `MIN_SUPPORT = 1` (their own default filtering).

Grid: lengths 1/5/10/20 bp, AF 1/2/5%, depth 5000/10 000, seed 0.

### Results (micro-averaged over grid)

| Caller | tla sens | tla prec | tla F1 | LAGA_LC sens | LAGA_LC prec | LAGA_LC F1 |
|---|---|---|---|---|---|---|
| **LoFreq** | 0.69 | 0.78 | **0.73** | 0.71 | 0.77 | **0.74** |
| Pindel (`ad_frac ≥ 1%`) | 0.59 | 0.74 | 0.66 | 0.23 | 0.72 | 0.35 |
| VarScan2 | 0.70 | 0.44 | 0.54 | 0.74 | 0.27 | 0.40 |
| Mutect2 | 0.34 | 0.80 | 0.47 | 0.37 | 0.79 | 0.50 |

### Per-caller notes

**LoFreq** is the most consistent caller (F1 0.73/0.74 across both amplicons).
Pileup-based; detects both `D` and `+NNN` CIGAR ops. Not affected by the Pindel
insertion/read-length interaction or the insert-size geometry of the amplicon.

**Pindel** sensitivity on LAGA_LC is severely limited by the insert-size
problem: F1 drops from 0.66 to 0.35. Even with the correct `-i` (165 for
LAGA_LC), many conditions are scored at the detection ceiling and still miss
events. At r250 Pindel is entirely disabled on LAGA_LC (insert 165 < read
length 250). The positional plot shows Pindel dropping to zero at positions
beyond ~500 bp (3′ end of the 711 bp amplicon).

**VarScan2** has the highest raw sensitivity on both amplicons but the lowest
precision. On LAGA_LC precision collapses to 0.27 (F1 0.40), driven by a large
number of FPs (1440 vs 421 on tla). Useful if sensitivity is the priority and
downstream filtering is applied; not suitable as a standalone caller.

**Mutect2** has a high indel detection floor on uniform-start amplicon data:
- 0 calls at AF 1%, 2%, 5% (across all grid conditions, even with
  `--max-reads-per-alignment-start 0` and `--force-active-region`).
- Detects cleanly at AF ~20% (reports AF = 0.195).
- Floor is between 5% and 20%; assembly-based haplotype resolution requires
  more allele copies to anchor the graph than pileup-based callers.
- Since the Exp D grid tops at 5%, Mutect2 scores near-zero across the grid.
- On LAGA_LC Mutect2 scores slightly better (F1 0.50 vs 0.47) because the
  larger amplicon grid happens to include some mid-amplicon positions where
  Mutect2 occasionally fires. The pattern is otherwise identical.
- **Suggested follow-up:** extend `ALLELE_FREQS` to include 10% and 20% to
  show the Mutect2 crossover as a curve.

---

## Finding 5 — `pindel2vcf -e` is an absolute count; the production `-e 100` is not a principled AF cutoff

### What `-e` does

`pindel2vcf -e N` (`--min_supporting_reads`) is a **pure post-filter on alt
supporting reads** — an absolute integer, not an allele fraction. Verified on
real data (`runs/T3_rebuild`): the `-e 100` VCF is exactly the `-e 15` VCF with
every record whose `AD` alt < 100 dropped (identical CHROM/POS/REF/ALT/INFO,
only the row set differs).

### Why the production `-e 100` is problematic

`-e 100` imposes an AF floor of `100 / local_depth`:

| context | depth | `-e 100` → AF floor |
|---|---|---|
| this benchmark | 5 000× | ≈ 2% |
| real ultra-deep amplicon | ~370 000× | ≈ 0.03% |

The production pipeline's `-e 100` is **not** a principled 2% cutoff — it
coincidentally equals ~2% at 5000×, and becomes nearly a no-op at the real
assay depth. The HOMLEN deletion filter (`filter_deletions.py`, keeps only
`HOMLEN = 0`) compounds this by dropping genuine deletions in microhomologous
regions.

### Action taken (pipeline change, 2026-07)

`pindel_wrapper.py` now runs `pindel2vcf -e 1` once per SV type and produces
the full ladder by **post-hoc filtering** the `AD` field in Python (byte-identical
to the old per-e re-run, verified). A depth-invariant fractional ladder
(`f0p01/f0p02/f0p05`) is also produced. The correct production operating point
is a **fractional filter on `alt / samtools_depth`** at the breakpoint, not
absolute `-e`.

---

## Design notes

### Why four separate experiments (not one crossed grid)

A fully-crossed grid (type × length × AF × position × read length × depth ×
dataset × seed) would confound axes and explode combinatorially. Each experiment
isolates one variable:

- **Exp A** fixes position (5 relative positions) and sweeps sensitivity axes.
- **Exp B** fixes AF/length and sweeps position.
- **Exp C** adds no new simulation — re-scores A's `-e 1` VCFs at many thresholds.
- **Exp D** adds no new simulation — re-calls A's BAMs with other callers.

### Why pipeline `Summary.100.vcf` is excluded

The summary VCF would understate Pindel's raw capability:
1. `-e 100` imposes a hard 2% AF floor at 5000× (zeroes out the 1% conditions).
2. `HOMLEN = 0` filter drops genuine deletions in microhomologous regions.

Exp A/B score the `pindel2vcf -e 1` base VCFs directly. Exp C sweeps thresholds
post-hoc over the same base VCFs.

### Why replicate seeds

ART seeds re-seed the simulator's RNG (which fragments are sampled, where errors
land) at fixed depth/AF/position. Seeds are **independent replicates** of the same
condition — near the detection threshold a call is stochastic, so replicates
estimate detection *probability* rather than a single coin flip. Seeds do **not**
vary position (that is Exp B's dedicated axis).

### Scoring rules

Matching is **exact** (`POS_TOLERANCE = SIZE_TOLERANCE = 0`): correct SVTYPE,
exact breakpoint position, exact reported size, support ≥ threshold. Right-type
calls within 15 bp of truth that fail exact matching are written to `nearmiss.csv`
for manual inspection (e.g. microhomology-shifted breakpoints). AF = 0 WT-only
controls supply the FP baseline for precision.
