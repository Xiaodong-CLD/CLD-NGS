#!/usr/bin/env python3
"""Compare VarScan variant calls between the full-plasmid and insert references.

The insert reference is an exact contiguous sub-region of the plasmid:
    plasmid_pos = insert_pos + INSERT_OFFSET   (1-based, + strand)
so insert-reference calls are translated into the plasmid coordinate frame and
matched against plasmid-reference calls on (position, Ref, VarAllele).

For each sample it reports:
  - MATCHED      : same variant from both references (shows VarFreq from each)
  - PLASMID-ONLY : in the insert-covered window but missing from the insert run
  - INSERT-ONLY  : called on insert ref but absent from plasmid run
  - BACKBONE     : plasmid-only calls OUTSIDE the insert window (expected)

Usage:
  python compare_refs.py PLASMID_RESULTS_DIR INSERT_RESULTS_DIR
e.g.
  python compare_refs.py <amp-pip-dev>/runs/T1_5_threads/Results <amp-pip-dev>/runs/T_insertref/Results
"""
import os
import sys
import glob

INSERT_OFFSET = 431          # plasmid_pos = insert_pos + 431 (0-based offset)
INSERT_LEN = 7196            # insert covers plasmid 432 .. 7627
INSERT_START = INSERT_OFFSET + 1
INSERT_END = INSERT_OFFSET + INSERT_LEN
FREQ_TOL = 0.5               # percentage-point tolerance when flagging VarFreq drift


def parse_varscan(path):
    """Return dict keyed by (pos, ref, varallele) -> VarFreq string."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            pos = int(f[idx["Position"]])
            ref = f[idx["Ref"]]
            alt = f[idx["VarAllele"]]
            freq = f[idx["VarFreq"]]
            out[(pos, ref, alt)] = freq
    return out


def freq_val(s):
    try:
        return float(s.rstrip("%"))
    except (ValueError, AttributeError):
        return None


def sample_names(results_dir):
    return sorted(
        d for d in os.listdir(results_dir)
        if os.path.isdir(os.path.join(results_dir, d))
    )


def compare_one(sample, plasmid_dir, insert_dir, kind):
    pfile = os.path.join(plasmid_dir, sample, "varscan-output", f"{sample}.{kind}.vcf")
    ifile = os.path.join(insert_dir, sample, "varscan-output", f"{sample}.{kind}.vcf")
    plasmid = parse_varscan(pfile)
    insert_raw = parse_varscan(ifile)
    # translate insert coords -> plasmid frame
    insert = {}
    for (pos, ref, alt), freq in insert_raw.items():
        insert[(pos + INSERT_OFFSET, ref, alt)] = freq

    pkeys = set(plasmid)
    ikeys = set(insert)
    matched = sorted(pkeys & ikeys)
    plasmid_only = sorted(pkeys - ikeys)
    insert_only = sorted(ikeys - pkeys)

    in_window = [k for k in plasmid_only if INSERT_START <= k[0] <= INSERT_END]
    backbone = [k for k in plasmid_only if not (INSERT_START <= k[0] <= INSERT_END)]

    print(f"\n=== {sample}  [{kind}] ===")
    print(f"  plasmid calls={len(plasmid)}  insert calls={len(insert)}  "
          f"matched={len(matched)}  plasmid-only-in-window={len(in_window)}  "
          f"insert-only={len(insert_only)}  backbone(outside)={len(backbone)}")

    drift = []
    for k in matched:
        pf, ifq = freq_val(plasmid[k]), freq_val(insert[k])
        if pf is not None and ifq is not None and abs(pf - ifq) > FREQ_TOL:
            drift.append((k, plasmid[k], insert[k]))
    if drift:
        print(f"  VarFreq drift > {FREQ_TOL}pp:")
        for (pos, ref, alt), pf, ifq in drift:
            print(f"    plasmid_pos {pos} {ref}->{alt}: plasmid {pf} vs insert {ifq}")

    if in_window:
        print("  PLASMID-ONLY inside insert window (investigate):")
        for pos, ref, alt in in_window:
            print(f"    plasmid_pos {pos} {ref}->{alt}  (insert_pos {pos - INSERT_OFFSET})  freq {plasmid[(pos,ref,alt)]}")
    if insert_only:
        print("  INSERT-ONLY (investigate):")
        for pos, ref, alt in insert_only:
            print(f"    plasmid_pos {pos} {ref}->{alt}  (insert_pos {pos - INSERT_OFFSET})  freq {insert[(pos,ref,alt)]}")
    return len(matched), len(in_window), len(insert_only)


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: compare_refs.py PLASMID_RESULTS_DIR INSERT_RESULTS_DIR")
    plasmid_dir, insert_dir = sys.argv[1], sys.argv[2]
    print(f"insert covers plasmid positions {INSERT_START}..{INSERT_END} "
          f"(plasmid_pos = insert_pos + {INSERT_OFFSET})")
    samples = sample_names(insert_dir)
    totals = {"matched": 0, "in_window": 0, "insert_only": 0}
    for sample in samples:
        for kind in ("snp", "indel"):
            m, w, io = compare_one(sample, plasmid_dir, insert_dir, kind)
            totals["matched"] += m
            totals["in_window"] += w
            totals["insert_only"] += io
    print("\n=== TOTALS ===")
    print(f"  matched={totals['matched']}  "
          f"plasmid-only-in-window={totals['in_window']}  "
          f"insert-only={totals['insert_only']}")
    print("  (plasmid-only-in-window and insert-only should ideally be 0 or "
          "explained by edge effects)")


if __name__ == "__main__":
    main()
