#!/usr/bin/env python3
"""Map insert-reference coordinates onto the full-plasmid reference.

The insert FASTA is a sub-region of the circular plasmid, possibly rotated
(different linear start) and possibly reverse-complemented. This script finds
how the insert sequence sits on the plasmid by:

  1. Doubling the plasmid (plasmid+plasmid) to absorb the circular origin.
  2. Searching for the insert (and its reverse complement) as a substring.
  3. If a single contiguous match is found, reporting the offset + strand so a
     1-based position on the insert can be translated to a 1-based position on
     the plasmid.

If no exact contiguous match exists (edits / SNPs between the two), it falls
back to anchoring on the first and last 40 bp of the insert independently and
reports both anchors so a piecewise/linear map can still be derived.

Usage:
  python ref_coord_map.py PLASMID.fa INSERT.fa
"""
import sys


def read_fasta_seq(path):
    seq = []
    name = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                continue
            seq.append(line.strip())
    return name, "".join(seq).upper()


def revcomp(s):
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(s))


def find_unique(haystack, needle):
    """Return list of 0-based start offsets of needle in haystack."""
    offsets = []
    start = 0
    while True:
        i = haystack.find(needle, start)
        if i < 0:
            break
        offsets.append(i)
        start = i + 1
    return offsets


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: ref_coord_map.py PLASMID.fa INSERT.fa")
    p_name, plasmid = read_fasta_seq(sys.argv[1])
    i_name, insert = read_fasta_seq(sys.argv[2])
    P = len(plasmid)
    I = len(insert)
    print(f"plasmid contig : {p_name}  length={P}")
    print(f"insert  contig : {i_name}  length={I}")

    doubled = plasmid + plasmid  # absorb circular origin

    for strand, query in (("+", insert), ("-", revcomp(insert))):
        hits = find_unique(doubled, query)
        # keep only hits that start within the first copy
        hits = [h for h in hits if h < P]
        if hits:
            off = hits[0]
            print()
            print(f"EXACT contiguous match on strand {strand}")
            print(f"  insert 1-based pos k  ->  plasmid 1-based pos:")
            if strand == "+":
                print(f"    plasmid_pos = ((k - 1 + {off}) % {P}) + 1")
            else:
                # insert base k corresponds to plasmid base at off + (I - k)
                print(f"    plasmid_pos = (({off} + {I} - k) % {P}) + 1   (minus strand)")
            print(f"  start offset (0-based) = {off}; wraps origin = {off + I > P}")
            return

    # Fallback: anchor on first/last 40 bp
    print()
    print("No exact contiguous match (insert and plasmid differ by edits/SNPs).")
    print("Falling back to end anchors (40 bp each):")
    for label, frag in (("5' anchor", insert[:40]), ("3' anchor", insert[-40:])):
        for strand, q in (("+", frag), ("-", revcomp(frag))):
            hits = [h for h in find_unique(doubled, q) if h < P]
            if hits:
                print(f"  {label} strand {strand}: plasmid 0-based offset {hits[0]}")


if __name__ == "__main__":
    main()
