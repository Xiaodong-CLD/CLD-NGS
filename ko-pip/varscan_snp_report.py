"""VarScan SNP amino-acid annotation + thin SNP report.

Two entry points:
  * ``annotate_mpileup2snp_snp()`` — amino-acid-annotates a VarScan
    ``mpileup2snp`` SNP table and writes ``<stem>.amino_annotated.vcf``.
  * ``build_snp_report()`` — emits ``<output_dir>/varscan_snp_report.tsv`` from
    that annotated table (falling back to the plain SNP table), with a per-row
    PASS_VARSCAN_SNP / FILTER status from the ``var_freq`` and
    ``total_informative`` thresholds.

"""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import subprocess
from typing import Iterable

SNP_REPORT_COLUMNS = [
    "sample_id",
    "chrom",
    "position",
    "ref",
    "alt",
    "consensus",
    "ref_support",
    "alt_support",
    "total_informative",
    "var_freq",
    "p_value",
    "qual_alt",
    "mapqual_alt",
    "alt_strands_plus",
    "alt_strands_minus",
    "amino_change",
    "mutation_type",
    "ref_amino",
    "var_amino",
    "ref_codon",
    "var_codon",
    "status",
    "filter_reasons",
    "notes",
]

OUTPUT_FILE = "varscan_snp_report.tsv"


# ---------------------------------------------------------------------------
# Local amino-acid annotation — modified copy of the legacy annotator
# (ngs/amino-anno/from_varscan.py, by Tongyu Liu), reworked to parse VarScan's
# packed mpileup2snp schema instead of the flat pileup2snp schema. Translation
# and missense/nonsense/nonstop logic are unchanged.
# ---------------------------------------------------------------------------
_AMINO_ANNO_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "..", "amino-anno")
_FA_TRANS = os.path.join(_AMINO_ANNO_DIR, "bin", "faTrans")

# Output column contract: flat 19-column VarScan layout + 6 amino columns,
# consumed downstream by name.
_FLAT_SNP_COLUMNS = [
    "Chrom", "Position", "Ref", "Cons", "Reads1", "Reads2", "VarFreq",
    "Strands1", "Strands2", "Qual1", "Qual2", "Pvalue", "MapQual1", "MapQual2",
    "Reads1Plus", "Reads1Minus", "Reads2Plus", "Reads2Minus", "VarAllele",
]
_AMINO_COLUMNS = ["AminoChange", "MutationType", "RefAmino", "VarAmino",
                  "RefCodon", "VarCodon"]

_AA_CODES = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly',
    'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn',
    'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser', 'T': 'Thr', 'V': 'Val',
    'W': 'Trp', 'Y': 'Tyr', 'Z': 'STOP',
}


def _read_fasta(filepath: str) -> dict:
    with open(filepath, "r") as f:
        name = f.readline().rstrip().lstrip(">").split(" ")[0]
        return {"name": name,
                "seq": "".join(line.rstrip().upper() for line in f)}


def _write_fasta(seq: str, seqname: str, out_dir: str) -> str:
    os.makedirs(out_dir, exist_ok=True)
    outpath = os.path.join(out_dir, seqname + ".fa")
    with open(outpath, "w") as f:
        f.write(f">{seqname}\n")
        f.write(seq + "\n")
    return outpath


def _translate(seq: str, tmp_dir: str) -> str:
    """NT -> AA via the amino-anno faTrans binary (same as the legacy script)."""
    nt_path = _write_fasta(seq, "nt_tmp", tmp_dir)
    aa_path = os.path.join(tmp_dir, "aa_tmp.fa")
    subprocess.call(shlex.split(f"{_FA_TRANS} {nt_path} {aa_path}"))
    aa_seq = _read_fasta(aa_path)["seq"]
    os.remove(nt_path)
    os.remove(aa_path)
    return aa_seq


def _nz(x: str) -> int:
    """1 if x is a positive integer string, else 0 (for strand-count fields)."""
    try:
        return 1 if int(x) > 0 else 0
    except (TypeError, ValueError):
        return 0


def annotate_mpileup2snp_snp(snp_path: str, reference_path: str,
                             output_dir: str) -> str:
    """Amino-acid-annotate a VarScan ``mpileup2snp`` SNP table.

    Packed fields parsed per row:
        col 3  Var                                   -> VarAllele
        col 4  Cons:Cov:Reads1:Reads2:Freq:P-value   -> Cons, Reads1, Reads2,
                                                        VarFreq, Pvalue
        col 5  StrandFilter:R1+:R1-:R2+:R2-:pval      -> Reads1/2 Plus/Minus,
                                                        Strands1/Strands2
    """
    os.makedirs(output_dir, exist_ok=True)
    ref = _read_fasta(reference_path)
    ref_seq = ref["seq"]
    ref_trans = _translate(ref_seq, output_dir)

    stem = os.path.basename(snp_path).rsplit(".", 1)[0]
    out_path = os.path.join(output_dir, stem + ".amino_annotated.vcf")
    header_cols = _FLAT_SNP_COLUMNS + _AMINO_COLUMNS
    n_written = 0
    with open(snp_path) as snp_f, open(out_path, "w") as of:
        of.write("\t".join(header_cols) + "\n")
        snp_f.readline()  # drop the mpileup2snp header
        for line in snp_f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6 or cols[0] != ref["name"]:
                continue
            packed = (cols[4].split(":") + [""] * 6)[:6]
            cons, _cov, reads1, reads2, var_freq, pvalue = packed
            sf = (cols[5].split(":") + [""] * 5)[:5]  # StrandFilter,R1+,R1-,R2+,R2-
            r1p, r1m, r2p, r2m = sf[1], sf[2], sf[3], sf[4]
            var_allele = cols[3]

            idx = int(cols[1]) - 1
            seq_l = list(ref_seq)
            seq_l[idx] = var_allele
            var_seq = "".join(seq_l)
            var_trans = _translate(var_seq, output_dir)

            amino_idx = idx // 3
            ref_aa = ref_trans[amino_idx]
            var_aa = var_trans[amino_idx]
            codon_start = idx - idx % 3
            if ref_aa != var_aa:
                amino_change = "true"
                if ref_aa != "Z" and var_aa == "Z":
                    mut_type = "nonsense"
                elif ref_aa == "Z" and var_aa != "Z":
                    mut_type = "nonstop"
                else:
                    mut_type = "missense"
            else:
                amino_change = "false"
                mut_type = "N/A"

            row = {
                "Chrom": cols[0], "Position": cols[1], "Ref": cols[2],
                "Cons": cons, "Reads1": reads1, "Reads2": reads2,
                "VarFreq": var_freq,
                "Strands1": str(_nz(r1p) + _nz(r1m)),
                "Strands2": str(_nz(r2p) + _nz(r2m)),
                "Qual1": "", "Qual2": "", "Pvalue": pvalue,
                "MapQual1": "", "MapQual2": "",
                "Reads1Plus": r1p, "Reads1Minus": r1m,
                "Reads2Plus": r2p, "Reads2Minus": r2m,
                "VarAllele": var_allele,
                "AminoChange": amino_change, "MutationType": mut_type,
                "RefAmino": _AA_CODES[ref_aa], "VarAmino": _AA_CODES[var_aa],
                "RefCodon": ref_seq[codon_start:codon_start + 3],
                "VarCodon": var_seq[codon_start:codon_start + 3],
            }
            of.write("\t".join(str(row[c]) for c in header_cols) + "\n")
            n_written += 1
    print(f"[Takeda Pipeline Message] annotate_mpileup2snp_snp wrote={n_written} "
          f"-> {out_path}")
    return out_path



def _try_int(x: str) -> int:
    try:
        return int(x)
    except (TypeError, ValueError):
        return 0


def _try_float(x: str) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def _parse_var_freq(token: str) -> float:
    """VarScan emits VarFreq as '12.34%'. Return a fraction in [0,1]."""
    s = (token or "").strip()
    if not s:
        return 0.0
    if s.endswith("%"):
        return _try_float(s[:-1]) / 100.0
    return _try_float(s)


def _read_varscan_rows(path: str) -> list[dict]:
    """Read the annotated VarScan SNP tab table. Returns [] if file is missing."""
    if not os.path.exists(path):
        return []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def _build_row(
    sample_id: str, chrom_filter: str | None, vs: dict,
    min_var_freq: float, min_total: int,
) -> dict | None:
    chrom = vs.get("Chrom", "")
    if chrom_filter and chrom != chrom_filter:
        return None

    ref_support = _try_int(vs.get("Reads1", "0"))
    alt_support = _try_int(vs.get("Reads2", "0"))
    total = ref_support + alt_support
    var_freq = _parse_var_freq(vs.get("VarFreq", ""))

    reasons: list[str] = []
    if total < min_total:
        reasons.append("low_total_informative")
    if var_freq < min_var_freq:
        reasons.append("low_var_freq")
    status = "PASS_VARSCAN_SNP" if not reasons else "FILTER"

    return {
        "sample_id": sample_id,
        "chrom": chrom,
        "position": vs.get("Position", ""),
        "ref": vs.get("Ref", ""),
        "alt": vs.get("VarAllele", ""),
        "consensus": vs.get("Cons", ""),
        "ref_support": ref_support,
        "alt_support": alt_support,
        "total_informative": total,
        "var_freq": f"{var_freq:.4f}",
        "p_value": vs.get("Pvalue", ""),
        "qual_alt": vs.get("Qual2", ""),
        "mapqual_alt": vs.get("MapQual2", ""),
        "alt_strands_plus": vs.get("Reads2Plus", ""),
        "alt_strands_minus": vs.get("Reads2Minus", ""),
        "amino_change": vs.get("AminoChange", ""),
        "mutation_type": vs.get("MutationType", ""),
        "ref_amino": vs.get("RefAmino", ""),
        "var_amino": vs.get("VarAmino", ""),
        "ref_codon": vs.get("RefCodon", ""),
        "var_codon": vs.get("VarCodon", ""),
        "status": status,
        "filter_reasons": ",".join(reasons),
        "notes": "",
    }


def _pick_input(varscan_dir: str, sample_id: str) -> tuple[str, str]:
    """Prefer the amino-annotated file; fall back to the plain SNP table.

    Returns (path, source_label).
    """
    annotated = os.path.join(varscan_dir, f"{sample_id}.snp.amino_annotated.vcf")
    plain = os.path.join(varscan_dir, f"{sample_id}.snp.vcf")
    if os.path.exists(annotated):
        return annotated, "amino_annotated"
    return plain, "plain_snp"


def build_snp_report(
    varscan_dir: str,
    sample_id: str,
    output_dir: str,
    *,
    chrom: str | None = None,
    min_var_freq: float = 0.05,
    min_total: int = 50,
) -> str:
    os.makedirs(output_dir, exist_ok=True)
    input_path, source_label = _pick_input(varscan_dir, sample_id)
    vs_rows = _read_varscan_rows(input_path)

    out_rows: list[dict] = []
    for vs in vs_rows:
        row = _build_row(sample_id, chrom, vs, min_var_freq, min_total)
        if row is None:
            continue
        out_rows.append(row)

    out_rows.sort(key=lambda r: (r["chrom"], _try_int(r["position"]), r["alt"]))

    out_path = os.path.join(output_dir, OUTPUT_FILE)
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=SNP_REPORT_COLUMNS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

    n_pass = sum(1 for r in out_rows if r["status"] == "PASS_VARSCAN_SNP")
    print(f"[Takeda Pipeline Message] varscan_snp_report: source={source_label} "
          f"in={len(vs_rows)} wrote={len(out_rows)} pass={n_pass} -> {out_path}")
    return out_path


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Thin VarScan SNP report (no overlap with indel pipeline).",
    )
    p.add_argument("--varscan-dir", required=True,
                   help="Directory containing <sample>.snp.vcf (and optionally <sample>.snp.amino_annotated.vcf).")
    p.add_argument("--sample-id", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--chrom", default=None,
                   help="Optional contig filter (e.g. Fut8). Default: include all chroms.")
    p.add_argument("--min-var-freq", type=float, default=0.05,
                   help="Minimum VarScan VarFreq (as fraction) for PASS. Default 0.05 (5%%).")
    p.add_argument("--min-total", type=int, default=50,
                   help="Minimum Reads1+Reads2 for PASS. Default 50.")
    return p


def main(argv: Iterable[str] | None = None) -> int:
    args = _build_parser().parse_args(list(argv) if argv is not None else None)
    build_snp_report(
        varscan_dir=args.varscan_dir,
        sample_id=args.sample_id,
        output_dir=args.output_dir,
        chrom=args.chrom,
        min_var_freq=args.min_var_freq,
        min_total=args.min_total,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
