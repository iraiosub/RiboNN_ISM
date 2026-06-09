#!/usr/bin/env python3
"""
prepare_scn2a_ism.py – Generate RiboNN prediction input for 5'UTR ISM.

Extracts a gene's 5'UTR / CDS / 3'UTR sequences from a genome FASTA + GTF,
then generates every SNV and deletion variant in the N bases immediately
upstream of the canonical AUG start codon.

Output: data/prediction_input.txt  (RiboNN TSV format)

Columns: tx_id | utr5_sequence | cds_sequence | utr3_sequence

After running this script, predict with:
    python3 run_ribonn_predict.py          # wrapper that calls src.predict directly
    # -- or --
    cp data/prediction_input.txt data/prediction_input1.txt && make predict_human

Then plot with:
    python plot_te_changes.py

Usage:
    python prepare_scn2a_ism.py \\
        --species human \\
        --gene-name SCN2A \\
        --fasta /path/to/GRCh38.primary_assembly.genome.fa \\
        --gtf   /path/to/gencode.v44.primary_assembly.annotation.gtf.gz \\
        [--upstream-bases 15] [--max-deletion 15] [--truncate-utr3]
"""

import argparse
import csv
import gzip
import os
from pathlib import Path

import numpy as np
import pyfaidx

# ── defaults matching the CAMP/NEMO shared reference location ────────────────
DEFAULT_REF_ROOT = "/camp/lab/ulej/home/shared/oscar_ira_riboloco/ref"
DEFAULT_REFS = {
    "human": {
        "fasta": f"{DEFAULT_REF_ROOT}/human/GRCh38.primary_assembly.genome.fa",
        "gtf": f"{DEFAULT_REF_ROOT}/human/gencode.v44.primary_assembly.annotation.longest_cds_transcripts.gtf.gz",
    },
    "mouse": {
        "fasta": f"{DEFAULT_REF_ROOT}/mouse/GRCm39.primary_assembly.genome.fa",
        "gtf": f"{DEFAULT_REF_ROOT}/mouse/gencode.vM33.primary_assembly.annotation.longest_cds_transcripts.gtf.gz",
    },
}

BASES      = ("A", "C", "G", "T")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# RiboNN hard limits (from src/predict.py and README)
MAX_UTR5_LEN     = 1_381
MAX_CDS_UTR3_LEN = 11_937


# ── GTF parsing (adapted from ag-play/ism_scn2a.py) ─────────────────────────

def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _parse_gtf_attrs(value):
    attrs = {}
    for item in value.strip().rstrip(";").split(";"):
        item = item.strip()
        if item and " " in item:
            k, v = item.split(" ", 1)
            attrs[k] = v.strip().strip('"')
    return attrs


def _merge(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def _interval_len(intervals):
    return sum(e - s for s, e in intervals)


def load_gene_transcript(gtf_path, gene_name, transcript_id=None):
    """Return the longest-CDS transcript for gene_name from the GTF."""
    transcripts = {}
    target_gene = gene_name.upper()
    keep = {"exon", "CDS", "start_codon", "stop_codon"}
    with _open(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _src, feature, start, end, _sc, strand, _fr, attrs = fields
            if feature not in keep:
                continue
            parsed = _parse_gtf_attrs(attrs)
            parsed_gene = parsed.get("gene_name", "")
            if parsed_gene.upper() != target_gene:
                continue
            tid = parsed.get("transcript_id")
            if not tid:
                continue
            if transcript_id and tid != transcript_id and tid.split(".")[0] != transcript_id.split(".")[0]:
                continue
            tx = transcripts.setdefault(
                tid,
                {"chrom": chrom, "strand": strand,
                 "gene_name": parsed.get("gene_name", ""),
                 "transcript_id": tid,
                 "transcript_name": parsed.get("transcript_name", ""),
                 "exons": [], "cds": [], "start_codons": [], "stop_codons": []},
            )
            s0, e0 = int(start) - 1, int(end)
            if feature == "exon":
                tx["exons"].append((s0, e0))
            elif feature == "CDS":
                tx["cds"].append((s0, e0))
            elif feature == "start_codon":
                tx["start_codons"].append((s0, e0))
            elif feature == "stop_codon":
                tx["stop_codons"].append((s0, e0))

    candidates = []
    for tx in transcripts.values():
        for k in ("exons", "cds", "start_codons", "stop_codons"):
            tx[k] = _merge(tx[k])
        if tx["exons"] and tx["cds"]:
            candidates.append(tx)
    if not candidates:
        req = f" transcript {transcript_id}" if transcript_id else ""
        raise ValueError(f"No {gene_name}{req} with exons+CDS found in {gtf_path}")
    candidates.sort(key=lambda t: (_interval_len(t["cds"]), _interval_len(t["exons"])), reverse=True)
    return candidates[0]


def _rc(seq):
    return seq.translate(COMPLEMENT)[::-1].upper()


def transcript_sequence(fasta, transcript):
    """Return (tx_seq, coords_array) for the transcript."""
    chrom = transcript["chrom"]
    if transcript["strand"] == "+":
        exons = sorted(transcript["exons"])
        seq_chunks   = [fasta[chrom][s:e].seq.upper() for s, e in exons]
        coord_chunks = [np.arange(s, e, dtype=np.int64) for s, e in exons]
    else:
        exons = sorted(transcript["exons"], reverse=True)
        seq_chunks   = [_rc(fasta[chrom][s:e].seq) for s, e in exons]
        coord_chunks = [np.arange(e - 1, s - 1, -1, dtype=np.int64) for s, e in exons]
    return "".join(seq_chunks), np.concatenate(coord_chunks)


def _interval_mask(coords, intervals):
    mask = np.zeros(len(coords), dtype=bool)
    for s, e in intervals:
        mask |= (coords >= s) & (coords < e)
    return mask


def transcript_regions(transcript, coords):
    tx_pos    = np.arange(len(coords), dtype=np.int64)
    cds_mask  = _interval_mask(coords, transcript["cds"])
    start_mask = _interval_mask(coords, transcript["start_codons"])
    stop_mask  = _interval_mask(coords, transcript["stop_codons"])
    cds_idx   = np.flatnonzero(cds_mask)
    if len(cds_idx) == 0:
        raise ValueError("No CDS bases found")
    canon_start = int(np.flatnonzero(start_mask).min()) if start_mask.any() else int(cds_idx.min())
    canon_stop_end = int(np.flatnonzero(stop_mask).max() + 1) if stop_mask.any() else int(cds_idx.max() + 1)
    return {"canonical_start": canon_start, "canonical_stop_end": canon_stop_end}


# ── variant generation ───────────────────────────────────────────────────────

def split_transcript(tx_seq, canonical_start, canonical_stop_end):
    """Split full transcript sequence into (utr5, cds, utr3)."""
    utr5 = tx_seq[:canonical_start]
    cds  = tx_seq[canonical_start:canonical_stop_end]
    utr3 = tx_seq[canonical_stop_end:]
    return utr5, cds, utr3


def generate_variants(utr5, cds, utr3, upstream_bases, max_deletion):
    """
    Yield (tx_id, utr5_seq, cds_seq, utr3_seq) rows for all ISM variants.

    SNV variants:   every base in the last `upstream_bases` of utr5 → 3 alt bases each
    del1 variants:  delete one base at each of the last `upstream_bases` positions
    del_N variants: delete the last N bases of utr5 (growing window, N = 2..max_deletion)
    """
    utr5_len = len(utr5)
    upstream_bases = min(upstream_bases, utr5_len)  # guard against values > UTR length

    # Reference
    yield "reference", utr5, cds, utr3

    # SNVs + single-base deletions at each upstream position
    for k in range(1, upstream_bases + 1):
        pos = utr5_len - k          # position in utr5 string (0-based)
        offset = -k                 # offset from canonical ATG
        ref_base = utr5[pos]

        # SNVs
        if ref_base in BASES:
            for alt in BASES:
                if alt == ref_base:
                    continue
                new_utr5 = utr5[:pos] + alt + utr5[pos + 1:]
                yield f"{offset:+d}_{ref_base}>{alt}", new_utr5, cds, utr3

        # Single-base deletion
        new_utr5 = utr5[:pos] + utr5[pos + 1:]
        yield f"del1_{offset:+d}_{ref_base}", new_utr5, cds, utr3

    # Growing deletions ending at the canonical AUG
    for n in range(2, max_deletion + 1):
        if n > utr5_len:
            break
        new_utr5 = utr5[:-n]
        yield f"del_{n}_to_-1", new_utr5, cds, utr3


# ── length checks ────────────────────────────────────────────────────────────

def check_and_truncate(utr5, cds, utr3, truncate_utr3):
    """
    Check RiboNN length limits. Optionally truncate utr3 to fit.
    Returns (utr5, cds, utr3, warning_str).
    """
    warn = ""
    cds_utr3 = len(cds) + len(utr3)

    if len(utr5) > MAX_UTR5_LEN:
        warn += f"[WARN] 5'UTR length {len(utr5)} exceeds RiboNN limit ({MAX_UTR5_LEN}); variants with longer UTR5 will be excluded.\n"

    if cds_utr3 > MAX_CDS_UTR3_LEN:
        if truncate_utr3:
            allowed_utr3 = MAX_CDS_UTR3_LEN - len(cds)
            if allowed_utr3 < 0:
                raise ValueError(
                    f"CDS alone ({len(cds)} nt) exceeds RiboNN CDS+3'UTR limit ({MAX_CDS_UTR3_LEN}). "
                    "Cannot truncate. Choose a shorter transcript."
                )
            utr3 = utr3[:allowed_utr3]
            warn += f"[INFO] 3'UTR truncated to {allowed_utr3} nt to fit RiboNN limit.\n"
        else:
            warn += (
                f"[WARN] CDS ({len(cds)} nt) + 3'UTR ({len(utr3)} nt) = {cds_utr3} nt "
                f"exceeds RiboNN limit ({MAX_CDS_UTR3_LEN}). "
                "All variants will be excluded! Re-run with --truncate-utr3.\n"
            )
    return utr5, cds, utr3, warn


# ── CLI ──────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--species", choices=sorted(DEFAULT_REFS), default="human",
                        help="Reference species for default FASTA/GTF paths (default: human)")
    parser.add_argument("--fasta", default=None,
                        help="Genome FASTA. Defaults to RIBONN_<SPECIES>_FASTA, RIBONN_FASTA, or the CAMP ref/<species> path.")
    parser.add_argument("--gtf", default=None,
                        help="GENCODE annotation GTF. Defaults to RIBONN_<SPECIES>_GTF, RIBONN_GTF, or the CAMP ref/<species> path.")
    parser.add_argument("--gene-name",     default="SCN2A")
    parser.add_argument("--transcript-id", default=None,
                        help="Specific transcript ID (default: longest CDS)")
    parser.add_argument("--upstream-bases", type=int, default=15,
                        help="Bases before AUG to test with SNVs and single-base deletions (default: 15)")
    parser.add_argument("--max-deletion",   type=int, default=15,
                        help="Max growing-deletion length (default: 15)")
    parser.add_argument("--output",  default="data/prediction_input.txt",
                        help="Output TSV path (default: data/prediction_input.txt)")
    parser.add_argument("--truncate-utr3", action="store_true",
                        help="Truncate 3'UTR to satisfy RiboNN's 11,937 nt CDS+3'UTR limit")
    parser.add_argument("--audit", action="store_true",
                        help="Print sequence info and exit without writing output")
    return parser.parse_args()


def discover_reference_path(species, kind):
    ref_root = Path(os.environ.get("RIBONN_REF_ROOT", DEFAULT_REF_ROOT))
    species_dir = ref_root / species

    if species == "human" and kind == "fasta":
        default_name = "GRCh38.primary_assembly.genome.fa"
    elif species == "human" and kind == "gtf":
        default_name = "gencode.v44.primary_assembly.annotation.longest_cds_transcripts.gtf.gz"
    elif species == "mouse" and kind == "fasta":
        default_name = "GRCm39.primary_assembly.genome.fa"
    else:
        default_name = "gencode.vM33.primary_assembly.annotation.longest_cds_transcripts.gtf.gz"

    default_path = species_dir / default_name
    if default_path.is_file():
        return str(default_path)

    patterns = (
        ("GRC*.primary_assembly.genome.fa", "*.primary_assembly.genome.fa", "*.fa")
        if kind == "fasta"
        else (
            "gencode.v*.primary_assembly.annotation.longest_cds_transcripts.gtf.gz",
            "*longest_cds_transcripts*.gtf.gz",
            "*.gtf.gz",
            "*.gtf",
        )
    )
    for pattern in patterns:
        for candidate in sorted(species_dir.glob(pattern)):
            if candidate.is_file():
                return str(candidate)

    return str(default_path)


def resolve_reference_paths(args):
    species_key = args.species.upper()
    fasta = (
        args.fasta
        or os.environ.get(f"RIBONN_{species_key}_FASTA")
        or os.environ.get("RIBONN_FASTA")
        or discover_reference_path(args.species, "fasta")
    )
    gtf = (
        args.gtf
        or os.environ.get(f"RIBONN_{species_key}_GTF")
        or os.environ.get("RIBONN_GTF")
        or discover_reference_path(args.species, "gtf")
    )
    return fasta, gtf


def normalize_gene_symbol(gene_name, species):
    if species == "human":
        return gene_name.upper()
    if species == "mouse":
        lower = gene_name.lower()
        return lower[:1].upper() + lower[1:]
    return gene_name


def main():
    args = parse_args()
    fasta_path, gtf_path = resolve_reference_paths(args)
    gene_name = normalize_gene_symbol(args.gene_name, args.species)

    print(f"Loading transcript for {gene_name} ({args.species}) ...")
    print(f"FASTA      : {fasta_path}")
    print(f"GTF        : {gtf_path}")
    fasta = pyfaidx.Fasta(fasta_path)
    transcript = load_gene_transcript(gtf_path, gene_name, args.transcript_id)
    tx_seq, coords = transcript_sequence(fasta, transcript)
    regions = transcript_regions(transcript, coords)
    fasta.close()

    canonical_start    = regions["canonical_start"]
    canonical_stop_end = regions["canonical_stop_end"]

    utr5_ref, cds, utr3 = split_transcript(tx_seq, canonical_start, canonical_stop_end)
    utr5_ref, cds, utr3, length_warning = check_and_truncate(
        utr5_ref, cds, utr3, args.truncate_utr3
    )

    print(f"Transcript : {transcript['transcript_id']} ({transcript['transcript_name']})")
    print(f"Strand     : {transcript['strand']}  Chrom: {transcript['chrom']}")
    print(f"5'UTR      : {len(utr5_ref)} nt")
    print(f"CDS        : {len(cds)} nt  (includes start + stop codons)")
    print(f"3'UTR      : {len(utr3)} nt")
    actual_upstream = min(args.upstream_bases, len(utr5_ref))

    print(f"Upstream bases requested: {args.upstream_bases}")
    print(f"Upstream bases to mutate/delete: {actual_upstream}")
    print(f"Max growing deletion: {args.max_deletion} bp")

    # Show the sequence window being mutated
    mutated_seq = utr5_ref[-actual_upstream:] if actual_upstream else ""
    print(f"\nLast {actual_upstream} bases of 5'UTR (ISM window):")
    for i, base in enumerate(mutated_seq):
        offset = -(actual_upstream - i)
        print(f"  offset {offset:+3d}  utr5_pos {len(utr5_ref) - actual_upstream + i}  base={base}")

    print(f"\nFirst 6 nt of CDS (should start ATG): {cds[:6]}")
    if not cds.startswith("ATG"):
        print("[WARN] CDS does not start with ATG — check transcript parsing.")

    if length_warning:
        print()
        print(length_warning.strip())

    if args.audit:
        print("\n--audit flag set; exiting without writing output.")
        return

    # Count variants
    snv_count  = sum(1 for k in range(1, actual_upstream + 1)
                     if utr5_ref[len(utr5_ref) - k] in BASES) * 3
    del1_count = actual_upstream
    delN_count = max(0, min(args.max_deletion, len(utr5_ref)) - 1)
    total      = 1 + snv_count + del1_count + delN_count

    print(f"\nVariants: 1 reference + {snv_count} SNVs + {del1_count} del1 + {delN_count} del_N = {total} rows")

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    written = 0
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["tx_id", "utr5_sequence", "cds_sequence", "utr3_sequence"])
        for tx_id, utr5, cds_row, utr3_row in generate_variants(
            utr5_ref, cds, utr3, args.upstream_bases, args.max_deletion
        ):
            writer.writerow([tx_id, utr5, cds_row, utr3_row])
            written += 1

    print(f"\nWrote {written} rows to {out_path}")
    print("\nNext steps:")
    print(f"  python run_ribonn_predict.py --species {args.species} --input {out_path}")
    print(f"  python plot_te_changes.py --gene-name {gene_name}")


if __name__ == "__main__":
    main()
