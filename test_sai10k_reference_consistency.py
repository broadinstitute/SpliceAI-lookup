#!/usr/bin/env python3
"""
Reference-consistency tests for SAI-10k predictions.

For each test variant, this script:
  1. Queries the production SpliceAI API for per-transcript delta scores.
  2. Selects the top-priority transcript (MS > MP > C > N, tie-break by sum
     of |delta scores|) -- matches sai10k_select_transcript.
  3. Runs THREE implementations on the SAME scores + transcript structure:
       a. google_cloud_run_services/sai10k_predictions.py (new Python impl;
          called via in-process Python API).
       b. SAI-10k-calc/spliceai_parser.py (reference Python parser; invoked
          as a subprocess on a synthesized VCF + refseq_table).
       c. SAI-10k-calc/spliceAI_parser.R (reference R parser; invoked via
          Rscript on the same inputs).
  4. Compares: aberration TYPE SET, affected region (exon/intron #), and
     frameshift YES/NO per aberration.

Known divergences:
  * The new impl computes frameshift from CDS-overlapping bases;
    spliceai_parser.* use the full affected-segment size. Mismatches on
    frameshift alone are reported but do not fail tests by default.
  * The reference parsers filter out indels (REF or ALT len > 1) at
    preprocessing; indel variants still run through our new impl but are
    reported as "ref skipped" in the comparison output.
  * The R reference requires BSgenome.Hsapiens.UCSC.hg19 / hg38 bioconductor
    packages. If not installed, the R comparison is skipped with a warning.
  * hg19.fa must be on disk for hg37 Python reference runs. If missing,
    the script falls back to hg38.fa (AA seqs will be wrong but
    aberration/region/frameshift flags are unaffected).

Usage:
  python3 test_sai10k_reference_consistency.py                 # run all tests
  python3 test_sai10k_reference_consistency.py --dump-only     # print per-variant diffs and exit
  python3 -m unittest test_sai10k_reference_consistency -v     # unittest mode
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

import requests


# ---- Paths -----------------------------------------------------------------

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "google_cloud_run_services"))
from sai10k_predictions import sai10k_compute_predictions, sai10k_select_transcript  # noqa: E402

REF_DIR_CANDIDATES = [
    os.path.expanduser("~/code/SpliceAI-lookup-dev/SAI-10k-calc"),
    os.path.expanduser("~/code/SAI-10k-calc"),
]
REF_DIR = next((p for p in REF_DIR_CANDIDATES if os.path.isdir(p)), None)
REF_PYTHON = os.path.join(REF_DIR, "spliceai_parser.py") if REF_DIR else None
REF_R = os.path.join(REF_DIR, "spliceAI_parser.R") if REF_DIR else None

# Locate FASTA files (prefer project-root symlinks, then home dir).
FASTA_SEARCH = {
    "37": [
        os.path.join(PROJECT_ROOT, "hg19.fa"),
        os.path.expanduser("~/hg19.fa"),
    ],
    "38": [
        os.path.join(PROJECT_ROOT, "hg38.fa"),
        os.path.expanduser("~/hg38.fa"),
    ],
}


_FASTA_FALLBACK_WARNED = set()


def _find_fasta(hg):
    for p in FASTA_SEARCH[hg]:
        if os.path.isfile(p):
            return p
    # Last resort: fall back to the other-assembly FASTA (AA seqs will be
    # wrong, but aberration/region/frameshift flags are computed without
    # reading sequence). Warn once per missing-assembly so users know.
    other = "38" if hg == "37" else "37"
    for p in FASTA_SEARCH[other]:
        if os.path.isfile(p):
            if hg not in _FASTA_FALLBACK_WARNED:
                print(
                    f"WARNING: no FASTA found for hg{hg}; falling back to {p} "
                    f"(hg{other}). Aberration/region/frameshift flags unaffected; "
                    "any AA-seq output from the reference parser will be wrong.",
                    file=sys.stderr,
                )
                _FASTA_FALLBACK_WARNED.add(hg)
            return p
    return None


REF_NAME = {"37": "hg19", "38": "hg38"}

RSCRIPT = shutil.which("Rscript")
HAS_R = RSCRIPT is not None

API_URL_TEMPLATE = (
    "https://spliceai-{hg}-xwkwwwxdwq-uc.a.run.app/spliceai/"
    "?hg={hg}&distance=500&mask=0&variant={variant}&raw={variant}"
)


# ---- Variant superset ------------------------------------------------------
# Collected from:
#   test_api_consistency.py (VARIANTS + VARIANTS_HG37_ONLY)
#   test_spliceai.py        (indel + SNV hg38 examples)
#   test_ui.py              (smoke-test variants; mostly dup with above)
#   google_cloud_run_services/test_sai10k_predictions.py  (BRCA1/2 hg19 natives)

# (variant, genome_versions_to_test) where genome_versions is a tuple.
# HGVS:-prefixed variants are resolved via the Ensembl VEP API at runtime.
VARIANT_SOURCES = [
    # -- test_api_consistency.py VARIANTS (both hg37 + hg38) --
    ("8-140300616-T-G",        ("37", "38"), "TRAPPC9 (index.html + test_ui)"),
    ("6-31740453-G-T",         ("37", "38"), "index.html space-separated"),
    ("1-930130-C-G",           ("37", "38"), "index.html overlapping genes"),
    ("1-1042601-A-AGAGAG",     ("37", "38"), "indel (insertion)"),
    ("1-1042466-GGGC-G",       ("37", "38"), "indel (deletion)"),
    ("17-58709857-A-C",        ("37", "38"), "RAD51C exon skipping + partial deletion"),
    ("22-28695127-C-T",        ("37", "38"), "CHEK2 exon skipping"),
    ("16-23641111-T-C",        ("37", "38"), "PALB2 partial exon deletion"),
    ("17-58696687-T-A",        ("37", "38"), "RAD51C partial intron retention"),
    ("22-28741768-C-A",        ("37", "38"), "CHEK2 non-coding partial exon deletion"),
    ("13-32345247-T-G",        ("37", "38"), "BRCA2 pseudoexon"),
    ("13-32363534-G-T",        ("37", "38"), "BRCA2 exon skipping"),
    ("2-47790924-C-CAGTTG",    ("37", "38"), "indel (insertion)"),
    ("1-55039916-G-A",         ("37", "38"), "index.html source comment"),
    ("17-43047638-C-G",        ("37", "38"), "BRCA1 c.5467+5G>C (hg38 liftover)"),
    ("17-43047642-C-T",        ("37", "38"), "BRCA1 c.5467+1G>A (hg38 liftover)"),
    ("17-43045714-G-C",        ("37", "38"), "BRCA1 region hg38"),
    ("X-154031414-C-T",        ("37", "38"), "MECP2 region"),
    ("HGVS:NM_000249.4:c.116G>A", ("37", "38"), "MLH1 HGVS"),

    # -- test_api_consistency.py VARIANTS_HG37_ONLY (hg19-native) --
    ("17-41199655-C-G",        ("37",),      "BRCA1 c.5467+5G>C (hg19 native)"),
    ("17-41199659-C-T",        ("37",),      "BRCA1 c.5467+1G>A (hg19 native)"),

    # -- test_sai10k_predictions.py BRCA2 hg19-native --
    ("13-32889805-G-A",        ("37",),      "BRCA2 c.-40+1G>A"),
    ("13-32890665-G-A",        ("37",),      "BRCA2 c.67+1G>A"),
    ("13-32890637-A-G",        ("37",),      "BRCA2 c.40A>G"),

    # -- test_spliceai.py (hg38 SNV + indel) --
    ("1-69091-A-AA",           ("38",),      "OR4F5 indel (insertion)"),
    ("1-69539-T-G",            ("38",),      "OR4F5 SNV"),

    # -- test_ui.py (mostly overlap; this one is unannotated / will skip) --
    ("1-10000-A-G",            ("38",),      "unannotated intergenic (expected: no scores)"),
]


def _is_indel(variant):
    parts = variant.split("-")
    if len(parts) < 4 or variant.startswith("HGVS:"):
        return False
    ref, alt = parts[-2], parts[-1]
    return len(ref) != 1 or len(alt) != 1


# ---- HGVS resolution ------------------------------------------------------
def normalize_hgvs(hgvs_string, hg):
    """Resolve an HGVS string to chrom-pos-ref-alt via Ensembl VEP."""
    prefix = "grch37." if hg == "37" else ""
    url = (
        f"https://{prefix}rest.ensembl.org/vep/human/hgvs/{hgvs_string}"
        "?content-type=application/json&vcf_string=1"
    )
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    data = resp.json()
    if not isinstance(data, list) or not data:
        raise ValueError(f"Unexpected VEP response for {hgvs_string}: {data}")
    vcf_string = data[0].get("vcf_string") if isinstance(data[0], dict) else None
    if not vcf_string:
        raise ValueError(f"VEP response missing vcf_string for {hgvs_string}: {data[0]}")
    return vcf_string


# ---- API query + transcript selection -------------------------------------
def query_spliceai_api(variant, hg):
    url = API_URL_TEMPLATE.format(hg=hg, variant=variant)
    resp = requests.get(url, timeout=120)
    resp.raise_for_status()
    return resp.json()


def parse_variant_coords(variant):
    """variant "chr-pos-ref-alt" -> (chrom, pos, ref, alt)."""
    c, p, r, a = variant.split("-", 3)
    if c.startswith("chr"):
        c = c[3:]
    return c, int(p), r, a


# ---- refseq_table construction (for reference parser input) ---------------
def _compute_exon_frames(exon_starts_1b, exon_ends_1b, cds_start_1b, cds_end_1b):
    """Compute per-exon frames (genePred order).

    Sets -1 for exons that don't overlap the CDS (UTR-only), 0 otherwise.
    0 is a safe default for coding exons -- it preserves the filter
    `eFrame != -1` used by the reference parser to drop UTR exons.
    (AA seqs will be wrong for frames that aren't actually 0, but
    aberration flags and region info are independent of eFrame.)
    """
    frames = []
    cds_valid = (
        cds_start_1b is not None and cds_end_1b is not None
        and cds_start_1b < cds_end_1b
    )
    for start, end in zip(exon_starts_1b, exon_ends_1b):
        if not cds_valid:
            frames.append(-1)
        elif end < cds_start_1b or start > cds_end_1b:
            frames.append(-1)
        else:
            frames.append(0)
    return frames


def write_refseq_table(path, chrom, transcript):
    """Write a 1-row refseq_table TSV for the selected transcript.

    Converts API's 1-based exon coords back to genePred 0-based starts.
    """
    exon_starts_1b = transcript["EXON_STARTS"]
    exon_ends_1b = transcript["EXON_ENDS"]
    cds_start_1b = transcript.get("CDS_START")
    cds_end_1b = transcript.get("CDS_END")

    starts_0b = [s - 1 for s in exon_starts_1b]
    ends_0b = list(exon_ends_1b)
    tx_start_0b = starts_0b[0]
    tx_end_0b = ends_0b[-1]
    # Mirror the new impl's coding check (sai10k_predictions.py: cds_start <
    # cds_end). An API that returns cds_start == cds_end signals non-coding;
    # encoding it as a 1bp CDS would make the ref parser treat it as coding
    # while the new impl treats it as non-coding -- an asymmetry that could
    # flip every frameshift comparison for non-coding transcripts.
    if (cds_start_1b is not None and cds_end_1b is not None
            and cds_start_1b < cds_end_1b):
        cds_start_0b = cds_start_1b - 1
        cds_end_0b = cds_end_1b
    else:
        # genePred non-coding convention: cdsStart == cdsEnd.
        cds_start_0b = tx_end_0b
        cds_end_0b = tx_end_0b

    name = transcript["t_id"]
    if transcript.get("t_refseq_ids"):
        name = transcript["t_refseq_ids"][0]

    exon_frames = _compute_exon_frames(
        exon_starts_1b, exon_ends_1b, cds_start_1b, cds_end_1b
    )
    gene = transcript.get("g_name") or name

    header = [
        "bin", "name", "refseq_nov", "refseq_version", "chrom", "strand",
        "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount",
        "exonStarts", "exonEnds", "score", "name2",
        "cdsStartStat", "cdsEndStat", "exonFrames",
        "Gene", "RefSeq_ID", "refseq_match",
    ]
    row = [
        "0", name,
        name.split(".")[0],
        name.split(".")[-1] if "." in name else "1",
        str(chrom), transcript["t_strand"],
        str(tx_start_0b), str(tx_end_0b),
        str(cds_start_0b), str(cds_end_0b),
        str(len(starts_0b)),
        ",".join(str(x) for x in starts_0b) + ",",
        ",".join(str(x) for x in ends_0b) + ",",
        "0", gene, "cmpl", "cmpl",
        ",".join(str(x) for x in exon_frames) + ",",
        gene, name, "TRUE",
    ]
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        f.write("\t".join(row) + "\n")


def write_variant_vcf(path, variant_id, chrom, pos, ref, alt, transcript, scores):
    """Write a 1-variant VCF with a SpliceAI INFO field (--include format).

    INFO field matches the --include layout from the reference parser:
       ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL|
       DS_AG_REF|DS_AG_ALT|DS_AL_REF|DS_AL_ALT|DS_DG_REF|DS_DG_ALT|DS_DL_REF|DS_DL_ALT
    """
    name = transcript["t_id"]
    if transcript.get("t_refseq_ids"):
        name = transcript["t_refseq_ids"][0]
    symbol = f"RefSeqTx-{name}"

    # `or 0` guards against explicit JSON null in the API response --
    # dict.get(k, 0) only returns the default when the key is absent.
    def f2(k):
        return f"{float(scores.get(k) or 0):.2f}"

    def i0(k):
        return str(int(scores.get(k) or 0))

    info = "SpliceAI={}".format("|".join([
        alt, symbol,
        f2("DS_AG"), f2("DS_AL"), f2("DS_DG"), f2("DS_DL"),
        i0("DP_AG"), i0("DP_AL"), i0("DP_DG"), i0("DP_DL"),
        f2("DS_AG_REF"), f2("DS_AG_ALT"),
        f2("DS_AL_REF"), f2("DS_AL_ALT"),
        f2("DS_DG_REF"), f2("DS_DG_ALT"),
        f2("DS_DL_REF"), f2("DS_DL_ALT"),
    ]))
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.1\n")
        f.write(f"##contig=<ID={chrom}>\n")
        f.write(
            '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAI '
            'with REF/ALT raw scores (--include format)">\n'
        )
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write(f"{chrom}\t{pos}\t{variant_id}\t{ref}\t{alt}\t.\tPASS\t{info}\n")


# ---- Reference parser invocation ------------------------------------------
def run_ref_python(vcf, refseq, out, hg, fasta):
    cmd = [
        sys.executable, REF_PYTHON,
        "-i", vcf, "-r", refseq, "-o", out,
        "--ref", REF_NAME[hg], "--fasta", fasta,
        "--include",
    ]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=180)


def run_ref_r(vcf, refseq, out, hg):
    cmd = [
        RSCRIPT, REF_R,
        "-i", vcf, "-r", refseq, "-o", out,
        "--ref", REF_NAME[hg], "--include",
    ]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=300)


def parse_ref_output(path):
    """Parse reference parser's TSV output into a normalized dict.

    Returns None if the file is missing or has no data rows.
    """
    if not os.path.exists(path):
        return None
    with open(path) as f:
        lines = f.read().splitlines()
    if len(lines) < 2:
        return None
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))

    def _na(v):
        return None if v in ("NA", "-", "", "nan", "NaN") else v

    type_map = [
        ("Exon_skipping",            "exon_skipping"),
        ("Intron_retention",         "whole_intron_retention"),
        ("Pseudoexon_activation",    "pseudoexon"),
        ("Increased_exon_inclusion", "increased_exon_inclusion"),
        ("Partial_exon_deletion",    "partial_exon_deletion"),
        ("Partial_intron_retention", "partial_intron_retention"),
    ]
    aberration_types = {our for ref, our in type_map if row.get(ref) == "YES"}

    def _fs(val):
        if val == "YES":
            return True
        if val == "NO":
            return False
        return None

    frameshifts = {
        "exon_skipping":             _fs(row.get("Exon_skipping_frameshift")),
        "whole_intron_retention":    _fs(row.get("Intron_retention_frameshift")),
        "pseudoexon":                _fs(row.get("Pseudoexon_frameshift")),
        # Ref parsers don't emit an Increased_exon_inclusion_frameshift column
        # (they report IEI frameshift in the aaseq column as text only), so we
        # can't frameshift-compare IEI. Leave as None and let the comparison
        # see matching None-vs-None (no divergence). Any new-impl frameshift
        # value for IEI will surface as a soft (bool, None) divergence.
        "increased_exon_inclusion":  None,
        # Reference parser has a single Partial_frameshift shared by the two
        # partial branches.
        "partial_exon_deletion":     _fs(row.get("Partial_frameshift")),
        "partial_intron_retention":  _fs(row.get("Partial_frameshift")),
    }

    def _to_int(val):
        v = _na(val)
        if v is None:
            return None
        try:
            return int(v)
        except ValueError:
            return None

    # Lost_exons is "N" or "N,M" (never a range). Normalize to a sorted
    # tuple of ints so it compares equal to the new impl's parsed tuple.
    def _lost_exons(val):
        v = _na(val)
        if v is None:
            return None
        try:
            return tuple(sorted(int(x) for x in v.split(",") if x.strip()))
        except ValueError:
            return None

    regions = {
        "exon_skipping":             _lost_exons(row.get("Lost_exons")),
        # whole_intron_retention: new impl doesn't structurally expose the
        # retained-intron number, so the region comparison is skipped
        # regardless of what ref reports. Keep the value for debugging only.
        "whole_intron_retention":    _to_int(row.get("Retained_intron")),
        "pseudoexon":                _to_int(row.get("Pseudoexon_intron")),
        "increased_exon_inclusion":  _to_int(row.get("Exon_with_increased_inclusion")),
        # Partial events: reference reports genomic start/end, not a region number.
        "partial_exon_deletion":     None,
        "partial_intron_retention":  None,
    }

    return {
        "any_aberration": row.get("Any_splicing_aberration") == "YES",
        "aberration_types": aberration_types,
        "frameshifts": frameshifts,
        "regions": regions,
        "raw_row": row,
    }


# ---- New-impl invocation --------------------------------------------------
def run_new_impl(variant_pos, transcript, scores):
    """Call sai10k_compute_predictions directly.

    Assembles the transcript_scores dict expected by the function from
    the API response (top transcript) + its delta-score fields.
    """
    # `or 0` guards against explicit JSON null in the API response --
    # dict.get(k, 0) only returns the default when the key is absent.
    ts = {
        "DS_AG": float(scores.get("DS_AG") or 0),
        "DS_AL": float(scores.get("DS_AL") or 0),
        "DS_DG": float(scores.get("DS_DG") or 0),
        "DS_DL": float(scores.get("DS_DL") or 0),
        "DP_AG": int(scores.get("DP_AG") or 0),
        "DP_AL": int(scores.get("DP_AL") or 0),
        "DP_DG": int(scores.get("DP_DG") or 0),
        "DP_DL": int(scores.get("DP_DL") or 0),
        "DS_AG_ALT": float(scores.get("DS_AG_ALT") or 0),
        "DS_AL_ALT": float(scores.get("DS_AL_ALT") or 0),
        "DS_DG_ALT": float(scores.get("DS_DG_ALT") or 0),
        "DS_DL_ALT": float(scores.get("DS_DL_ALT") or 0),
        "EXON_STARTS": transcript["EXON_STARTS"],
        "EXON_ENDS": transcript["EXON_ENDS"],
        "CDS_START": transcript.get("CDS_START"),
        "CDS_END": transcript.get("CDS_END"),
        "STRAND": transcript["t_strand"],
        "t_strand": transcript["t_strand"],
        "t_id": transcript["t_id"],
        "g_name": transcript.get("g_name"),
        "t_priority": transcript.get("t_priority"),
    }
    predictions = sai10k_compute_predictions(ts, variant_pos)

    aberration_types = {a["aberration_type"] for a in predictions["aberrations"]}

    # The new impl can emit multiple aberrations with the same type for a
    # single variant (e.g., two partial_exon_deletion from both gain_B_acceptor
    # and gain_B_donor -- sai10k_predictions.py:554-580). A naive
    # dict[type]=value would drop all but the last. Aggregate instead:
    #   frameshift: collapse to True if ANY is True, False if ALL are False,
    #               else None (matches the ref parser's single Partial_frameshift
    #               column which also flips on any true).
    #   region: sorted tuple of distinct values (None excluded).
    frameshifts_list = {}
    regions_list = {}
    # Track aberrations where the new impl emitted its explicit
    # "could not be mapped" fallback -- distinguishing benign fallbacks from
    # regex mismatches on unanticipated frameshift_description formats.
    unmappable = set()
    for a in predictions["aberrations"]:
        ab_type = a["aberration_type"]
        frameshifts_list.setdefault(ab_type, []).append(a.get("frameshift"))
        regions_list.setdefault(ab_type, []).append(_region_for_comparison(a, ab_type))
        desc = a.get("frameshift_description") or ""
        if "could not be mapped" in desc:
            unmappable.add(ab_type)

    frameshifts = {}
    for ab_type, vals in frameshifts_list.items():
        if any(v is True for v in vals):
            frameshifts[ab_type] = True
        elif all(v is False for v in vals):
            frameshifts[ab_type] = False
        else:
            frameshifts[ab_type] = None

    regions = {}
    for ab_type, vals in regions_list.items():
        distinct = [v for v in vals if v is not None]
        if not distinct:
            regions[ab_type] = None
        elif len(distinct) == 1:
            regions[ab_type] = distinct[0]
        else:
            regions[ab_type] = tuple(sorted(set(distinct), key=str))

    return {
        "any_aberration": bool(aberration_types),
        "aberration_types": aberration_types,
        "frameshifts": frameshifts,
        "regions": regions,
        "unmappable": unmappable,
        "raw_aberrations": predictions["aberrations"],
    }


# Matches the new impl's frameshift_description formats for exon_skipping, e.g.
#   "Exon 11 skipping (123bp coding seq.) - frameshift"
#   "Exons 2, 3 skipping (..."
#   "Potential exon 11 skipping (..."
#   "Potential exons 2, 3 skipping (..."
_SKIPPED_EXONS_RE = re.compile(
    r"(?:Potential\s+)?[Ee]xons?\s+([\d,\s]+?)\s+skipping",
)


def _region_for_comparison(aberration, ab_type):
    """Return the region identifier the new impl exposes, in a form comparable
    to the reference parser's TSV output.

    For exon_skipping, the new impl doesn't expose the skipped-exon numbers
    on the aberration dict structurally -- it embeds them in
    frameshift_description. Parse them out as a sorted tuple of ints.

    For whole_intron_retention, the new impl's affected_region is the
    VARIANT's containing region (intron N for intronic variants, but for
    exonic-boundary variants it's an exon). That's not the same semantic as
    the reference's Retained_intron, so return None and the comparison skips.

    For pseudoexon / increased_exon_inclusion, the new impl's affected_region
    coincides semantically with the reference fields (variant is in / near
    the affected exon / host intron).
    """
    if ab_type == "exon_skipping":
        desc = aberration.get("frameshift_description") or ""
        match = _SKIPPED_EXONS_RE.search(desc)
        if not match:
            return None
        nums = [int(x) for x in match.group(1).split(",") if x.strip()]
        if not nums:
            return None
        return tuple(sorted(nums))

    if ab_type == "whole_intron_retention":
        # No structural retained-intron field on the aberration dict; skip.
        return None

    region = aberration.get("affected_region") or {}
    return region.get("region_number")


# ---- Full per-variant pipeline --------------------------------------------
def prepare_variant(variant, hg):
    """Resolve HGVS, query API, pick top transcript. Returns a dict or None."""
    actual_variant = variant
    if variant.startswith("HGVS:"):
        actual_variant = normalize_hgvs(variant[len("HGVS:"):], hg)

    response = query_spliceai_api(actual_variant, hg)
    if response.get("error"):
        return {"skip_reason": f"API error: {response['error']}"}

    scores_list = response.get("scores") or []
    if not scores_list:
        return {"skip_reason": "No transcripts returned by API"}

    top = sai10k_select_transcript(scores_list)
    if not top:
        return {"skip_reason": "No top transcript selected"}

    chrom, pos, ref, alt = parse_variant_coords(actual_variant)

    # Sanity: API returns t_strand for Ensembl transcripts but some records
    # may omit structure fields when the variant is outside transcribed regions.
    if not top.get("EXON_STARTS") or not top.get("EXON_ENDS"):
        return {"skip_reason": "Top transcript has no exon structure"}
    if top.get("t_strand") not in ("+", "-"):
        return {"skip_reason": f"Top transcript has invalid t_strand: {top.get('t_strand')!r}"}

    return {
        "actual_variant": actual_variant,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "top_transcript": top,
        "is_indel": _is_indel(actual_variant),
    }


def compare_variant(variant, hg):
    """Run all three implementations for one (variant, hg) and return a diff dict."""
    prep = prepare_variant(variant, hg)
    if prep.get("skip_reason"):
        return {"status": "skipped", "reason": prep["skip_reason"]}

    top = prep["top_transcript"]
    scores = top  # (the API nests scores directly into the transcript dict)

    new_impl_result = run_new_impl(prep["pos"], top, scores)

    with tempfile.TemporaryDirectory(prefix="sai10k_ref_") as tmp:
        vcf = os.path.join(tmp, "var.vcf")
        refseq = os.path.join(tmp, "refseq.tsv")
        out_py = os.path.join(tmp, "out_py.tsv")
        out_r = os.path.join(tmp, "out_r.tsv")

        variant_id = f"v_{prep['actual_variant']}"
        write_variant_vcf(vcf, variant_id, prep["chrom"], prep["pos"],
                          prep["ref"], prep["alt"], top, scores)
        write_refseq_table(refseq, prep["chrom"], top)

        # Reference Python
        fasta = _find_fasta(hg)
        ref_py_result = None
        ref_py_err = None
        if REF_PYTHON is None:
            ref_py_err = f"Reference python parser not found. Searched: {REF_DIR_CANDIDATES}"
        elif fasta is None:
            ref_py_err = f"No FASTA available for hg{hg}"
        else:
            try:
                proc = run_ref_python(vcf, refseq, out_py, hg, fasta)
            except subprocess.TimeoutExpired:
                ref_py_err = "Python ref timed out"
            except (OSError, subprocess.SubprocessError) as e:
                ref_py_err = f"Python ref failed to run: {type(e).__name__}: {e}"
            else:
                if proc.returncode != 0:
                    ref_py_err = f"Python ref exit {proc.returncode}: {proc.stderr[-500:]}"
                else:
                    ref_py_result = parse_ref_output(out_py)
                    if ref_py_result is None:
                        ref_py_err = "Python ref produced no output rows"

        # Reference R. Runs regardless of whether Python ref succeeded so a
        # timeout/crash on one side doesn't silently drop the other.
        ref_r_result = None
        ref_r_err = None
        if not HAS_R:
            ref_r_err = "Rscript not installed on this machine"
        elif REF_R is None:
            ref_r_err = f"Reference R parser not found. Searched: {REF_DIR_CANDIDATES}"
        else:
            try:
                proc = run_ref_r(vcf, refseq, out_r, hg)
            except subprocess.TimeoutExpired:
                ref_r_err = "R ref timed out"
            except (OSError, subprocess.SubprocessError) as e:
                ref_r_err = f"R ref failed to run: {type(e).__name__}: {e}"
            else:
                if proc.returncode != 0:
                    ref_r_err = f"R ref exit {proc.returncode}: {proc.stderr[-500:]}"
                else:
                    ref_r_result = parse_ref_output(out_r)
                    if ref_r_result is None:
                        ref_r_err = "R ref produced no output rows"

    # -- Compare results --
    diff = {
        "status": "ok",
        "variant": prep["actual_variant"],
        "hg": hg,
        "is_indel": prep["is_indel"],
        "transcript": top.get("t_id"),
        "gene": top.get("g_name"),
        "new_impl": new_impl_result,
        "ref_py": ref_py_result,
        "ref_py_err": ref_py_err,
        "ref_r": ref_r_result,
        "ref_r_err": ref_r_err,
    }

    diff["discrepancies"] = _collect_discrepancies(new_impl_result, ref_py_result, ref_r_result,
                                                   is_indel=prep["is_indel"])
    return diff


# Region comparison is skipped for aberration types where the new impl and
# the reference parser don't share a common structural field. Mismatches here
# aren't meaningful; they're reported as soft-divergence info, not hard fails.
#
# * partial_exon_deletion / partial_intron_retention: ref reports genomic
#   start/end, not a region number.
# * whole_intron_retention: new impl doesn't structurally expose the
#   retained-intron number (only `variant's containing region`, which is
#   NOT the same thing for exonic-edge variants).
# * increased_exon_inclusion: ref reports the INCLUDED exon number; new
#   impl exposes the VARIANT's containing region. For splice-region variants
#   in the flanking intron these disagree by construction.
_SKIP_REGION_COMPARE = frozenset({
    "partial_exon_deletion",
    "partial_intron_retention",
    "whole_intron_retention",
    "increased_exon_inclusion",
})


# Frameshift comparison is skipped for aberration types where the ref parser
# doesn't emit a structural frameshift value. Ref reports IEI frameshift only
# as a free-text aaseq label, so we can't machine-compare it.
_SKIP_FRAMESHIFT_COMPARE = frozenset({
    "increased_exon_inclusion",
})


def _compare_one_ref(new, ref, ref_label, is_indel):
    """Compare the new impl's output against one reference impl (Python or R).

    Returns a list of human-readable discrepancy strings. Messages tagged
    "(KNOWN DIVERGENCE ...)" are soft divergences that don't fail tests.

    Indels: reference parsers NaN-out scores for indels, so ref['aberration_types']
    is always empty. If the new impl predicted aberrations, that's the expected
    case (not a failure) -- ref simply can't evaluate the indel.
    """
    out = []
    prefix = f"[{ref_label}]"

    if is_indel:
        if ref["aberration_types"]:
            # Unexpected: the ref parser usually drops indels at preprocessing.
            # If it did emit aberrations, surface them.
            out.append(
                f"{prefix} indel: unexpected aberrations from ref: "
                f"{sorted(ref['aberration_types'])} (KNOWN DIVERGENCE: indel)"
            )
        elif new["aberration_types"]:
            out.append(
                f"{prefix} indel: new impl predicted "
                f"{sorted(new['aberration_types'])}, ref skipped "
                "(KNOWN DIVERGENCE: indel)"
            )
        return out

    if new["aberration_types"] != ref["aberration_types"]:
        out.append(
            f"{prefix} aberration_types: "
            f"new={sorted(new['aberration_types'])} "
            f"ref={sorted(ref['aberration_types'])}"
        )

    for ab_type in new["aberration_types"] & ref["aberration_types"]:
        nr = new["regions"].get(ab_type)
        rr = ref["regions"].get(ab_type)

        if ab_type not in _SKIP_REGION_COMPARE:
            # For exon_skipping, treat None as a soft divergence ONLY when we
            # can attribute it to one side's explicit "could not map" fallback:
            #   * new impl: frameshift_description literally contains
            #     "could not be mapped" (tracked via new["unmappable"])
            #   * ref: Lost_exons column is NA (rr is None)
            # A None that doesn't match either condition means the regex
            # failed on an unexpected frameshift_description format -- which
            # is a test bug, not a benign divergence, so it should fail hard.
            new_is_fallback = ab_type in new.get("unmappable", set())
            if nr is None or rr is None:
                if nr != rr:
                    if ab_type == "exon_skipping" and (new_is_fallback or rr is None):
                        out.append(
                            f"{prefix} {ab_type} region: new={nr} ref={rr} "
                            "(KNOWN DIVERGENCE: one side could not map)"
                        )
                    else:
                        out.append(f"{prefix} {ab_type} region: new={nr} ref={rr}")
            elif nr != rr:
                out.append(f"{prefix} {ab_type} region: new={nr} ref={rr}")

        if ab_type not in _SKIP_FRAMESHIFT_COMPARE:
            new_fs = new["frameshifts"].get(ab_type)
            ref_fs = ref["frameshifts"].get(ab_type)
            # Python's `None != None` is False, so both-None is silently OK.
            # (bool, None) pairs DO differ and are surfaced as soft divergences
            # (ref parser doesn't compute frameshift for non-coding transcripts;
            # new impl applies its own CDS-based predicate).
            if new_fs != ref_fs:
                out.append(
                    f"{prefix} {ab_type} frameshift: new={new_fs} ref={ref_fs} "
                    "(KNOWN DIVERGENCE: CDS vs full-size)"
                )

    return out


def _collect_discrepancies(new, ref_py, ref_r, is_indel):
    """List human-readable mismatches across the three impls.

    Soft divergences are tagged "(KNOWN DIVERGENCE: ...)" and don't fail
    tests. Hard discrepancies (aberration type set mismatch, region mismatch
    for comparable types) fail.
    """
    out = []
    if ref_py is not None:
        out.extend(_compare_one_ref(new, ref_py, "ref_py", is_indel))
    if ref_r is not None:
        out.extend(_compare_one_ref(new, ref_r, "ref_r", is_indel))
    return out


def _expand_cases():
    """Expand VARIANT_SOURCES into [(variant, hg, label), ...] dedup'd."""
    seen = set()
    out = []
    for v, hgs, label in VARIANT_SOURCES:
        for hg in hgs:
            key = (v, hg)
            if key in seen:
                continue
            seen.add(key)
            out.append((v, hg, label))
    return out


TEST_CASES = _expand_cases()


# ---- unittest integration -------------------------------------------------
class TestSAI10kReferenceConsistency(unittest.TestCase):

    def _run_case(self, variant, hg):
        diff = compare_variant(variant, hg)
        if diff["status"] == "skipped":
            self.skipTest(diff["reason"])
        # If BOTH reference parsers failed to run, there's nothing to compare
        # against; skip rather than silently pass with an empty discrepancy list.
        if diff["ref_py"] is None and diff["ref_r"] is None:
            reasons = [r for r in (diff["ref_py_err"], diff["ref_r_err"]) if r]
            self.skipTest(
                "Both reference parsers unavailable: " + " | ".join(reasons)
            )
        # KNOWN DIVERGENCE messages are soft: CDS-vs-full-size frameshift,
        # indel handling, and region fields the new impl doesn't structurally
        # expose. Hard discrepancies (type-set mismatch, comparable regions)
        # fail the test.
        hard_discrepancies = [d for d in diff["discrepancies"]
                              if "KNOWN DIVERGENCE" not in d]
        if hard_discrepancies:
            msg = (
                f"\nVariant: {diff['variant']} (hg{hg}) "
                f"gene={diff['gene']} tx={diff['transcript']}\n  "
                + "\n  ".join(hard_discrepancies)
            )
            self.fail(msg)


def _make_test(variant, hg):
    def test_fn(self):
        self._run_case(variant, hg)
    safe = variant.replace("-", "_").replace(":", "_").replace(".", "_").replace(">", "_")
    test_fn.__name__ = f"test_{safe}_hg{hg}"
    return test_fn


for _v, _hg, _label in TEST_CASES:
    _fn = _make_test(_v, _hg)
    setattr(TestSAI10kReferenceConsistency, _fn.__name__, _fn)


# ---- CLI entrypoint -------------------------------------------------------
def dump_all(json_out=None):
    """Run every case and print a summary diff. Useful for triaging."""
    all_diffs = []
    fail = 0
    skip = 0
    ok = 0
    for i, (variant, hg, label) in enumerate(TEST_CASES, 1):
        print(f"\n[{i}/{len(TEST_CASES)}] {variant} (hg{hg}) -- {label}")
        try:
            diff = compare_variant(variant, hg)
        except Exception as e:
            print(f"  ERROR: {type(e).__name__}: {e}")
            all_diffs.append({"variant": variant, "hg": hg, "error": str(e)})
            fail += 1
            continue

        if diff["status"] == "skipped":
            print(f"  SKIP: {diff['reason']}")
            skip += 1
            all_diffs.append({"variant": variant, "hg": hg, "skipped": diff["reason"]})
            continue

        print(f"  gene={diff['gene']} tx={diff['transcript']}")
        print(f"  new_impl aberrations: {sorted(diff['new_impl']['aberration_types'])}")
        if diff["ref_py"] is not None:
            print(f"  ref_py   aberrations: {sorted(diff['ref_py']['aberration_types'])}")
        elif diff["ref_py_err"]:
            print(f"  ref_py   ERROR: {diff['ref_py_err']}")
        if diff["ref_r"] is not None:
            print(f"  ref_r    aberrations: {sorted(diff['ref_r']['aberration_types'])}")
        elif diff["ref_r_err"]:
            print(f"  ref_r    ERROR: {diff['ref_r_err']}")

        hard = [d for d in diff["discrepancies"] if "KNOWN DIVERGENCE" not in d]
        soft = [d for d in diff["discrepancies"] if "KNOWN DIVERGENCE" in d]
        if hard:
            for d in hard:
                print(f"  DIFF: {d}")
            fail += 1
        else:
            if soft:
                for d in soft:
                    print(f"  (known) {d}")
            ok += 1

        all_diffs.append({
            "variant": variant, "hg": hg, "gene": diff["gene"],
            "transcript": diff["transcript"],
            "new_impl_types": sorted(diff["new_impl"]["aberration_types"]),
            "ref_py_types": (sorted(diff["ref_py"]["aberration_types"])
                             if diff["ref_py"] else None),
            "ref_r_types": (sorted(diff["ref_r"]["aberration_types"])
                            if diff["ref_r"] else None),
            "ref_py_err": diff["ref_py_err"],
            "ref_r_err": diff["ref_r_err"],
            "discrepancies": diff["discrepancies"],
        })

    print(f"\n=== Summary: {ok} ok | {fail} mismatched | {skip} skipped ===")

    if json_out:
        with open(json_out, "w") as f:
            json.dump(all_diffs, f, indent=2, sort_keys=True)
        print(f"Wrote {json_out}")

    return fail


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dump-only", action="store_true",
                        help="Print per-variant diffs and exit (no unittest).")
    parser.add_argument("--json", help="Write dump results to this JSON file.")
    args, remaining = parser.parse_known_args()

    if args.dump_only:
        sys.exit(1 if dump_all(args.json) else 0)
    else:
        unittest.main(argv=[sys.argv[0]] + remaining)
