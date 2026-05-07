"""
SAI-10k-calc predictions module.

Implements the SAI-10k-calc algorithm (Canson et al. 2023, Bioinformatics 39(4),
btad179) to classify splicing aberrations from SpliceAI outputs. Terminal paper
classes produced by this module:

    - exon_skipping
    - whole_intron_retention
    - pseudoexon (pseudoexon activation)
    - increased_exon_inclusion
    - partial_exon_deletion
    - partial_intron_retention

When a variant doesn't match any paper terminal class, an empty aberrations
list is returned (equivalent to Any_splicing_aberration=NO in the reference
parser). Non-paper placeholder labels (donor_shift, *_loss, *_gain, unknown)
are no longer emitted.

The flowchart's three branches evaluate independently so a single variant can
produce combinations (e.g. exon_skipping + partial_intron_retention per
Canson Figure 1 D/E).

Divergence from the reference parser - "Potential" loss labels:
    The reference parser (spliceAI_parser.R / spliceai_parser.py) emits a
    binary PASS/FAIL on the loss pair using fixed MIN/MAX thresholds. This
    module additionally flags loss calls where min(DS_AL, DS_DL) is above
    DS_ALDL_MIN but below DS_ALDL_MAX as "Potential" (e.g. "Potential exon
    skipping", "Potential whole intron retention") in the user-visible
    `frameshift_description`. The aberration_type field itself is unchanged.
    This third tier was added at reviewer request (2026-04-17) to surface
    asymmetric loss-pair calls where one delta score is weak; it is not
    present in the reference implementation.
"""

import time


# Reference-parser defaults (spliceai_parser.py: --DS_ALDL_MIN_T,
# --DS_ALDL_MAX_T, --DS_AGDG_MIN_T, --DS_AGDG_MAX_T, all default 0.02 / 0.2).
#
# NOTE: paper main text (btad179 §2) and Supplementary File 1 state the
# DS_AGDG_MAX default as 0.05 for pseudoexon gain. The reference Python/R
# implementation instead defaults both MAX thresholds to 0.2. This module
# follows the reference implementation so users see the same calls they'd
# get running SAI-10k-calc locally. To match the paper text, change
# DS_AGDG_MAX to 0.05.
TRANSCRIPT_PRIORITY_ORDER = {
    "MS": 3,  # MANE select transcript
    "MP": 2,  # MANE plus clinical transcript
    "C": 1,   # canonical transcript
    "N": 0,
}

# Loss pair (exon skipping / whole intron retention):
DS_ALDL_MIN = 0.02
DS_ALDL_MAX = 0.2
# Gain pair (pseudoexon activation / increased exon inclusion):
DS_AGDG_MIN = 0.02
DS_AGDG_MAX = 0.2

# Single-score threshold used in the cryptic acceptor/donor subflow.
MIN_DELTA_SCORE = 0.2
# Raw alt-allele score threshold used by the cryptic subflow's conditions.
ALT_SCORE_STRONG = 0.9
# Loss-partner score threshold used in the cryptic subflow's conditions.
LOSS_PARTNER_MIN = 0.1
# Delta between ALT-allele scores required for cryptic activation.
ALT_DIFF_HIGH = 0.2
ALT_DIFF_LOW = -0.2

# Size constraints on predicted gained exons (pseudoexon).
GEX_SIZE_MIN = 25
GEX_SIZE_MAX = 500

# d_50bp in Canson et al.: loss-branch aberrations require the variant to be
# within 50 bp of an exon-intron junction (exonic variants trivially qualify).
# Pseudoexon activation requires the variant to be > 50 bp into an intron.
DISTANCE_EXON_MAX_LOSS = 50
DISTANCE_NATIVE_SITE_MIN_PSEUDOEXON = 50
# The pseudoexon branch also requires AG and DG to sit comfortably inside the
# host intron. spliceAI_parser.R lines 211-214 (find_pseudoexon_position) use
# directional checks against the host intron's two boundaries: on + strand,
# `intronStart + 50 < AGpos` (AG vs the 5' boundary) and
# `intronEndAdj - 50 > DGpos` (DG vs the 3' boundary); on - strand the roles
# are swapped (DG vs 5', AG vs 3'). `intronStart` is the first intronic base,
# so a distance measured from the flanking exon boundary (last exonic base)
# is one larger and the threshold becomes `> 51`.
DISTANCE_AG_DG_MIN_FROM_INTRON_BOUNDARY = 51
# Partial intron retention requires both the variant and the retained segment
# to be within 250 bp of an exon-intron junction.
DISTANCE_PARTIAL_RETENTION_MAX = 250



def _coerce_float(value):
    if value is None or value == '':
        return 0.0
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _coerce_int(value):
    if value is None or value == '':
        return 0
    try:
        return int(value)
    except (TypeError, ValueError):
        try:
            return int(float(value))
        except (TypeError, ValueError):
            return 0


def sai10k_find_affected_region(variant_pos, exon_starts, exon_ends, strand='+'):
    """
    Locate the variant in an exon, intron, or transcript-flanking region, with
    biological (transcript-order) numbering on reverse-strand genes.

    Returns dict with:
        region_type: 'exon', 'intron', or 'flank'
                     ('flank' = upstream of first exon or downstream of last;
                      matches the R reference's behavior of associating
                      out-of-bounds variants with the closest exon, R:851-874)
        region_number: 1-based biological exon/intron number
        _genomic_index: 0-based index into exon_starts/exon_ends (internal use)
        distance_to_boundary: bp to nearest exon-intron junction
                              (0 for exonic variants)
        nearest_boundary: 'donor' or 'acceptor' in transcript semantics
    """
    if not exon_starts or not exon_ends:
        return None

    num_exons = len(exon_starts)
    is_reverse = strand == '-'

    for i in range(num_exons):
        e_start = exon_starts[i]
        e_end = exon_ends[i]
        if e_start <= variant_pos <= e_end:
            dist_to_start = variant_pos - e_start
            dist_to_end = e_end - variant_pos
            exon_num = num_exons - i if is_reverse else i + 1
            if is_reverse:
                nearest = 'donor' if dist_to_start < dist_to_end else 'acceptor'
            else:
                nearest = 'acceptor' if dist_to_start < dist_to_end else 'donor'
            return {
                'region_type': 'exon',
                'region_number': exon_num,
                '_genomic_index': i,
                'distance_to_boundary': 0,  # variant is IN the exon
                'nearest_boundary': nearest,
            }

    for i in range(num_exons - 1):
        intron_start = exon_ends[i] + 1
        intron_end = exon_starts[i + 1] - 1
        if intron_start <= variant_pos <= intron_end:
            # Distance to nearest EXON boundary (not intron boundary).
            # Reference parser (R:857-871) computes d_from_exon as
            # POS - eEnd / eStart - POS. Using exon_ends[i] and
            # exon_starts[i+1] directly avoids the off-by-one that
            # measuring from intron_start (= exon_ends[i]+1) introduces.
            dist_to_upstream_exon = variant_pos - exon_ends[i]
            dist_to_downstream_exon = exon_starts[i + 1] - variant_pos
            intron_num = num_exons - 1 - i if is_reverse else i + 1
            # Tie-break with `<=` so the upstream genomic exon wins on an
            # exact midpoint, matching `_get_native_splice_sites`'s assoc_idx
            # tie-break (`assoc_idx = idx if dist_to_upstream <= dist_to_downstream`).
            # Without this, the loss branch's `nearest_boundary`-derived exon
            # and the gain-B branch's `_native.assoc_idx` could diverge for
            # an equidistant intronic variant.
            if is_reverse:
                nearest = 'acceptor' if dist_to_upstream_exon <= dist_to_downstream_exon else 'donor'
            else:
                nearest = 'donor' if dist_to_upstream_exon <= dist_to_downstream_exon else 'acceptor'
            return {
                'region_type': 'intron',
                'region_number': intron_num,
                '_genomic_index': i,
                'distance_to_boundary': min(dist_to_upstream_exon, dist_to_downstream_exon),
                'nearest_boundary': nearest,
            }

    # Variant is outside the transcript span (5' or 3' flank). The R reference
    # (spliceAI_parser.R:851-874) computes d_from_exon against the closest exon
    # regardless of containment, so loss / pseudoexon / partial gates can still
    # fire near terminal exons. Associate the variant with the nearest exon.
    if variant_pos < exon_starts[0]:
        idx = 0
        dist = exon_starts[0] - variant_pos
        # Variant is on the genomic-low side of exon 0. On + strand this is
        # 5'-flank: variant precedes the exon's acceptor (at exon_starts[0]).
        # On - strand this is 3'-flank: variant lies past the exon's biological
        # donor (also at exon_starts[0]).
        nearest = 'donor' if is_reverse else 'acceptor'
        exon_num = num_exons if is_reverse else 1
    else:
        idx = num_exons - 1
        dist = variant_pos - exon_ends[idx]
        # Genomic-high side of last exon. On + strand: 3'-flank past donor at
        # exon_ends[-1]. On - strand: 5'-flank before biological acceptor at
        # exon_ends[-1].
        nearest = 'acceptor' if is_reverse else 'donor'
        exon_num = 1 if is_reverse else num_exons
    return {
        'region_type': 'flank',
        'region_number': exon_num,
        '_genomic_index': idx,
        'distance_to_boundary': dist,
        'nearest_boundary': nearest,
    }


def _get_native_splice_sites(affected_region, variant_pos, exon_starts, exon_ends, strand):
    """
    Return the native splice-site genomic positions and their delta-positions
    (offset from the variant) for the exon "associated" with the variant.

    For an exonic variant, the associated exon is the containing exon.
    For an intronic variant, the associated exon is the CLOSEST flanking exon
    (smaller distance_to_boundary), matching the reference parser's behavior
    of filtering to the single closest exon (spliceAI_parser.R:865-867,
    spliceai_parser.py:1218 — both keep `dist_exon_closest_abs == min(...)`).
    The actual difference from the reference is containment-based selection
    (this module locates the variant inside an exon or internal intron) vs.
    the reference's distance-only selection over all exons.

    Both acceptor and donor come from the SAME exon, so the increased_exon_inclusion
    condition (DP_AG == DP_NA AND DP_DG == DP_ND) is geometrically satisfiable.

    Also returns the adjacent exons' boundaries (prev_eend, next_estart)
    in genomic order, used by the cryptic subflow's orientation checks.

    Returns dict with:
        geo_na, geo_nd: 1-based genomic positions of native acceptor / donor
                        of the associated exon
        dp_na, dp_nd:   offsets from variant_pos
        native_exon_size: size of the associated native exon (bp)
        assoc_idx:      0-based genomic index of the associated exon
        prev_eend, next_estart: adjacent exon boundaries
                        (genomic order), or None at the transcript ends
    """
    if not affected_region or not exon_starts:
        return None

    idx = affected_region['_genomic_index']
    is_reverse = strand == '-'
    num_exons = len(exon_starts)

    region_type = affected_region['region_type']
    if region_type == 'exon' or region_type == 'flank':
        # 'flank' variants (outside the transcript span) are already associated
        # with the single closest exon in sai10k_find_affected_region.
        assoc_idx = idx
    else:
        # Intronic: pick the closer flanking exon. Measure distance from each
        # exon's nearest boundary (the last exonic base for the upstream exon,
        # the first exonic base for the downstream exon), consistent with
        # sai10k_find_affected_region:142-143 and the Canson R parser's
        # d_from_exon formulation (R:857-871). Prior versions measured from
        # the first intronic base, introducing a 1-bp upstream bias at ties.
        dist_to_upstream = variant_pos - exon_ends[idx]
        dist_to_downstream = exon_starts[idx + 1] - variant_pos
        assoc_idx = idx if dist_to_upstream <= dist_to_downstream else idx + 1

    exon_start = exon_starts[assoc_idx]
    exon_end = exon_ends[assoc_idx]
    # On '+' strand: acceptor at exon_start (5' end of exon), donor at exon_end (3' end).
    # On '-' strand: roles flip (biological 5' end is at exon_end, 3' end at exon_start).
    if is_reverse:
        geo_na = exon_end
        geo_nd = exon_start
    else:
        geo_na = exon_start
        geo_nd = exon_end

    prev_eend = exon_ends[assoc_idx - 1] if assoc_idx > 0 else None
    next_estart = exon_starts[assoc_idx + 1] if assoc_idx + 1 < num_exons else None

    return {
        'geo_na': geo_na,
        'geo_nd': geo_nd,
        'dp_na': geo_na - variant_pos,
        'dp_nd': geo_nd - variant_pos,
        'native_exon_size': exon_end - exon_start + 1,
        'assoc_idx': assoc_idx,
        'prev_eend': prev_eend,
        'next_estart': next_estart,
    }


def _find_exons_between_loss_positions(geo_al, geo_dl, strand, exon_starts, exon_ends):
    """
    Given the predicted acceptor-loss and donor-loss genomic positions, find
    the exon(s) whose native splice sites they match, and compute the total
    lost exon size + biological exon numbers for display.

    The reference parser's `find_exon_lost_gain()` matches GEO_AL to an
    exon_start (on '+' strand) or exon_end (on '-' strand), and GEO_DL to the
    opposite boundary. For multi-exon skipping, GEO_AL is at a more-5' exon's
    acceptor and GEO_DL is at a more-3' exon's donor (in transcript order).

    Returns dict with:
        genomic_indices: list of 0-based genomic exon indices bracketed by the
                         losses, or [] if no match could be made
        biological_numbers: list of 1-based biological exon numbers (sorted in
                            transcript order) for display
        total_size: sum of exon sizes (bp)
        match: True if at least one exon was identified
    """
    if not exon_starts or not exon_ends:
        return {'genomic_indices': [], 'biological_numbers': [], 'total_size': 0, 'match': False}

    num_exons = len(exon_starts)
    is_reverse = strand == '-'
    # Exact match (reference parser behavior). SpliceAI returns splice-site
    # positions exactly at exon boundaries; a wider tolerance would risk
    # matching the wrong exon in compact gene structures. Cases where DP_AL
    # or DP_DL lands interior to the affected exon — or only one of the two
    # lands at a native site, with the other on a cryptic position — fall
    # through to the sai10k_annotate_frameshift fallback (single-exon
    # skipping of the affected exon, for both exonic AND intronic variants).
    # On '+' strand: acceptor_loss is at exon_start (1-based closed), donor_loss
    # at exon_end. On '-' strand: acceptor_loss is at exon_end, donor_loss at
    # exon_start.
    acceptor_candidates = []
    donor_candidates = []
    for i in range(num_exons):
        if is_reverse:
            acceptor_pos = exon_ends[i]
            donor_pos = exon_starts[i]
        else:
            acceptor_pos = exon_starts[i]
            donor_pos = exon_ends[i]
        if acceptor_pos == geo_al:
            acceptor_candidates.append(i)
        if donor_pos == geo_dl:
            donor_candidates.append(i)

    if not acceptor_candidates or not donor_candidates:
        return {'genomic_indices': [], 'biological_numbers': [], 'total_size': 0, 'match': False}

    # Pick the pair that makes biological sense. In transcript order, the
    # acceptor should bound the 5' end of the skipped block and the donor
    # should bound the 3' end; single-exon skipping has acceptor_idx == donor_idx.
    # Try to find the tightest bracket.
    best = None
    for a_idx in acceptor_candidates:
        for d_idx in donor_candidates:
            lo = min(a_idx, d_idx)
            hi = max(a_idx, d_idx)
            span = hi - lo
            if best is None or span < best[0]:
                best = (span, lo, hi)

    if best is None:
        return {'genomic_indices': [], 'biological_numbers': [], 'total_size': 0, 'match': False}

    _, lo, hi = best
    # Divergence from reference parser: this implementation enumerates EVERY
    # exon between the bracketing acceptor and donor (range(lo, hi + 1)) and
    # sums their sizes for cds_size / frameshift / displayed exon list. The
    # reference's find_exon_lost_gain (spliceAI_parser.R:227-246,
    # spliceai_parser.py:293-329) instead applies a disjunctive boundary
    # filter (eStartAdj == ALAGposition | eEnd == DLDGposition), which only
    # matches the two boundary exons. For 3+-exon skips this means the
    # reference shows e.g. "2,4" while this module shows "2,3,4", and the
    # two parsers can disagree on the frameshift call when the omitted
    # middle exons' summed size mod 3 != 0. The full-bracket sum is
    # biologically more correct (the entire skipped block is removed from
    # the mature transcript); this divergence is intentional.
    genomic_indices = list(range(lo, hi + 1))
    biological_numbers = [num_exons - i if is_reverse else i + 1 for i in genomic_indices]
    biological_numbers.sort()
    total_size = sum(exon_ends[i] - exon_starts[i] + 1 for i in genomic_indices)

    return {
        'genomic_indices': genomic_indices,
        'biological_numbers': biological_numbers,
        'total_size': total_size,
        'match': True,
    }


def _whole_intron_retention_native_match(geo_al, geo_dl, strand, exon_starts, exon_ends):
    """
    Verify that geo_dl matches a native donor and geo_al matches a native
    acceptor of a SINGLE intron in the reference transcript (i.e. they bound
    the same intron between adjacent exons i and i+1). Returns the matched
    intron's biological (1-based, transcript-orientation) number, or None if
    no match.

    Without this check, the loss branch can emit "(potential) intron N
    retention" with a misleading size when SpliceAI's max_distance window is
    too narrow to capture both native sites and the weaker partner peak lands
    on a cryptic position. Example: BRCA2 c.-40+5G>T at max_distance=500
    yields DS_AL=0.03 with DP_AL pointing to a cryptic site 269bp downstream
    rather than the native acceptor of exon 2 (754bp downstream).
    """
    if len(exon_starts) < 2 or len(exon_ends) < 2:
        return None
    is_reverse = strand == '-'
    num_exons = len(exon_starts)
    for i in range(num_exons - 1):
        if is_reverse:
            native_donor = exon_starts[i + 1]
            native_acceptor = exon_ends[i]
        else:
            native_donor = exon_ends[i]
            native_acceptor = exon_starts[i + 1]
        if geo_dl == native_donor and geo_al == native_acceptor:
            return num_exons - 1 - i if is_reverse else i + 1
    return None


def sai10k_determine_aberrations(
        ds_ag, ds_al, ds_dg, ds_dl,
        dp_ag, dp_al, dp_dg, dp_dl,
        ds_ag_alt, ds_al_alt, ds_dg_alt, ds_dl_alt,
        variant_pos, exon_starts, exon_ends, strand,
        has_alt_scores=True,
):
    """
    Classify splice aberration(s) per Canson et al. 2023 flowchart v1.1.

    Returns a list of aberration dicts (possibly empty). Each dict contains:
        aberration_type: one of the six paper terminal classes
            (exon_skipping, whole_intron_retention, pseudoexon,
             increased_exon_inclusion, partial_exon_deletion,
             partial_intron_retention).
        affected_region: dict from sai10k_find_affected_region, or None.
        geo_ag, geo_al, geo_dg, geo_dl: 1-based genomic positions of predicted
            acceptor-gain, acceptor-loss, donor-gain, donor-loss sites.
        _branch: internal tag ('loss', 'gain_A', 'gain_B_acceptor', 'gain_B_donor')
        _native: dict with geo_na/geo_nd/dp_na/dp_nd/native_exon_size (or None)

    Args:
        has_alt_scores: True when the caller supplied real DS_*_ALT site scores
            (matches the reference parser's --include flag); False when only
            standard SpliceAI delta scores are available. With False, cryptic
            acceptor/donor activation falls back to the reference's simpler
            DS_AG/DS_DG dominance logic (spliceai_parser.py:1331-1338) instead
            of the ALT-score side-selection logic.

    The paper flowchart's three branches (LOSS, GAIN Subflow A, GAIN Subflow B)
    evaluate INDEPENDENTLY; a single variant can produce up to one aberration
    from each. The LOSS branch additionally requires whole_intron_retention to
    have geo_dl/geo_al match native splice sites bounding a single intron;
    otherwise the call is suppressed (see _whole_intron_retention_native_match).
    """
    ds_ag = _coerce_float(ds_ag)
    ds_al = _coerce_float(ds_al)
    ds_dg = _coerce_float(ds_dg)
    ds_dl = _coerce_float(ds_dl)
    ds_ag_alt = _coerce_float(ds_ag_alt)
    ds_al_alt = _coerce_float(ds_al_alt)
    ds_dg_alt = _coerce_float(ds_dg_alt)
    ds_dl_alt = _coerce_float(ds_dl_alt)
    dp_ag = _coerce_int(dp_ag)
    dp_al = _coerce_int(dp_al)
    dp_dg = _coerce_int(dp_dg)
    dp_dl = _coerce_int(dp_dl)

    max_score = max(ds_ag, ds_al, ds_dg, ds_dl)
    # Short-circuit only when every score is below the lowest paper threshold
    # (loss/gain paired MIN = 0.02). Below that, no paper branch can fire.
    if max_score < DS_ALDL_MIN:
        return []

    affected_region = sai10k_find_affected_region(variant_pos, exon_starts, exon_ends, strand)
    native = _get_native_splice_sites(affected_region, variant_pos, exon_starts, exon_ends, strand)

    geo_ag = variant_pos + dp_ag
    geo_al = variant_pos + dp_al
    geo_dg = variant_pos + dp_dg
    geo_dl = variant_pos + dp_dl

    strand_sign = -1 if strand == '-' else 1

    # d_from_exon: variant's distance from nearest exon-intron junction. 0 for
    # exonic variants; positive inside an intron.
    if affected_region and affected_region['region_type'] == 'exon':
        d_from_exon = 0
    elif affected_region:
        # Both 'intron' and 'flank' regions store distance-to-nearest-exon-edge.
        d_from_exon = affected_region['distance_to_boundary']
    else:
        d_from_exon = None  # couldn't locate variant

    def make(ab_type, branch):
        return {
            'aberration_type': ab_type,
            'affected_region': affected_region,
            'geo_ag': geo_ag,
            'geo_al': geo_al,
            'geo_dg': geo_dg,
            'geo_dl': geo_dl,
            '_branch': branch,
            '_native': native,
            '_strand': strand,
        }

    aberrations = []

    # -------------------------------------------------------------------
    # LOSS branch: exon_skipping / whole_intron_retention
    # Gated by paired min/max thresholds AND d_50bp (variant within 50 bp of
    # an exon-intron junction). The d_50bp gate suppresses deep-intronic
    # calls, consistent with spliceai_parser.py:1416-1421.
    # -------------------------------------------------------------------
    loss_paired_pass = (
        min(ds_dl, ds_al) >= DS_ALDL_MIN
        and max(ds_dl, ds_al) >= DS_ALDL_MAX
    )
    loss_d50_pass = d_from_exon is not None and d_from_exon <= DISTANCE_EXON_MAX_LOSS

    if loss_paired_pass and loss_d50_pass:
        # Canson reviewer feedback (2026-04-17): when the weaker score in the
        # loss pair is in [0.02, 0.2), flag the prediction as "Potential". The
        # paper's MIN threshold (0.02) was chosen for GWAS-scale sensitivity,
        # but high-impact variants can still land here (e.g. BRCA2 13:32363534 G>T,
        # DS_AL=0.13/DS_DL=0.93, experimentally 81% exon skipping). The
        # "Potential" prefix communicates the asymmetric pair to the reader.
        is_potential = min(ds_dl, ds_al) < DS_ALDL_MAX
        # Canson orientation gate:
        #   DP_AL*Strand < DP_DL*Strand  -> exon skipping
        #   else                          -> whole intron retention
        if dp_al * strand_sign < dp_dl * strand_sign:
            ab = make('exon_skipping', 'loss')
            ab['_potential'] = is_potential
            aberrations.append(ab)
        elif (matched_intron_num := _whole_intron_retention_native_match(
                geo_al, geo_dl, strand, exon_starts, exon_ends)) is not None:
            # Only emit whole_intron_retention when the predicted lost donor
            # and acceptor map to the native splice sites bounding a single
            # native intron. Suppresses misleading calls with a wrong intron
            # size when SpliceAI's max_distance window doesn't reach both
            # native sites (e.g. BRCA2 c.-40+5G>T at max_distance=500).
            ab = make('whole_intron_retention', 'loss')
            ab['_potential'] = is_potential
            ab['affected_intron_number'] = matched_intron_num
            aberrations.append(ab)

    # -------------------------------------------------------------------
    # GAIN Subflow A: pseudoexon / increased_exon_inclusion
    # -------------------------------------------------------------------
    gain_paired_pass = (
        min(ds_dg, ds_ag) >= DS_AGDG_MIN
        and max(ds_dg, ds_ag) >= DS_AGDG_MAX
    )
    orientation_pass = dp_ag * strand_sign < dp_dg * strand_sign

    if gain_paired_pass and orientation_pass:
        gex_size = (dp_dg - dp_ag) * strand_sign + 1

        # Identify which intron AG and DG fall into (same-intron check).
        def intron_of(position):
            for i in range(len(exon_starts) - 1):
                if exon_ends[i] < position < exon_starts[i + 1]:
                    return i
            return None

        ag_intron = intron_of(geo_ag)
        dg_intron = intron_of(geo_dg)

        # Directional distance check matching reference find_pseudoexon_position
        # (spliceai_parser.py:270-279): on + strand AG is checked against the 5'
        # (genomically-upstream) intron boundary and DG against the 3' boundary;
        # on - strand the roles are swapped. Only meaningful when AG and DG are
        # in the same intron.
        if ag_intron is not None and dg_intron is not None and ag_intron == dg_intron:
            host_5p_exon_end = exon_ends[ag_intron]          # last exonic base before intron
            host_3p_exon_start = exon_starts[ag_intron + 1]  # first exonic base after intron
            if strand_sign > 0:
                ag_dist_to_native = geo_ag - host_5p_exon_end
                dg_dist_to_native = host_3p_exon_start - geo_dg
            else:
                ag_dist_to_native = host_3p_exon_start - geo_ag
                dg_dist_to_native = geo_dg - host_5p_exon_end
        else:
            ag_dist_to_native = None
            dg_dist_to_native = None

        pseudoexon_conditions = (
            GEX_SIZE_MIN <= gex_size <= GEX_SIZE_MAX
            and d_from_exon is not None
            and d_from_exon > DISTANCE_NATIVE_SITE_MIN_PSEUDOEXON
            and ag_dist_to_native is not None
            and ag_dist_to_native > DISTANCE_AG_DG_MIN_FROM_INTRON_BOUNDARY
            and dg_dist_to_native is not None
            and dg_dist_to_native > DISTANCE_AG_DG_MIN_FROM_INTRON_BOUNDARY
        )

        if pseudoexon_conditions:
            aberrations.append(make('pseudoexon', 'gain_A'))
        elif native and native['native_exon_size'] and gex_size == native['native_exon_size']:
            # Increased exon inclusion: the gained exon exactly matches a
            # native exon AND the gain positions coincide with native sites.
            if dp_ag == native['dp_na'] and dp_dg == native['dp_nd']:
                aberrations.append(make('increased_exon_inclusion', 'gain_A'))

    # -------------------------------------------------------------------
    # GAIN Subflow B: cryptic acceptor / donor activation
    # Produces partial_exon_deletion or partial_intron_retention.
    # Uses raw alt-allele site scores (DS_*_ALT) per flowchart conditions.
    # Evaluates INDEPENDENTLY of Subflow A (flowchart: parallel subflows);
    # the reference parser likewise emits pseudoexon and partial calls
    # independently (spliceai_parser.py:1392-1463).
    # -------------------------------------------------------------------
    # Top-level entry gate from the reference parser: DS_AGDG_MAX > AGDG_T
    # (default AGDG_T = 0). Require at least one gain score to be positive.
    if max(ds_ag, ds_dg) > 0 and native is not None:
        # "Match" checks: does the predicted gain position exactly coincide
        # with the associated exon's native splice site? Reference parser
        # names these DG_native_match (donor-at-native-donor) and
        # DG_acceptor_match (acceptor-at-native-acceptor).
        dg_native_match = geo_dg == native['geo_nd']
        dg_acceptor_match = geo_ag == native['geo_na']

        if has_alt_scores:
            acceptor_diff = ds_ag_alt - ds_al_alt
            donor_diff = ds_dg_alt - ds_dl_alt

            # ALT-score conditions — two-branch OR, one for "cryptic overwhelms a
            # weak native" and one for "strong delta with no native suppression".
            cryptic_acceptor_check = (
                (ds_ag < MIN_DELTA_SCORE and ds_al > LOSS_PARTNER_MIN
                 and ds_ag_alt >= ALT_SCORE_STRONG and acceptor_diff >= ALT_DIFF_HIGH)
                or
                (ds_ag >= MIN_DELTA_SCORE and acceptor_diff > ALT_DIFF_LOW)
            )
            cryptic_donor_check = (
                (ds_dg < MIN_DELTA_SCORE and ds_dl > LOSS_PARTNER_MIN
                 and ds_dg_alt >= ALT_SCORE_STRONG and donor_diff >= ALT_DIFF_HIGH)
                or
                (ds_dg >= MIN_DELTA_SCORE and donor_diff > ALT_DIFF_LOW)
            )

            # Reference parser's side-selection logic (spliceai_parser.py:1312-1328):
            #   Cryptic acceptor fires iff:
            #       (DG_native_match=PASS AND DG_acceptor_match=FAIL)            # donor at native, acceptor is cryptic
            #     OR (DG_native_match=FAIL AND DG_acceptor_match=FAIL AND DS_AG > DS_DG)
            #   Cryptic donor fires iff:
            #       (DG_native_match=FAIL AND DG_acceptor_match=PASS)            # acceptor at native, donor is cryptic
            #     OR (DG_native_match=FAIL AND DG_acceptor_match=FAIL AND DS_AG < DS_DG)
            #   Neither fires when both positions match native (=,=) or when DS_AG == DS_DG at (≠,≠).
            activate_acceptor = cryptic_acceptor_check and (
                (dg_native_match and not dg_acceptor_match)
                or (not dg_native_match and not dg_acceptor_match and ds_ag > ds_dg)
            )
            activate_donor = cryptic_donor_check and (
                (not dg_native_match and dg_acceptor_match)
                or (not dg_native_match and not dg_acceptor_match and ds_ag < ds_dg)
            )
        else:
            # No DS_*_ALT scores supplied (standard SpliceAI output). Mirror the
            # reference parser's --include=False path (spliceai_parser.py:1331-1338):
            # simple dominance — a side activates only if its DS_*G score is over
            # threshold AND strictly greater than the other side's DS_*G. This
            # prevents the include-mode side-selection logic from firing donor
            # activation when DS_AG actually dominates DS_DG (or vice versa).
            activate_acceptor = ds_ag >= MIN_DELTA_SCORE and ds_ag > ds_dg
            activate_donor = ds_dg >= MIN_DELTA_SCORE and ds_dg > ds_ag

        dp_na = native['dp_na']
        dp_nd = native['dp_nd']

        if activate_acceptor and _cryptic_acceptor_orientation_ok(geo_ag, native, is_reverse=strand == '-'):
            # partial_size = (DP_AG*Strand) - (DP_NA*Strand)
            partial_size = (dp_ag - dp_na) * strand_sign
            if partial_size > 0 and d_from_exon is not None and d_from_exon <= DISTANCE_EXON_MAX_LOSS:
                aberrations.append(make('partial_exon_deletion', 'gain_B_acceptor'))
                aberrations[-1]['_partial_size'] = partial_size
                # When the cryptic site lands past the opposite native exon
                # boundary, the reference parser still emits Partial_exon_deletion
                # but its sequence step records "deletion greater than exon size"
                # (spliceai_parser.py:552-554). Flag here so the annotation step
                # surfaces the same outcome instead of computing a misleading
                # cds_size against a segment that overruns the exon.
                if partial_size >= native['native_exon_size']:
                    aberrations[-1]['_oversize_partial_deletion'] = True
            elif (partial_size < 0 and partial_size >= -DISTANCE_PARTIAL_RETENTION_MAX
                  and d_from_exon is not None and d_from_exon <= DISTANCE_PARTIAL_RETENTION_MAX):
                aberrations.append(make('partial_intron_retention', 'gain_B_acceptor'))
                # Store signed partial_size (negative = retention); formatters take abs() for display.
                aberrations[-1]['_partial_size'] = partial_size

        if activate_donor and _cryptic_donor_orientation_ok(geo_dg, native, is_reverse=strand == '-'):
            # partial_size = (DP_ND*Strand) - (DP_DG*Strand)
            partial_size = (dp_nd - dp_dg) * strand_sign
            if partial_size > 0 and d_from_exon is not None and d_from_exon <= DISTANCE_EXON_MAX_LOSS:
                aberrations.append(make('partial_exon_deletion', 'gain_B_donor'))
                aberrations[-1]['_partial_size'] = partial_size
                if partial_size >= native['native_exon_size']:
                    aberrations[-1]['_oversize_partial_deletion'] = True
            elif (partial_size < 0 and partial_size >= -DISTANCE_PARTIAL_RETENTION_MAX
                  and d_from_exon is not None and d_from_exon <= DISTANCE_PARTIAL_RETENTION_MAX):
                aberrations.append(make('partial_intron_retention', 'gain_B_donor'))
                aberrations[-1]['_partial_size'] = partial_size

    return aberrations


def _cryptic_acceptor_orientation_ok(geo_ag, native, is_reverse):
    """
    Reference R parser (spliceAI_parser.R:953):
        + strand: prev_eEnd - GEO_AG < 0      (GEO_AG past prev exon's donor)
        - strand: GEO_AG - prev_eStart < 0     (in R biological order)

    The R parser sorts '-' strand exons in descending genomic order, so its
    `prev_` / `next_` follow transcript (biological) order. Our `native` dict
    stores adjacent exon coords in ASCENDING GENOMIC order. Therefore on '-'
    strand we must swap: use `next_*` (genomically next = biologically previous).

    Returns False when the biological previous exon is absent (variant's
    associated exon is the first biological exon). This matches the R
    reference's NA-propagation, which suppresses partial events at the 5' end
    of the transcript.
    """
    if native is None:
        return False
    if is_reverse:
        # Biologically "prev" = genomically "next" for - strand.
        bio_prev_estart = native.get('next_estart')
        if bio_prev_estart is None:
            return False
        return geo_ag - bio_prev_estart < 0
    else:
        prev_eend = native.get('prev_eend')
        if prev_eend is None:
            return False
        return prev_eend - geo_ag < 0


def _cryptic_donor_orientation_ok(geo_dg, native, is_reverse):
    """
    Reference R parser (spliceAI_parser.R:954):
        + strand: GEO_DG - next_eStart < 0    (GEO_DG before next exon's acceptor)
        - strand: next_eEnd - GEO_DG < 0       (in R biological order)

    Same prev/next swap as _cryptic_acceptor_orientation_ok. Returns False when
    the biological next exon is absent (variant's associated exon is the last
    biological exon), matching the R reference's NA-propagation suppression of
    partial events at the 3' end of the transcript.
    """
    if native is None:
        return False
    if is_reverse:
        # Biologically "next" = genomically "prev" for - strand.
        bio_next_eend = native.get('prev_eend')
        if bio_next_eend is None:
            return False
        return bio_next_eend - geo_dg < 0
    else:
        next_estart = native.get('next_estart')
        if next_estart is None:
            return False
        return geo_dg - next_estart < 0


def sai10k_annotate_frameshift(aberration, cds_start, cds_end, exon_starts, exon_ends):
    """
    Mutate `aberration` in place with size / frameshift / coding fields.

    Adds:
        affects_coding: True | False | None
        frameshift: True | False | None
        size: integer bp of the full affected segment (for reference / display)
        cds_size: integer bp of the CDS-overlapping portion of that segment
        size_type: 'exon' | 'intron' | 'pseudoexon' | 'partial'
        frameshift_description: formatted string for display

    INTENTIONAL DIVERGENCE FROM spliceAI_parser.R: the reference parser computes
    frameshift as (full_affected_size % 3). That is incorrect for exons / introns /
    pseudoexons / partial events whose range straddles a UTR/CDS boundary — only
    the CDS-overlapping bases affect the reading frame of the translated protein.
    This implementation computes frameshift as (cds_size % 3), where cds_size is
    the portion of the affected segment that overlaps [cds_start, cds_end].
    The display label correspondingly shows "<N>bp coding seq." for the CDS-based
    size in coding cases, and "<N>bp" with a "non-coding" suffix otherwise.
    """
    # Treat cds_start==cds_end as non-coding (genePred encodes non-coding
    # transcripts with coincident cds start/end; after 0-based-to-1-based
    # conversion the range collapses or inverts).
    is_coding_transcript = (
        cds_start is not None and cds_end is not None and cds_start < cds_end
    )

    aberration_type = aberration['aberration_type']

    geo_al = aberration.get('geo_al')
    geo_dl = aberration.get('geo_dl')
    geo_ag = aberration.get('geo_ag')
    geo_dg = aberration.get('geo_dg')
    affected_region = aberration.get('affected_region')

    # Helper: number of bases in the closed segment [a, b] that overlap the CDS.
    # Returns 0 for non-coding transcripts and for segments that do not overlap [cds_start, cds_end].
    def cds_overlap_size(a, b):
        if not is_coding_transcript or a is None or b is None:
            return 0
        lo = min(a, b)
        hi = max(a, b)
        overlap_lo = max(lo, cds_start)
        overlap_hi = min(hi, cds_end)
        if overlap_lo > overlap_hi:
            return 0
        return overlap_hi - overlap_lo + 1

    # genePred convention: cds_start and cds_end are GENOMIC bounds, so
    # cds_start < cds_end regardless of strand. Biologically the start codon
    # (ATG) and stop codon each span 3 bases at different genomic endpoints
    # depending on strand:
    #   + strand: ATG = [cds_start, cds_start+2], stop = [cds_end-2, cds_end]
    #   − strand: ATG = [cds_end-2, cds_end],     stop = [cds_start, cds_start+2]
    # Deleting ANY base of a codon destroys it, so the "contains" check tests
    # whether the deleted segment overlaps the 3-base codon range.
    strand = aberration.get('_strand', '+')
    if is_coding_transcript:
        if strand == '-':
            start_codon_lo, start_codon_hi = cds_end - 2, cds_end
            stop_codon_lo, stop_codon_hi = cds_start, cds_start + 2
        else:
            start_codon_lo, start_codon_hi = cds_start, cds_start + 2
            stop_codon_lo, stop_codon_hi = cds_end - 2, cds_end
    else:
        start_codon_lo = start_codon_hi = stop_codon_lo = stop_codon_hi = None

    def _segment_overlaps(a, b, codon_lo, codon_hi):
        if (not is_coding_transcript or a is None or b is None
                or codon_lo is None or codon_hi is None):
            return False
        seg_lo = min(a, b)
        seg_hi = max(a, b)
        return not (seg_hi < codon_lo or seg_lo > codon_hi)

    # Helpers: does the closed segment [a, b] overlap the 3-base start/stop
    # codon range? (Partial overlap of the codon is enough to destroy it.)
    def contains_start_codon(a, b):
        return _segment_overlaps(a, b, start_codon_lo, start_codon_hi)

    def contains_stop_codon(a, b):
        return _segment_overlaps(a, b, stop_codon_lo, stop_codon_hi)

    # Set the five coding-consequence fields (affects_coding, cds_size, frameshift,
    # start_codon_lost, stop_codon_lost) on the aberration dict, and return the
    # status string for the display label.
    #
    # Priority of the status label: start-codon loss > stop-codon loss > frameshift
    # > in-frame > non-coding. The flags match the displayed label (at most one of
    # start_codon_lost / stop_codon_lost is True), so when a single skipped exon
    # contains BOTH cds_start and cds_end, only start_codon_lost is True.
    # Types that structurally can't delete a codon (whole_intron_retention,
    # pseudoexon, increased_exon_inclusion, partial_intron_retention) pass
    # start_lost=False and stop_lost=False.
    def set_coding_fields(aberration, cds_size, start_lost, stop_lost):
        affects_coding = cds_size > 0
        aberration['cds_size'] = cds_size
        aberration['affects_coding'] = affects_coding
        if not affects_coding:
            aberration['frameshift'] = None
            aberration['start_codon_lost'] = None
            aberration['stop_codon_lost'] = None
            return 'non-coding change'
        if start_lost:
            aberration['frameshift'] = None
            aberration['start_codon_lost'] = True
            aberration['stop_codon_lost'] = False
            return 'start codon lost'
        if stop_lost:
            aberration['frameshift'] = None
            aberration['start_codon_lost'] = False
            aberration['stop_codon_lost'] = True
            return 'stop codon lost'
        aberration['frameshift'] = (cds_size % 3 != 0)
        aberration['start_codon_lost'] = False
        aberration['stop_codon_lost'] = False
        return 'frameshift' if aberration['frameshift'] else 'in-frame'

    is_potential = aberration.get('_potential', False)

    # Uniform output-dict contract: every aberration dict carries these fields
    # (value = None when not applicable / unknown). Individual branches below
    # overwrite with type-specific values. Types that structurally cannot lose
    # the start or stop codon (whole_intron_retention, pseudoexon,
    # increased_exon_inclusion, partial_intron_retention) set the codon-lost
    # fields to False when the transcript is coding and the aberration has
    # CDS overlap, and to None on non-coding / unknown branches.
    aberration.setdefault('cds_size', None)
    aberration.setdefault('start_codon_lost', None)
    aberration.setdefault('stop_codon_lost', None)
    aberration.setdefault('size', None)
    aberration.setdefault('size_type', None)

    if aberration_type == 'exon_skipping':
        result = _find_exons_between_loss_positions(geo_al, geo_dl, strand, exon_starts, exon_ends)
        # Deliberate deviation from the reference parser. Both
        # spliceAI_parser.R:227-246 and spliceai_parser.py:293-334 require
        # BOTH geo_al and geo_dl to sit exactly at exon boundaries; when
        # SpliceAI's DP_AL/DP_DL peak lands interior to the affected exon, or
        # one peak lands on a cryptic position while the other matches a
        # native site, the reference's exact-match lookup returns "NA|NA".
        # We instead anchor on whichever Δ-score position lands on a
        # canonical native site of the variant's affected exon and treat
        # that single exon as the skipped one (per developer reply to
        # PDF feedback comment #6: "anchor on whichever Δ-score position
        # lands on a canonical native site and ignore the cryptic position
        # for localization purposes"). Covers both exonic variants
        # (e.g. RAD51C c.404G>A — DP_DL=0 anchors on exon 2's native donor)
        # and intronic splice-site variants (e.g. RAD51C c.404+2T>C —
        # DP_DL=-2 anchors on exon 2's native donor while DP_AL points to
        # a cryptic acceptor position).
        if not result['match'] and affected_region:
            is_reverse = strand == '-'
            region_type = affected_region.get('region_type')
            if region_type == 'exon':
                idx = affected_region.get('_genomic_index')
            elif region_type == 'intron':
                # _genomic_index is the upstream-genomic-exon index for the
                # intron containing the variant. nearest_boundary names the
                # transcript-orientation splice site the variant is closest
                # to. Map (strand, nearest_boundary) -> affected exon index:
                #   + strand, donor    -> intron_idx     (upstream exon i)
                #   + strand, acceptor -> intron_idx + 1 (downstream exon i+1)
                #   - strand, donor    -> intron_idx + 1 (donor of exon i+1 sits at exon_starts[i+1])
                #   - strand, acceptor -> intron_idx     (acceptor of exon i sits at exon_ends[i])
                intron_idx = affected_region.get('_genomic_index')
                nearest = affected_region.get('nearest_boundary')
                if intron_idx is None:
                    idx = None
                elif nearest == 'donor':
                    idx = intron_idx if not is_reverse else intron_idx + 1
                elif nearest == 'acceptor':
                    idx = intron_idx + 1 if not is_reverse else intron_idx
                else:
                    idx = None
            else:
                idx = None
            if idx is not None and 0 <= idx < len(exon_starts):
                acc_pos = exon_ends[idx] if is_reverse else exon_starts[idx]
                don_pos = exon_starts[idx] if is_reverse else exon_ends[idx]
                if geo_al == acc_pos or geo_dl == don_pos:
                    num_exons = len(exon_starts)
                    bio_num = num_exons - idx if is_reverse else idx + 1
                    result = {
                        'genomic_indices': [idx],
                        'biological_numbers': [bio_num],
                        'total_size': exon_ends[idx] - exon_starts[idx] + 1,
                        'match': True,
                    }
        if result['match']:
            # Persist the skipped genomic indices so the post-annotation premature-stop
            # detector in sai10k_compute_predictions can rebuild the altered mRNA. The
            # _* prefix marks this as internal — the strip block at the end of
            # sai10k_compute_predictions removes it before the response is serialized.
            aberration['_skipped_indices'] = list(result['genomic_indices'])
            total_size = result['total_size']
            cds_size = sum(
                cds_overlap_size(exon_starts[i], exon_ends[i])
                for i in result['genomic_indices']
            )
            start_lost = any(
                contains_start_codon(exon_starts[i], exon_ends[i])
                for i in result['genomic_indices']
            )
            stop_lost = any(
                contains_stop_codon(exon_starts[i], exon_ends[i])
                for i in result['genomic_indices']
            )
            status = set_coding_fields(aberration, cds_size, start_lost, stop_lost)
            aberration['size'] = total_size
            aberration['size_type'] = 'exon'
            if is_potential:
                label = 'Potential exon' if len(result['biological_numbers']) == 1 else 'Potential exons'
            else:
                label = 'Exon' if len(result['biological_numbers']) == 1 else 'Exons'
            nums = ', '.join(str(n) for n in result['biological_numbers'])
            size_text = f'{cds_size}bp coding seq.' if aberration['affects_coding'] else f'{total_size}bp'
            aberration['frameshift_description'] = f'{label} {nums} skipping ({size_text}) - {status}'
        else:
            aberration['affects_coding'] = None
            aberration['frameshift'] = None
            aberration['frameshift_description'] = 'Exon skipping predicted but skipped exon(s) could not be mapped'
        return

    if aberration_type == 'whole_intron_retention':
        if geo_al is not None and geo_dl is not None:
            # max(0, ...) guards the degenerate geo_al == geo_dl case where
            # abs(0) - 1 would otherwise produce a negative size in the emitted
            # aberration.
            intron_size = max(0, abs(geo_al - geo_dl) - 1)
            # Segment = retained intron: from min(geo_al, geo_dl)+1 to max(...)-1
            seg_lo = min(geo_al, geo_dl) + 1
            seg_hi = max(geo_al, geo_dl) - 1
            cds_size = cds_overlap_size(seg_lo, seg_hi)
            # Intron retention adds bases; it cannot delete a codon.
            status = set_coding_fields(aberration, cds_size, start_lost=False, stop_lost=False)
            aberration['size'] = intron_size
            aberration['size_type'] = 'intron'
            # affected_intron_number is set at emission (loss branch) when the
            # geo_dl/geo_al pair matches the native sites of a single intron;
            # the emit path in sai10k_determine_aberrations only appends the
            # aberration when that match returns a non-None intron number.
            intron_num = aberration['affected_intron_number']
            prefix = 'Potential intron' if is_potential else 'Intron'
            label = f'{prefix} {intron_num} retention'
            size_text = f'{cds_size}bp coding seq.' if aberration['affects_coding'] else f'{intron_size}bp'
            aberration['frameshift_description'] = f'{label} ({size_text}) - {status}'
        else:
            aberration['affects_coding'] = None
            aberration['frameshift'] = None
            aberration['frameshift_description'] = 'Whole intron retention - size unclear'
        return

    if aberration_type == 'pseudoexon':
        if geo_ag is not None and geo_dg is not None:
            # GEX_size formula: (DP_DG*Strand) - (DP_AG*Strand) + 1. Orientation
            # was verified before firing this branch, so abs() is equivalent.
            gex_size = abs(geo_dg - geo_ag) + 1
            seg_lo = min(geo_ag, geo_dg)
            seg_hi = max(geo_ag, geo_dg)
            cds_size = cds_overlap_size(seg_lo, seg_hi)
            # Pseudoexon activation adds bases; it cannot delete a codon.
            status = set_coding_fields(aberration, cds_size, start_lost=False, stop_lost=False)
            aberration['size'] = gex_size
            aberration['size_type'] = 'pseudoexon'
            size_text = f'{cds_size}bp coding seq.' if aberration['affects_coding'] else f'{gex_size}bp'
            aberration['frameshift_description'] = f'Pseudoexon activation ({size_text}) - {status}'
        else:
            aberration['affects_coding'] = None
            aberration['frameshift'] = None
            aberration['frameshift_description'] = 'Pseudoexon activation - size unclear'
        return

    if aberration_type == 'increased_exon_inclusion':
        # Mirrors the reference parser's Exon_with_increased_inclusion column
        # (spliceAI_parser.R / spliceai_parser.py). Default to None and override
        # below when assoc_idx is in range.
        aberration['included_exon_number'] = None
        native = aberration.get('_native') or {}
        native_exon_size = native.get('native_exon_size')
        if native_exon_size is not None:
            aberration['size'] = native_exon_size
            aberration['size_type'] = 'exon'
            # Including a native exon verbatim cannot change the frame or
            # delete a codon, so we write the five fields directly here
            # rather than going through set_coding_fields (which would set
            # frameshift based on cds_size % 3, which is irrelevant for IEI).
            idx = native.get('assoc_idx')
            if idx is not None and 0 <= idx < len(exon_starts):
                num_exons = len(exon_starts)
                aberration['included_exon_number'] = (
                    num_exons - idx if strand == '-' else idx + 1
                )
                cds_size = cds_overlap_size(exon_starts[idx], exon_ends[idx])
                affects_coding = cds_size > 0
                aberration['affects_coding'] = affects_coding
                aberration['cds_size'] = cds_size
                aberration['frameshift'] = False if affects_coding else None
                aberration['start_codon_lost'] = False if affects_coding else None
                aberration['stop_codon_lost'] = False if affects_coding else None
                label_status = 'coding' if affects_coding else 'non-coding change'
            else:
                aberration['affects_coding'] = None
                aberration['cds_size'] = None
                aberration['frameshift'] = None
                aberration['start_codon_lost'] = None
                aberration['stop_codon_lost'] = None
                label_status = 'unknown'
            aberration['frameshift_description'] = (
                f'Increased exon inclusion ({native_exon_size}bp) - {label_status}')
        else:
            aberration['affects_coding'] = None
            aberration['frameshift'] = None
            aberration['frameshift_description'] = 'Increased exon inclusion - size unclear'
        return

    if aberration_type in ('partial_exon_deletion', 'partial_intron_retention'):
        partial_size = aberration.get('_partial_size')  # signed: +ve deletion, -ve retention
        if partial_size is None:
            aberration['affects_coding'] = None
            aberration['frameshift'] = None
            aberration['frameshift_description'] = f'{aberration_type.replace("_", " ")} - size unclear'
            return
        native = aberration.get('_native') or {}
        # Resolve the human-readable affected exon number (1-based, transcript
        # orientation) so partial_exon_deletion descriptions can read
        # "partial exon N deletion". Falls back to no number if assoc_idx is
        # unavailable.
        assoc_idx = native.get('assoc_idx')
        if assoc_idx is not None and 0 <= assoc_idx < len(exon_starts):
            num_exons = len(exon_starts)
            affected_exon_number = num_exons - assoc_idx if strand == '-' else assoc_idx + 1
        else:
            affected_exon_number = None
        if aberration_type == 'partial_exon_deletion':
            aberration['affected_exon_number'] = affected_exon_number
            consequence = (f'partial exon {affected_exon_number} deletion'
                           if affected_exon_number is not None
                           else 'partial exon deletion')
        else:
            # partial_intron_retention: gain_B_acceptor activates a cryptic
            # acceptor in the upstream intron (in transcript order) of the
            # variant's associated exon; gain_B_donor activates a cryptic
            # donor in the downstream intron. Biological intron N sits
            # between exon N and exon N+1, so:
            #   acceptor branch -> intron (bio_exon - 1)
            #   donor branch    -> intron bio_exon
            # The orientation gates upstream guarantee the relevant neighbor
            # exon exists, so intron numbering won't underflow.
            if affected_exon_number is not None:
                affected_intron_number = (affected_exon_number - 1
                                          if aberration.get('_branch') == 'gain_B_acceptor'
                                          else affected_exon_number)
            else:
                affected_intron_number = None
            aberration['affected_intron_number'] = affected_intron_number
            consequence = (f'partial intron {affected_intron_number} retention'
                           if affected_intron_number is not None
                           else 'partial intron retention')
        # Cryptic acceptor/donor site lands past the opposite native exon
        # boundary. The reference parser still emits the call but its
        # get_partial_seq step returns "deletion greater than exon size"
        # (spliceai_parser.py:552-554), nullifying the cds_size/frameshift
        # downstream. Mirror that behavior here.
        if aberration.get('_oversize_partial_deletion'):
            aberration['size'] = abs(partial_size)
            aberration['size_type'] = 'partial'
            aberration['affects_coding'] = None
            aberration['frameshift'] = None
            aberration['start_codon_lost'] = None
            aberration['stop_codon_lost'] = None
            shift_type = 'Acceptor' if aberration.get('_branch') == 'gain_B_acceptor' else 'Donor'
            aberration['frameshift_description'] = (
                f'{shift_type} shift ({abs(partial_size)}bp {consequence}) '
                f'- deletion greater than exon size')
            return
        display_size = abs(partial_size)
        # Identify the native site and the predicted cryptic site for the active side.
        branch = aberration.get('_branch', '')
        if branch == 'gain_B_acceptor':
            native_site = native.get('geo_na')
            cryptic_site = geo_ag
            shift_type = 'Acceptor'
        else:
            native_site = native.get('geo_nd')
            cryptic_site = geo_dg
            shift_type = 'Donor'
        # Compute the inclusive genomic range of bases actually affected:
        #   - partial_exon_deletion: bases deleted from the exon. Runs from the
        #     native site toward the cryptic site but EXCLUDES the cryptic site,
        #     which is the new splice boundary and remains in the mature mRNA.
        #   - partial_intron_retention: intronic bases now retained in the mRNA.
        #     Runs from the cryptic site toward the native site but EXCLUDES the
        #     native site (which is the original splice boundary).
        # The computed range has length == display_size by construction.
        if native_site is not None and cryptic_site is not None:
            if aberration_type == 'partial_exon_deletion':
                if native_site < cryptic_site:
                    seg_lo, seg_hi = native_site, cryptic_site - 1
                else:
                    seg_lo, seg_hi = cryptic_site + 1, native_site
            else:  # partial_intron_retention
                if cryptic_site < native_site:
                    seg_lo, seg_hi = cryptic_site, native_site - 1
                else:
                    seg_lo, seg_hi = native_site + 1, cryptic_site
            cds_size = cds_overlap_size(seg_lo, seg_hi)
            # Start/stop-codon-lost only applies to partial_exon_deletion
            # (partial_intron_retention adds bases to the mRNA rather than
            # removing any; it can't delete the start or stop codon).
            if aberration_type == 'partial_exon_deletion':
                start_lost = contains_start_codon(seg_lo, seg_hi)
                stop_lost = contains_stop_codon(seg_lo, seg_hi)
            else:
                start_lost = False
                stop_lost = False
        else:
            cds_size = 0
            start_lost = False
            stop_lost = False
        status = set_coding_fields(aberration, cds_size, start_lost, stop_lost)
        aberration['size'] = display_size
        aberration['size_type'] = 'partial'
        size_text = f'{cds_size}bp coding seq.' if aberration['affects_coding'] else f'{display_size}bp'
        aberration['frameshift_description'] = f'{shift_type} shift ({size_text} {consequence}) - {status}'
        return

    # Unknown aberration type — leave blank (should not happen with current paper-only emission).
    aberration['affects_coding'] = None
    aberration['frameshift'] = None
    aberration['frameshift_description'] = f'{aberration_type} - frameshift status unclear'


def sai10k_compute_predictions(transcript_scores, variant_pos, chrom=None, ref=None, alt=None, fasta=None):
    """
    Compute splice aberration predictions for a single transcript.

    Args:
        transcript_scores: dict with SpliceAI delta scores, delta positions,
            raw alt-allele site scores, and transcript structure.
        variant_pos: 1-based genomic position of the variant.
        chrom: variant chromosome (e.g. "17"). Optional; required for premature-
            stop detection.
        ref: variant reference allele. Optional; required for premature-stop
            detection.
        alt: variant alternate allele. Optional; required for premature-stop
            detection.
        fasta: pyfastx Fasta handle for the reference genome. Optional; when
            None, premature-stop detection is skipped and stop_codon_introduced
            / aa_change are set to None on every aberration.

    Returns:
        dict {
            'aberrations': [list of aberration dicts, possibly empty],
            'transcript_info': {...}
        }
    """
    ds_ag = transcript_scores.get('DS_AG', 0)
    ds_al = transcript_scores.get('DS_AL', 0)
    ds_dg = transcript_scores.get('DS_DG', 0)
    ds_dl = transcript_scores.get('DS_DL', 0)

    dp_ag = transcript_scores.get('DP_AG', 0)
    dp_al = transcript_scores.get('DP_AL', 0)
    dp_dg = transcript_scores.get('DP_DG', 0)
    dp_dl = transcript_scores.get('DP_DL', 0)

    # Standard SpliceAI output omits raw alt-allele site scores (DS_*_ALT);
    # they're only present when the caller supplies the pre-computed ALT
    # annotations. When absent, `sai10k_determine_aberrations` falls back to
    # the reference parser's --include=False dominance logic instead of
    # treating the missing values as zeros (which would silently route through
    # the ALT-score side-selection branch).
    has_alt_scores = all(
        k in transcript_scores
        for k in ('DS_AG_ALT', 'DS_AL_ALT', 'DS_DG_ALT', 'DS_DL_ALT')
    )
    ds_ag_alt = transcript_scores.get('DS_AG_ALT', 0)
    ds_al_alt = transcript_scores.get('DS_AL_ALT', 0)
    ds_dg_alt = transcript_scores.get('DS_DG_ALT', 0)
    ds_dl_alt = transcript_scores.get('DS_DL_ALT', 0)

    exon_starts = transcript_scores.get('EXON_STARTS', [])
    exon_ends = transcript_scores.get('EXON_ENDS', [])
    cds_start = transcript_scores.get('CDS_START')
    cds_end = transcript_scores.get('CDS_END')
    strand = transcript_scores.get('STRAND') or transcript_scores.get('t_strand')
    if strand not in ('+', '-'):
        # Falling back to '+' on a '-' strand transcript would silently flip
        # the start/stop codon lost flags (biological ATG lives at cds_end on
        # '-' strand). Warn loudly if we have to guess.
        print(
            f"WARNING: transcript has no STRAND/t_strand field; defaulting to '+'. "
            f"Codon-lost flags may be wrong if the transcript is actually on '-'. "
            f"Transcript id = {transcript_scores.get('t_id') or transcript_scores.get('NAME')}"
        )
        strand = '+'

    t0 = time.perf_counter()
    aberrations = sai10k_determine_aberrations(
        ds_ag, ds_al, ds_dg, ds_dl,
        dp_ag, dp_al, dp_dg, dp_dl,
        ds_ag_alt, ds_al_alt, ds_dg_alt, ds_dl_alt,
        variant_pos, exon_starts, exon_ends, strand,
        has_alt_scores=has_alt_scores,
    )
    determine_ms = (time.perf_counter() - t0) * 1000

    # Pick the overall delta_type (the Δ branch with the highest Δ score). This
    # is constant per variant and labels which SpliceAI score drove the call.
    _, overall_delta_type = max(
        ((_coerce_float(ds_al), 'acceptor loss'),
         (_coerce_float(ds_dl), 'donor loss'),
         (_coerce_float(ds_ag), 'acceptor gain'),
         (_coerce_float(ds_dg), 'donor gain')),
        key=lambda item: item[0],
    )

    t1 = time.perf_counter()
    for aberration in aberrations:
        sai10k_annotate_frameshift(
            aberration,
            cds_start, cds_end,
            exon_starts, exon_ends,
        )
        aberration['delta_type'] = overall_delta_type
    annotate_ms = (time.perf_counter() - t1) * 1000

    # Premature-stop detection (PDF feedback comment #5). Runs after annotation
    # so it can read affects_coding / frameshift / start_codon_lost /
    # stop_codon_lost from the annotation step. Mutates each aberration:
    #   stop_codon_introduced:    True | False | None
    #   aa_change:                formatted "ABC[XYZ*]" string, the
    #                             EXTENDS_PAST_NATIVE_STOP sentinel for a
    #                             frameshift extension, or None.
    #   wt_protein_window:        ±15-aa context around the change for the WT
    #                             protein, with truncation/total markers (PTC
    #                             cases only; None otherwise).
    #   altered_protein_window:   same window for the predicted altered protein.
    # When the detector signals a PTC, we extend frameshift_description inline.
    # The detector silently returns (None, None, None, None) when fasta=None
    # (no genome FASTA available) or when chrom/ref/alt aren't supplied —
    # preserves backward compatibility for any caller that doesn't supply them.
    t2 = time.perf_counter()
    n_premature_stop_calls = 0
    # Compute the per-transcript consensus context lazily on the first
    # qualifying aberration, then reuse for the rest. Avoids re-fetching the
    # consensus CDS sequence and re-translating the full consensus K times for
    # transcripts with K aberrations (e.g. TTN with K=4).
    consensus_ctx = None
    consensus_ctx_computed = False
    for aberration in aberrations:
        if (aberration.get('affects_coding') is True
                and not aberration.get('start_codon_lost')
                and not aberration.get('stop_codon_lost')):
            if not consensus_ctx_computed:
                consensus_ctx = _compute_consensus_context(transcript_scores, fasta, chrom)
                consensus_ctx_computed = True
            n_premature_stop_calls += 1
            stop_introduced, aa_change, wt_window, altered_window = _detect_premature_stop(
                transcript_scores, aberration, fasta, chrom, ref, alt, variant_pos,
                consensus_ctx=consensus_ctx,
            )
            aberration['stop_codon_introduced'] = stop_introduced
            aberration['aa_change'] = aa_change
            aberration['wt_protein_window'] = wt_window
            aberration['altered_protein_window'] = altered_window
            if aberration.get('frameshift_description'):
                if stop_introduced is True and aberration.get('frameshift') is not True:
                    aberration['frameshift_description'] += ' but introduces a stop codon'
                elif stop_introduced is True and aberration.get('frameshift') is True:
                    aberration['frameshift_description'] += ', introduces a stop codon'
                elif (aberration.get('frameshift') is True
                        and aa_change == EXTENDS_PAST_NATIVE_STOP):
                    aberration['frameshift_description'] += ' and extends past the native stop'
        else:
            aberration['stop_codon_introduced'] = None
            aberration['aa_change'] = None
            aberration['wt_protein_window'] = None
            aberration['altered_protein_window'] = None
    premature_stop_ms = (time.perf_counter() - t2) * 1000

    # Strip internal-only fields before returning to callers (and ultimately the
    # client). Keys prefixed with '_' are implementation details.
    for aberration in aberrations:
        for key in list(aberration.keys()):
            if key.startswith('_'):
                del aberration[key]
        region = aberration.get('affected_region')
        if region:
            for key in list(region.keys()):
                if key.startswith('_'):
                    del region[key]

    return {
        'aberrations': aberrations,
        'transcript_info': {
            'strand': strand,
            'tx_start': transcript_scores.get('TX_START'),
            'tx_end': transcript_scores.get('TX_END'),
            'cds_start': cds_start,
            'cds_end': cds_end,
            'num_exons': len(exon_starts) if exon_starts else 0,
            'is_coding': (
                cds_start is not None and cds_end is not None and cds_start < cds_end
            ),
        },
        # Internal timing breakdown for the server's SAI10K_TIMING log line.
        # Stripped by sai10k_get_transcript_predictions' caller before the
        # response is serialized to the client.
        '_timing_ms': {
            'determine': determine_ms,
            'annotate': annotate_ms,
            'premature_stop': premature_stop_ms,
            'n_aberrations': len(aberrations),
            'n_premature_stop_calls': n_premature_stop_calls,
        },
    }


def sai10k_select_transcript(scores_list):
    """
    Select the highest-priority transcript from a list.

    Selection rule (single source of truth for SpliceAI-lookup):
        1. Highest transcript priority (MS > MP > C > N).
        2. Tie-break: highest SUM of |DS_AG| + |DS_AL| + |DS_DG| + |DS_DL|.

    Returns the selected transcript dict, or None if scores_list is empty.
    """
    if not scores_list:
        return None

    best = None
    best_priority = -1
    best_sum = -1.0

    for transcript_scores in scores_list:
        priority = TRANSCRIPT_PRIORITY_ORDER.get(transcript_scores.get('t_priority', 'N'), 0)
        score_sum = sum(
            abs(_coerce_float(transcript_scores.get(key, 0)))
            for key in ('DS_AG', 'DS_AL', 'DS_DG', 'DS_DL')
        )
        if priority > best_priority or (priority == best_priority and score_sum > best_sum):
            best = transcript_scores
            best_priority = priority
            best_sum = score_sum

    return best


def sai10k_get_transcript_predictions(transcript_scores, variant_pos, chrom=None, ref=None, alt=None, fasta=None):
    """Compute predictions for one transcript and attach its metadata.

    Optional chrom/ref/alt/fasta enable the premature-stop detector. When fasta
    is None the detector is skipped (stop_codon_introduced / aa_change set to
    None on every aberration), preserving behavior for callers that don't
    supply a genome FASTA.
    """
    if not transcript_scores:
        return None

    predictions = sai10k_compute_predictions(
        transcript_scores, variant_pos, chrom=chrom, ref=ref, alt=alt, fasta=fasta,
    )
    predictions['transcript_id'] = transcript_scores.get('t_id') or transcript_scores.get('NAME')
    predictions['transcript_priority'] = transcript_scores.get('t_priority', 'N')
    predictions['gene_name'] = transcript_scores.get('g_name')
    return predictions


# ===========================================================================
# Premature-stop detection helpers (PDF feedback comment #5).
#
# These mirror the reference parser at
#   ~/code/SpliceAI-lookup-dev/SAI-10k-calc/spliceai_parser.py
# (functions add_variant, find_difference_point, determine_aa_seq, and the
# per-type get_*_seq functions). The reference operates on pandas DataFrames
# with refseq-specific columns; here we work directly off the existing
# EXON_STARTS / EXON_ENDS / CDS_START / CDS_END / STRAND fields plus the
# variant chrom/pos/ref/alt and a pyfastx FASTA handle.
#
# See also ~/code/SpliceAI-lookup-dev/.memory/PLAN_stop_codon_introduced.md.
# ===========================================================================

_CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

_COMPLEMENT_TRANS = str.maketrans('ACGTacgt', 'TGCAtgca')


def _translate_dna(dna_seq):
    """Translate DNA to single-letter AA. Stop codons -> '*', unknowns -> 'X'."""
    dna_seq = dna_seq.upper()
    protein = []
    for i in range(0, len(dna_seq) - 2, 3):
        protein.append(_CODON_TABLE.get(dna_seq[i:i+3], 'X'))
    return ''.join(protein)


def _reverse_complement(seq):
    """Reverse-complement a DNA sequence (preserves case, passes through any non-ACGT base)."""
    return seq.translate(_COMPLEMENT_TRANS)[::-1]


def _find_difference_point(seq_a, seq_b, direction):
    """Return 1-based position of the first differing character, or 'identical'.

    direction='forward' walks from the start; direction='reverse' walks from the
    end and returns the position back-translated into forward coordinates of
    seq_a. When one sequence is a prefix of the other, returns the length of
    the shorter sequence + 1. Mirrors spliceai_parser.py:136-173.
    """
    if seq_a == seq_b:
        return 'identical'

    a, b = seq_a, seq_b
    if direction == 'reverse':
        a = a[::-1]
        b = b[::-1]

    min_len = min(len(a), len(b))
    point = None
    for i in range(min_len):
        if a[i] != b[i]:
            point = i + 1
            break
    if point is None:
        point = min_len + 1

    if direction == 'reverse':
        point = len(seq_a) - point + 1

    return point


# Cache: id(fasta_handle) -> bool. True = FASTA contigs are named with a "chr"
# prefix; False = bare numeric / X / Y. Probed once per handle on first fetch.
# Keyed by id() so we don't accidentally hold a strong reference that prevents
# the FASTA handle from being garbage-collected; the cache entries leak with
# the handle but that's fine — these are long-lived singletons in server.py.
_FASTA_CHROM_PREFIX_CACHE = {}


def _fetch_seq(fasta, chrom, start, end):
    """Fetch [start, end] (1-based, inclusive) from the FASTA, normalizing chrom prefix.

    server.py's VARIANT_RE strips the leading 'chr' from incoming variants, but
    the hg19/hg38 FASTAs in the container index contigs as 'chr1..chrX'. This
    helper probes the FASTA's keys once on first call and caches whether
    'chr' should be prepended; all FASTA reads in the new helpers must go
    through this function.
    """
    fasta_id = id(fasta)
    if fasta_id not in _FASTA_CHROM_PREFIX_CACHE:
        try:
            sample_keys = list(fasta.keys())[:5]
        except Exception:
            sample_keys = []
        _FASTA_CHROM_PREFIX_CACHE[fasta_id] = any(str(k).startswith('chr') for k in sample_keys)

    needs_chr = _FASTA_CHROM_PREFIX_CACHE[fasta_id]
    chrom_str = str(chrom)
    if needs_chr and not chrom_str.startswith('chr'):
        chrom_str = 'chr' + chrom_str
    elif not needs_chr and chrom_str.startswith('chr'):
        chrom_str = chrom_str[3:]
    return str(fasta[chrom_str][start - 1:end])


def _apply_variant_to_altered_exon_seqs(altered_exon_spans, exon_seqs, var_pos, ref, alt,
                                        fasta, chrom, cds_start, cds_end):
    """Apply variant to the affected exon sequence (in genomic orientation).

    Returns one of:
        "<exon_idx>_<modified_seq>"      — success; caller writes back into exon_seqs.
        "reference mismatch"             — `ref` doesn't match the genome.
        "impacts native start or stop site" — variant footprint hits the 6 bp of
                                             the start codon or stop codon.
        "variant straddles splice boundary" — footprint partially overlaps a
                                             kept span (indels only). Caller skips.
        "cannot determine"               — footprint is entirely outside every
                                             kept span (variant in a non-retained
                                             intron). Caller proceeds with the
                                             unmodified altered sequences.

    Generalizes spliceai_parser.py:89-133 to support indels: spliceing
        `affected_exon[:adjustment_size] + alt + affected_exon[adjustment_size + len(ref):]`
    instead of hardcoding `len(ref) == 1`.
    """
    var_end = var_pos + len(ref) - 1

    # Reference check.
    genome_ref = _fetch_seq(fasta, chrom, var_pos, var_end).upper()
    if genome_ref != ref.upper():
        return 'reference mismatch'

    # Native start / stop codon footprint check (6 bp total). The two 3-bp
    # CDS-terminal windows are checked symmetrically: on + strand the low
    # window is the start codon and the high window is the stop, on - strand
    # the roles flip. Since both branches return the same sentinel, the
    # naming stays strand-neutral so the symmetry is obvious.
    if cds_start is not None and cds_end is not None:
        cds_lo_codon_lo, cds_lo_codon_hi = cds_start, cds_start + 2
        cds_hi_codon_lo, cds_hi_codon_hi = cds_end - 2, cds_end
        if not (var_end < cds_lo_codon_lo or var_pos > cds_lo_codon_hi):
            return 'impacts native start or stop site'
        if not (var_end < cds_hi_codon_lo or var_pos > cds_hi_codon_hi):
            return 'impacts native start or stop site'

    # Find a kept span that fully contains the variant footprint.
    fully_contained_idx = None
    partial_overlap = False
    for i, (span_start, span_end) in enumerate(altered_exon_spans):
        if var_pos >= span_start and var_end <= span_end:
            fully_contained_idx = i
            break
        if not (var_end < span_start or var_pos > span_end):
            partial_overlap = True

    if fully_contained_idx is None:
        if partial_overlap:
            return 'variant straddles splice boundary'
        return 'cannot determine'

    span_start, span_end = altered_exon_spans[fully_contained_idx]

    # VCF-anchored insertion at the trailing edge of a kept span (e.g.
    # `pos=20, ref='C', alt='CTAA'` for a span ending at 20): the anchor
    # is the last exonic base, so by VCF convention the inserted bases sit
    # AFTER position 20 — i.e., past the splice donor in the intron. We
    # cannot faithfully apply that to the altered mRNA without modeling
    # the splicing details, so treat it the same as an indel that
    # straddles a splice boundary and skip detection.
    if len(alt) > len(ref) and var_end == span_end:
        return 'variant straddles splice boundary'

    affected_exon = exon_seqs[fully_contained_idx]
    adjustment_size = var_pos - span_start
    adjusted = affected_exon[:adjustment_size] + alt + affected_exon[adjustment_size + len(ref):]
    return f'{fully_contained_idx}_{adjusted}'


def _consensus_filtered_in_tx_order(transcript_scores):
    """Return the consensus exon list in transcription order, CDS-clamped, sans UTR-only exons.

    Each entry is a dict {orig_start, orig_end, cds_lo, cds_hi, eNum} with the
    biological (1-based, transcription-order) exon number. UTR-only exons
    (those with zero CDS overlap) are dropped — matches the reference's
    `eFrame != -1` filter at spliceai_parser.py:208.
    """
    exon_starts = transcript_scores.get('EXON_STARTS', [])
    exon_ends = transcript_scores.get('EXON_ENDS', [])
    cds_start = transcript_scores.get('CDS_START')
    cds_end = transcript_scores.get('CDS_END')
    strand = transcript_scores.get('STRAND') or transcript_scores.get('t_strand')

    if not exon_starts or cds_start is None or cds_end is None:
        return []

    n = len(exon_starts)
    out = []
    for gi in range(n):
        cds_lo = max(exon_starts[gi], cds_start)
        cds_hi = min(exon_ends[gi], cds_end)
        if cds_lo > cds_hi:
            continue
        e_num = (n - gi) if strand == '-' else (gi + 1)
        out.append({
            'orig_start': exon_starts[gi],
            'orig_end': exon_ends[gi],
            'cds_lo': cds_lo,
            'cds_hi': cds_hi,
            'eNum': e_num,
        })

    if strand == '-':
        out.sort(key=lambda x: -x['orig_start'])
    else:
        out.sort(key=lambda x: x['orig_start'])
    return out


def _compute_consensus_context(transcript_scores, fasta, chrom):
    """Precompute per-transcript state shared across all aberrations.

    Returns a dict with:
        consensus:         _consensus_filtered_in_tx_order list (CDS-clamped, tx order).
        cum_cds_bases:     prefix-sum array; cum_cds_bases[i] = total CDS bases in
                           consensus[:i]. Length len(consensus)+1.
        consensus_aa_full: full consensus AA translated from row_anchor=0 (eFrame=0).
                           None if `fasta` or `chrom` is missing.
        seq_cache:         {(start, end): genomic_seq} populated by the consensus
                           translation. Threaded into the per-aberration altered
                           translation so unchanged exon spans aren't re-fetched.

    For any per-aberration row_anchor `i`, the consensus AA tail equals
    consensus_aa_full[(cum_cds_bases[i] + 2) // 3:] (integer ceil(cum/3) — accounts
    for the (3 - first_eframe) % 3 leading-base trim that the eFrame block of
    `_translate_spans` applies). Hoisting this out of the per-aberration loop
    avoids re-fetching the consensus CDS sequence and re-translating it K times
    (where K is the number of qualifying aberrations on the transcript).

    Returns None if the transcript has no consensus exons (e.g. non-coding) or
    if a FASTA fetch raises — detection is then silently skipped.
    """
    try:
        consensus = _consensus_filtered_in_tx_order(transcript_scores)
        if not consensus:
            return None

        cum = [0]
        for ex in consensus:
            cum.append(cum[-1] + ex['cds_hi'] - ex['cds_lo'] + 1)

        consensus_aa_full = None
        seq_cache = {}
        if fasta is not None and chrom is not None:
            strand = transcript_scores.get('STRAND') or transcript_scores.get('t_strand')
            consensus_spans = [(ex['cds_lo'], ex['cds_hi']) for ex in consensus]
            consensus_aa_full = _translate_spans(
                consensus_spans, 0, fasta, chrom, strand, seq_cache=seq_cache,
            )

        return {
            'consensus': consensus,
            'cum_cds_bases': cum,
            'consensus_aa_full': consensus_aa_full,
            'seq_cache': seq_cache,
        }
    except Exception as exc:
        print(f'WARNING: SAI-10k consensus precompute raised '
              f'{type(exc).__name__}: {exc}')
        return None


def _build_altered_exon_spans(transcript_scores, aberration, consensus=None):
    """Build the per-type altered exon list for premature-stop detection.

    Returns (altered_spans, row_anchor, first_eframe), where altered_spans is
    a list of (1-based, inclusive, CDS-clamped) genomic exon spans for the slice
    starting at row_anchor in transcription order, OR None to signal "skip
    detection" (caller leaves stop_codon_introduced = None).

    Mirrors per-type table builders in spliceai_parser.py — see the per-type
    table in PLAN_stop_codon_introduced.md.

    `consensus`, when provided, is a precomputed _consensus_filtered_in_tx_order
    result; passing it in avoids redundant rebuilds when this is called for
    multiple aberrations on the same transcript.
    """
    if consensus is None:
        consensus = _consensus_filtered_in_tx_order(transcript_scores)
    if not consensus:
        return None

    aberration_type = aberration.get('aberration_type')
    cds_start = transcript_scores.get('CDS_START')
    cds_end = transcript_scores.get('CDS_END')
    strand = transcript_scores.get('STRAND') or transcript_scores.get('t_strand')

    altered_full = None
    row_anchor = None

    if aberration_type == 'exon_skipping':
        skipped_indices = aberration.get('_skipped_indices')
        if not skipped_indices:
            return None
        n_total = len(transcript_scores.get('EXON_STARTS', []))
        skipped_eNums = set()
        for gi in skipped_indices:
            if strand == '-':
                skipped_eNums.add(n_total - gi)
            else:
                skipped_eNums.add(gi + 1)
        altered_full = [ex for ex in consensus if ex['eNum'] not in skipped_eNums]
        if not altered_full:
            return None
        anchor_eNum = min(skipped_eNums) - 1
        if anchor_eNum < 1:
            row_anchor = 0
        else:
            row_anchor = next(
                (i for i, ex in enumerate(altered_full) if ex['eNum'] == anchor_eNum),
                None,
            )
            if row_anchor is None:
                row_anchor = 0

    elif aberration_type == 'whole_intron_retention':
        geo_al = aberration.get('geo_al')
        geo_dl = aberration.get('geo_dl')
        if geo_al is None or geo_dl is None:
            return None
        intron_lo = min(geo_al, geo_dl) + 1
        intron_hi = max(geo_al, geo_dl) - 1
        if intron_lo > intron_hi:
            return None
        # CDS-clamp the retained intron; non-coding intron retention -> None.
        intron_cds_lo = max(intron_lo, cds_start)
        intron_cds_hi = min(intron_hi, cds_end)
        if intron_cds_lo > intron_cds_hi:
            return None
        # Find the consensus exon pair that flanks the retained intron in tx order.
        flanking_idx = None
        for i in range(len(consensus) - 1):
            cur, nxt = consensus[i], consensus[i + 1]
            if strand == '-':
                if cur['orig_start'] == intron_hi + 1 and nxt['orig_end'] == intron_lo - 1:
                    flanking_idx = i
                    break
            else:
                if cur['orig_end'] == intron_lo - 1 and nxt['orig_start'] == intron_hi + 1:
                    flanking_idx = i
                    break
        if flanking_idx is None:
            return None
        intron_row = {
            'orig_start': intron_lo,
            'orig_end': intron_hi,
            'cds_lo': intron_cds_lo,
            'cds_hi': intron_cds_hi,
            'eNum': None,
        }
        altered_full = consensus[:flanking_idx + 1] + [intron_row] + consensus[flanking_idx + 1:]
        # row_number = flanking_idx + 1 (position of intron_row); row_anchor = max(0, row_number - 1).
        row_anchor = flanking_idx

    elif aberration_type in ('partial_intron_retention', 'partial_exon_deletion'):
        native = aberration.get('_native') or {}
        branch = aberration.get('_branch', '')
        if branch == 'gain_B_acceptor':
            native_site = native.get('geo_na')
            cryptic_site = aberration.get('geo_ag')
        elif branch == 'gain_B_donor':
            native_site = native.get('geo_nd')
            cryptic_site = aberration.get('geo_dg')
        else:
            return None
        if native_site is None or cryptic_site is None:
            return None
        # Find the affected exon by matching native_site to a tx-orientation boundary.
        affected_idx = None
        for i, ex in enumerate(consensus):
            if branch == 'gain_B_acceptor':
                target = ex['orig_end'] if strand == '-' else ex['orig_start']
            else:
                target = ex['orig_start'] if strand == '-' else ex['orig_end']
            if target == native_site:
                affected_idx = i
                break
        if affected_idx is None:
            return None
        ex = consensus[affected_idx]
        new_orig_start = ex['orig_start']
        new_orig_end = ex['orig_end']
        if branch == 'gain_B_acceptor':
            if strand == '-':
                new_orig_end = cryptic_site
            else:
                new_orig_start = cryptic_site
        else:
            if strand == '-':
                new_orig_start = cryptic_site
            else:
                new_orig_end = cryptic_site
        new_cds_lo = max(new_orig_start, cds_start)
        new_cds_hi = min(new_orig_end, cds_end)
        if new_cds_lo > new_cds_hi:
            return None
        altered_full = list(consensus)
        altered_full[affected_idx] = {
            'orig_start': new_orig_start,
            'orig_end': new_orig_end,
            'cds_lo': new_cds_lo,
            'cds_hi': new_cds_hi,
            'eNum': ex['eNum'],
        }
        row_anchor = max(0, affected_idx - 1)

    elif aberration_type == 'pseudoexon':
        geo_ag = aberration.get('geo_ag')
        geo_dg = aberration.get('geo_dg')
        if geo_ag is None or geo_dg is None:
            return None
        pseudo_lo = min(geo_ag, geo_dg)
        pseudo_hi = max(geo_ag, geo_dg)
        # Stricter than the upstream affects_coding check: pseudoexon must lie
        # entirely within CDS, otherwise translation of the inserted bases is
        # ill-defined.
        if pseudo_lo < cds_start or pseudo_hi > cds_end:
            return None
        flanking_idx = None
        for i in range(len(consensus) - 1):
            cur, nxt = consensus[i], consensus[i + 1]
            if strand == '-':
                if cur['orig_start'] > pseudo_hi and nxt['orig_end'] < pseudo_lo:
                    flanking_idx = i
                    break
            else:
                if cur['orig_end'] < pseudo_lo and nxt['orig_start'] > pseudo_hi:
                    flanking_idx = i
                    break
        if flanking_idx is None:
            return None
        pseudo_row = {
            'orig_start': pseudo_lo,
            'orig_end': pseudo_hi,
            'cds_lo': pseudo_lo,
            'cds_hi': pseudo_hi,
            'eNum': None,
        }
        altered_full = consensus[:flanking_idx + 1] + [pseudo_row] + consensus[flanking_idx + 1:]
        row_anchor = flanking_idx

    elif aberration_type == 'increased_exon_inclusion':
        native = aberration.get('_native') or {}
        assoc_idx = native.get('assoc_idx')
        if assoc_idx is None:
            return None
        n_total = len(transcript_scores.get('EXON_STARTS', []))
        target_eNum = (n_total - assoc_idx) if strand == '-' else (assoc_idx + 1)
        included_row = next(
            (i for i, ex in enumerate(consensus) if ex['eNum'] == target_eNum),
            None,
        )
        if included_row is None:
            return None
        altered_full = list(consensus)
        row_anchor = max(0, included_row - 1)

    else:
        return None

    if not altered_full or row_anchor is None or row_anchor >= len(altered_full):
        return None

    altered_spans = [(ex['cds_lo'], ex['cds_hi']) for ex in altered_full[row_anchor:]]

    # first_eframe: cumulative CDS bases of consensus entries strictly preceding
    # row_anchor, mod 3. Walking the consensus list (vs. altered_full) is
    # equivalent — for every aberration type, the entries strictly preceding
    # row_anchor are unchanged from consensus by construction.
    first_eframe = sum(ex['cds_hi'] - ex['cds_lo'] + 1 for ex in consensus[:row_anchor]) % 3

    return (altered_spans, row_anchor, first_eframe)


def _translate_spans(spans, first_eframe, fasta, chrom, strand, apply_variant_args=None,
                      seq_cache=None):
    """Translate a list of (1-based, inclusive) genomic spans to a single-letter AA string.

    Order of operations (must match spliceai_parser.py:374-434 — applying the
    variant after reverse-complement would mutate the wrong base on - strand):
        1. Fetch each span in genomic orientation.
        2. If apply_variant_args is given, apply the variant to the genomic-
           orientation sequences. On a non-success sentinel the sentinel is
           returned unchanged.
        3. Reverse-complement per-exon for - strand.
        4. Trim leading bases off the first exon to align reading frame.
        5. Join and translate.

    apply_variant_args, when supplied, is the tuple
    (var_pos, ref, alt, cds_start, cds_end). Returns the AA string on success,
    or one of the sentinels {"reference mismatch", "impacts native start or
    stop site", "variant straddles splice boundary"} when variant application
    fails in a way the caller should skip.

    seq_cache, when supplied, is a write-through dict {(start, end): genomic_seq}
    that lets repeated calls on the same transcript skip pyfastx slices for
    spans already fetched. Only genomic-orientation sequences are stored — the
    reverse-complement, eFrame trim, and variant application all happen on
    per-call copies so cache entries are not mutated.
    """
    if not spans:
        return ''

    if seq_cache is None:
        exon_seqs = [_fetch_seq(fasta, chrom, s, e) for s, e in spans]
    else:
        exon_seqs = []
        for span in spans:
            seq = seq_cache.get(span)
            if seq is None:
                seq = _fetch_seq(fasta, chrom, span[0], span[1])
                seq_cache[span] = seq
            exon_seqs.append(seq)

    if apply_variant_args is not None:
        var_pos, ref, alt, cds_start, cds_end = apply_variant_args
        result = _apply_variant_to_altered_exon_seqs(
            spans, exon_seqs, var_pos, ref, alt, fasta, chrom, cds_start, cds_end,
        )
        if result in ('reference mismatch', 'impacts native start or stop site',
                      'variant straddles splice boundary'):
            return result
        if result != 'cannot determine':
            parts = result.split('_', 1)
            if len(parts) == 2:
                exon_seqs[int(parts[0])] = parts[1]

    if strand == '-':
        exon_seqs = [_reverse_complement(s) for s in exon_seqs]

    # eFrame 1 -> trim 2 leading bases; eFrame 2 -> trim 1; eFrame 0 -> no trim.
    # (Reference convention from spliceai_parser.py:417-426.)
    if exon_seqs and first_eframe == 1 and len(exon_seqs[0]) >= 2:
        exon_seqs[0] = exon_seqs[0][2:]
    elif exon_seqs and first_eframe == 2 and len(exon_seqs[0]) >= 1:
        exon_seqs[0] = exon_seqs[0][1:]

    return _translate_dna(''.join(exon_seqs))


EXTENDS_PAST_NATIVE_STOP = 'protein sequence extends beyond native stop site'


def _build_protein_windows(altered_aa, consensus_aa, frameshift,
                           prefix_offset_aa=0, total_wt_aa=None, flank=15):
    """Build human-readable WT vs. altered protein windows for the popup tooltip.

    Caller contract: only invoked when a PTC has been detected (i.e.
    `_format_aa_change_and_detect` returned `stop_introduced=True`), so
    `altered_aa` is guaranteed to contain '*'. Non-PTC cases are not handled.

    Returns (wt_display, altered_display) — both with the format
    "[NAA]...VISIBLE...[NAA] (TotalAA total)", where the leading/trailing
    "[NAA]..." markers indicate residues truncated from the displayed window
    and "(TotalAA total)" gives the full-protein residue count (excluding the
    natural stop). The altered line terminates at the bracketed PTC and has no
    trailing residues. Returns (None, None) if the two translations are
    identical (defensive — shouldn't occur given the caller contract).

    Inputs `altered_aa` and `consensus_aa` are the per-aberration translations
    starting at row_anchor; `prefix_offset_aa` is the count of WT residues that
    precede consensus_aa in the full protein; `total_wt_aa` is the full-protein
    WT residue count (excluding terminator). `flank` controls how many residues
    of context are shown on either side of the change.

    Bracket convention in altered_display: '[' opens at the first divergent
    residue (forward_diff) and ']' closes after the premature stop. This is
    intentionally aligned to the biological change boundary, not to the fixed
    3-residue prefix used by `_format_aa_change_and_detect`'s compact
    `aa_change`.
    """
    # Mirror _format_aa_change_and_detect's trailing-* trim. Consensus always
    # strips the natural stop; altered only strips it when not frameshift
    # (frameshift extension reads the natural stop's bases out of frame, so a
    # trailing '*' there is meaningful — irrelevant for PTC cases since the
    # first '*' is always the PTC, but kept for parity with the format helper).
    if consensus_aa.endswith('*'):
        consensus_aa = consensus_aa[:-1]
    if not frameshift and altered_aa.endswith('*'):
        altered_aa = altered_aa[:-1]

    forward_diff_raw = _find_difference_point(altered_aa, consensus_aa, 'forward')
    if forward_diff_raw == 'identical' or not isinstance(forward_diff_raw, int):
        return (None, None)
    # _find_difference_point returns 1-based forward positions; convert to a
    # 0-based index of the first divergent residue for slicing.
    diff_idx = forward_diff_raw - 1

    # Caller contract guarantees altered_aa contains '*'; find('*') returns the
    # 0-based PTC position. The changed region ends one past the PTC; the
    # consensus has no symmetric end so we just render flank residues
    # downstream of diff_idx in the WT.
    altered_stop_local = altered_aa.find('*')
    altered_change_end = altered_stop_local + 1
    consensus_change_end = diff_idx

    if total_wt_aa is None:
        total_wt_aa = prefix_offset_aa + len(consensus_aa)
    altered_total_aa = prefix_offset_aa + altered_stop_local

    pre_start = max(0, diff_idx - flank)

    # WT line: leading-truncation marker, plain residues, trailing-truncation marker.
    wt_post_end = min(len(consensus_aa), consensus_change_end + flank)
    wt_pre_count = prefix_offset_aa + pre_start
    wt_post_count = max(0, total_wt_aa - (prefix_offset_aa + wt_post_end))
    wt_display = (
        (f'[{wt_pre_count}AA]...' if wt_pre_count > 0 else '')
        + consensus_aa[pre_start:wt_post_end]
        + (f'...[{wt_post_count}AA]' if wt_post_count > 0 else '')
        + f' ({total_wt_aa}AA total)'
    )

    # Altered line: leading-truncation marker, prefix identical to WT,
    # [changed region terminating at PTC]. The predicted protein terminates at
    # the PTC, so no trailing residues are shown.
    altered_pre_count = prefix_offset_aa + pre_start
    altered_prefix = altered_aa[pre_start:diff_idx]
    altered_changed = altered_aa[diff_idx:altered_change_end]
    altered_display = (
        (f'[{altered_pre_count}AA]...' if altered_pre_count > 0 else '')
        + altered_prefix
        + '[' + altered_changed + ']'
        + f' ({altered_total_aa}AA total)'
    )

    return (wt_display, altered_display)


def _format_aa_change_and_detect(altered_aa, consensus_aa, frameshift=False):
    """Return (stop_introduced, formatted_aa_change) for an aberration.

    Mirrors spliceai_parser.py:436-518. Handles both in-frame (frameshift=False)
    and frameshift (frameshift=True) cases.

    For frameshift=True with no `*` in the altered translation, returns
    (False, EXTENDS_PAST_NATIVE_STOP) — matches the reference's
    "protein sequence extends beyond native stop site" sentinel.

    Trailing-`*` strip: for in-frame, strip from BOTH proteins so that
    `find('*')` uniquely signals a premature stop (without the strip, an
    in-frame deletion that shifts the natural stop earlier would be
    misclassified as introduced). For frameshift, only the consensus is
    stripped — the altered translation reads the natural stop's 3 bp out of
    frame, so any '*' in altered_aa is necessarily a PTC.
    """
    if consensus_aa.endswith('*'):
        consensus_aa = consensus_aa[:-1]
    if not frameshift and altered_aa.endswith('*'):
        altered_aa = altered_aa[:-1]

    forward_diff = _find_difference_point(altered_aa, consensus_aa, 'forward')
    if forward_diff == 'identical':
        return (False, None)

    # diff_length must be computed BEFORE altered_aa is trimmed, mirroring the
    # reference's :438-:445 ordering. Computing it after the trim flips a
    # net-gain insertion into the net-loss branch and corrupts the output.
    diff_length = len(altered_aa) - len(consensus_aa)

    if isinstance(forward_diff, int) and forward_diff >= 4:
        altered_aa = altered_aa[forward_diff - 4:]

    if not frameshift:
        if diff_length <= 0:
            # Net loss / equal-length divergence.
            reverse_diff = _find_difference_point(altered_aa, consensus_aa, 'reverse')
            if reverse_diff == 'identical':
                return (False, None)
            if isinstance(reverse_diff, int):
                if reverse_diff < 3:
                    altered_aa = altered_aa[:6]
                else:
                    altered_aa = altered_aa[:reverse_diff + 3]
        else:
            # Net gain. Mirrors :459-:484.
            altered_seq_check = altered_aa[3:3 + diff_length] if len(altered_aa) > 3 else ''
            uniq_chars = set(altered_seq_check)
            if len(uniq_chars) == 1 and altered_seq_check:
                altered_seq_check = altered_seq_check[0]
                consensus_seq_check = (
                    consensus_aa[forward_diff - 2:forward_diff - 1]
                    if isinstance(forward_diff, int) else ''
                )
            else:
                start_pos = forward_diff - 2 if isinstance(forward_diff, int) else 0
                end_pos = start_pos + diff_length
                if end_pos > len(consensus_aa):
                    consensus_seq_check = consensus_aa[start_pos:]
                else:
                    consensus_seq_check = consensus_aa[start_pos:end_pos]
            if altered_seq_check == consensus_seq_check:
                altered_aa = altered_aa[:3 + diff_length + 3]
            else:
                reverse_diff = _find_difference_point(altered_aa, consensus_aa, 'reverse')
                if reverse_diff == 'identical':
                    return (False, None)
                if isinstance(reverse_diff, int):
                    if reverse_diff < 3:
                        altered_aa = altered_aa[:6]
                    else:
                        altered_aa = altered_aa[:reverse_diff + 3]

    stop_pos = altered_aa.find('*')

    # Frameshift extension: no PTC reached within the translated CDS-clamped
    # spans. Mirrors spliceai_parser.py:489-490.
    if frameshift and stop_pos == -1:
        return (False, EXTENDS_PAST_NATIVE_STOP)

    stop_introduced = stop_pos != -1
    if stop_introduced:
        altered_aa = altered_aa[:stop_pos + 1]

    pre_out_info = altered_aa
    length_out_info = len(pre_out_info)
    if isinstance(forward_diff, int) and forward_diff <= 3:
        prefix_stop = forward_diff - 1
        suffix_start = prefix_stop
    else:
        prefix_stop = 3
        suffix_start = 3

    if not stop_introduced:
        prefix = pre_out_info[:prefix_stop]
        if length_out_info >= 6:
            suffix = pre_out_info[length_out_info - 3:]
            middle = pre_out_info[suffix_start:length_out_info - 3]
        else:
            suffix = pre_out_info[3:] if length_out_info > 3 else ''
            middle = ''
        formatted = f'{prefix}[{middle}]{suffix}'
    else:
        formatted = f'{pre_out_info[:prefix_stop]}[{pre_out_info[suffix_start:]}]'

    return (stop_introduced, formatted)


def _detect_premature_stop(transcript_scores, aberration, fasta, chrom, ref, alt, var_pos,
                            consensus_ctx=None):
    """Detect whether an aberration introduces a premature stop codon.

    Handles both in-frame and frameshift aberrations (reads
    aberration['frameshift'] to decide). For frameshift cases without a PTC
    in the translated CDS-clamped spans, returns
    (False, EXTENDS_PAST_NATIVE_STOP).

    Returns (stop_introduced, aa_change_string, wt_protein_window,
    altered_protein_window). `stop_introduced` is True/False/None (None =
    detection skipped, e.g. mitochondrial transcript, missing FASTA, indel
    straddling a splice boundary, reference mismatch). The two window strings
    are populated when a PTC is introduced and None otherwise.

    `consensus_ctx`, when provided, is a `_compute_consensus_context` result;
    callers detecting over multiple aberrations on the same transcript should
    compute it once and pass it in. When omitted, it is computed lazily here.
    """
    # Mitochondrial guard: mt transcripts use a non-standard codon table and
    # the chrM/chrMT contig naming differs across FASTAs. Skip rather than
    # mis-translate.
    if chrom is not None:
        chrom_normalized = str(chrom).upper()
        if chrom_normalized.startswith('CHR'):
            chrom_normalized = chrom_normalized[3:]
        if chrom_normalized in ('M', 'MT'):
            return (None, None, None, None)

    if fasta is None or chrom is None or ref is None or alt is None or var_pos is None:
        return (None, None, None, None)

    cds_start = transcript_scores.get('CDS_START')
    cds_end = transcript_scores.get('CDS_END')
    strand = transcript_scores.get('STRAND') or transcript_scores.get('t_strand')

    try:
        if consensus_ctx is None:
            consensus_ctx = _compute_consensus_context(transcript_scores, fasta, chrom)
        if consensus_ctx is None or consensus_ctx.get('consensus_aa_full') is None:
            return (None, None, None, None)

        built = _build_altered_exon_spans(
            transcript_scores, aberration, consensus=consensus_ctx['consensus'],
        )
        if built is None:
            return (None, None, None, None)
        altered_spans, row_anchor, first_eframe = built
        if not altered_spans:
            return (None, None, None, None)

        # eFrame=1 wants to trim 2 leading bases off the first altered span,
        # but `_translate_spans` only does so when the span has >= 2 bases. A
        # 1-bp first altered span at first_eframe=1 silently skips the trim,
        # so the altered tail comes out frame-shifted relative to the
        # consensus slice (which the slice formula assumes was trimmed).
        # Bail just for this aberration. eFrame=2 trims a single base, which
        # is always available since CDS-clamped spans are non-empty.
        if first_eframe == 1 and (altered_spans[0][1] - altered_spans[0][0] + 1) < 2:
            return (None, None, None, None)

        cum = consensus_ctx['cum_cds_bases']
        prefix_offset_aa = (cum[row_anchor] + 2) // 3
        consensus_aa = consensus_ctx['consensus_aa_full'][prefix_offset_aa:]
        consensus_aa_full = consensus_ctx['consensus_aa_full']
        total_wt_aa = (len(consensus_aa_full) - 1
                       if consensus_aa_full.endswith('*') else len(consensus_aa_full))

        altered_aa = _translate_spans(
            altered_spans, first_eframe, fasta, chrom, strand,
            apply_variant_args=(var_pos, ref, alt, cds_start, cds_end),
            seq_cache=consensus_ctx.get('seq_cache'),
        )

        if altered_aa in ('reference mismatch', 'impacts native start or stop site',
                          'variant straddles splice boundary'):
            print(f'WARNING: SAI-10k stop-codon detection skipped: {altered_aa}')
            return (None, None, None, None)

        is_frameshift = (aberration.get('frameshift') is True)
        stop_introduced, formatted = _format_aa_change_and_detect(
            altered_aa, consensus_aa, frameshift=is_frameshift,
        )
        # Only build windows for the cases the popup actually renders (PTC
        # introduced). Skipping non-PTC cases keeps the API surface minimal
        # without precomputing strings the frontend won't show.
        if stop_introduced is True:
            wt_window, altered_window = _build_protein_windows(
                altered_aa, consensus_aa, frameshift=is_frameshift,
                prefix_offset_aa=prefix_offset_aa, total_wt_aa=total_wt_aa,
            )
        else:
            wt_window, altered_window = (None, None)
        return (stop_introduced, formatted, wt_window, altered_window)
    except Exception as exc:
        # Defensive: detection must never break the rest of the SAI-10k response.
        print(f'WARNING: SAI-10k stop-codon detection raised '
              f'{type(exc).__name__}: {exc}')
        return (None, None, None, None)
