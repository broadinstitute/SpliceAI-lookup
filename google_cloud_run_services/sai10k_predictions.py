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
"""


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
    Locate the variant in an exon or intron, with biological (transcript-order)
    numbering on reverse-strand genes.

    Returns dict with:
        region_type: 'exon' or 'intron'
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
            if is_reverse:
                nearest = 'acceptor' if dist_to_upstream_exon < dist_to_downstream_exon else 'donor'
            else:
                nearest = 'donor' if dist_to_upstream_exon < dist_to_downstream_exon else 'acceptor'
            return {
                'region_type': 'intron',
                'region_number': intron_num,
                '_genomic_index': i,
                'distance_to_boundary': min(dist_to_upstream_exon, dist_to_downstream_exon),
                'nearest_boundary': nearest,
            }

    return None


def _get_native_splice_sites(affected_region, variant_pos, exon_starts, exon_ends, strand):
    """
    Return the native splice-site genomic positions and their delta-positions
    (offset from the variant) for the exon "associated" with the variant.

    For an exonic variant, the associated exon is the containing exon.
    For an intronic variant, the associated exon is the CLOSEST flanking exon
    (smaller distance_to_boundary). The reference parser iterates per-exon and
    checks each; this simplification picks one exon per variant.

    Both acceptor and donor come from the SAME exon, so the increased_exon_inclusion
    condition (DP_AG == DP_NA AND DP_DG == DP_ND) is geometrically satisfiable.

    Also returns the adjacent exons' boundaries (prev_eend, prev_estart,
    next_estart, next_eend) in genomic order, used by the cryptic subflow's
    orientation checks.

    Returns dict with:
        geo_na, geo_nd: 1-based genomic positions of native acceptor / donor
                        of the associated exon
        dp_na, dp_nd:   offsets from variant_pos
        native_exon_size: size of the associated native exon (bp)
        assoc_idx:      0-based genomic index of the associated exon
        prev_eend, prev_estart, next_estart, next_eend: adjacent exon
                        boundaries (genomic order), or None at the transcript ends
    """
    if not affected_region or not exon_starts:
        return None

    idx = affected_region['_genomic_index']
    is_reverse = strand == '-'
    num_exons = len(exon_starts)

    if affected_region['region_type'] == 'exon':
        assoc_idx = idx
    else:
        # Intronic: pick the closer flanking exon.
        # Intron `idx` lies between genomic exons `idx` and `idx + 1`.
        intron_start = exon_ends[idx] + 1
        intron_end = exon_starts[idx + 1] - 1
        dist_to_upstream = variant_pos - intron_start   # from exon idx's 3' end
        dist_to_downstream = intron_end - variant_pos   # from exon idx+1's 5' end
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

    prev_estart = exon_starts[assoc_idx - 1] if assoc_idx > 0 else None
    prev_eend = exon_ends[assoc_idx - 1] if assoc_idx > 0 else None
    next_estart = exon_starts[assoc_idx + 1] if assoc_idx + 1 < num_exons else None
    next_eend = exon_ends[assoc_idx + 1] if assoc_idx + 1 < num_exons else None

    return {
        'geo_na': geo_na,
        'geo_nd': geo_nd,
        'dp_na': geo_na - variant_pos,
        'dp_nd': geo_nd - variant_pos,
        'native_exon_size': exon_end - exon_start + 1,
        'assoc_idx': assoc_idx,
        'prev_estart': prev_estart,
        'prev_eend': prev_eend,
        'next_estart': next_estart,
        'next_eend': next_eend,
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
    # positions exactly at exon boundaries; a wider tolerance risks matching
    # the wrong exon in compact gene structures.
    tolerance = 0

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
        if abs(acceptor_pos - geo_al) <= tolerance:
            acceptor_candidates.append(i)
        if abs(donor_pos - geo_dl) <= tolerance:
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


def sai10k_determine_aberrations(
        ds_ag, ds_al, ds_dg, ds_dl,
        dp_ag, dp_al, dp_dg, dp_dl,
        ds_ag_alt, ds_al_alt, ds_dg_alt, ds_dl_alt,
        variant_pos, exon_starts, exon_ends, strand,
):
    """
    Classify splice aberration(s) per Canson et al. 2023 flowchart v1.1.

    Returns a list of aberration dicts (possibly empty). Each dict contains:
        aberration_type: one of the six paper terminal classes
            (exon_skipping, whole_intron_retention, pseudoexon,
             increased_exon_inclusion, partial_exon_deletion,
             partial_intron_retention).
        description: short human-readable description.
        max_delta_score: driving delta score.
        affected_region: dict from sai10k_find_affected_region, or None.
        geo_ag, geo_al, geo_dg, geo_dl: 1-based genomic positions of predicted
            acceptor-gain, acceptor-loss, donor-gain, donor-loss sites.
        _branch: internal tag ('loss', 'gain_A', 'gain_B_acceptor', 'gain_B_donor')
        _native: dict with geo_na/geo_nd/dp_na/dp_nd/native_exon_size (or None)

    The paper flowchart's three branches (LOSS, GAIN Subflow A, GAIN Subflow B)
    evaluate INDEPENDENTLY; a single variant can produce up to one aberration
    from each.
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
        d_from_exon = affected_region['distance_to_boundary']
    else:
        d_from_exon = None  # couldn't locate variant

    def make(ab_type, description, driving_score, branch):
        return {
            'aberration_type': ab_type,
            'description': description,
            'max_delta_score': driving_score,
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
        driving = max(ds_dl, ds_al)
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
            ab = make(
                'exon_skipping',
                'Predicted exon skipping (donor and acceptor loss)',
                driving, 'loss')
        else:
            ab = make(
                'whole_intron_retention',
                'Predicted whole intron retention (donor and acceptor loss)',
                driving, 'loss')
        ab['_potential'] = is_potential
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
        driving = max(ds_dg, ds_ag)

        # Identify which intron AG and DG fall into (same-intron check).
        def intron_of(position):
            for i in range(len(exon_starts) - 1):
                if exon_ends[i] < position < exon_starts[i + 1]:
                    return i
            return None

        ag_intron = intron_of(geo_ag)
        dg_intron = intron_of(geo_dg)

        # Distances from AG and DG to the boundaries of their host intron.
        # The reference parser checks distance from the intron boundaries
        # (i.e. the flanking exon edges), not from all splice sites in the
        # transcript. Using only the host intron avoids false negatives in
        # compact gene structures where a non-flanking exon boundary
        # happens to be within 50 bp.
        if ag_intron is not None:
            ag_dist_to_native = min(
                abs(geo_ag - exon_ends[ag_intron]),
                abs(geo_ag - exon_starts[ag_intron + 1]),
            )
        else:
            ag_dist_to_native = None
        if dg_intron is not None:
            dg_dist_to_native = min(
                abs(geo_dg - exon_ends[dg_intron]),
                abs(geo_dg - exon_starts[dg_intron + 1]),
            )
        else:
            dg_dist_to_native = None

        pseudoexon_conditions = (
            GEX_SIZE_MIN <= gex_size <= GEX_SIZE_MAX
            and d_from_exon is not None
            and d_from_exon > DISTANCE_NATIVE_SITE_MIN_PSEUDOEXON
            and ag_dist_to_native is not None
            and ag_dist_to_native > DISTANCE_NATIVE_SITE_MIN_PSEUDOEXON
            and dg_dist_to_native is not None
            and dg_dist_to_native > DISTANCE_NATIVE_SITE_MIN_PSEUDOEXON
            and ag_intron is not None
            and dg_intron is not None
            and ag_intron == dg_intron
        )

        if pseudoexon_conditions:
            aberrations.append(make(
                'pseudoexon',
                'Predicted pseudoexon activation',
                driving, 'gain_A'))
        elif native and native['native_exon_size'] and gex_size == native['native_exon_size']:
            # Increased exon inclusion: the gained exon exactly matches a
            # native exon AND the gain positions coincide with native sites.
            if dp_ag == native['dp_na'] and dp_dg == native['dp_nd']:
                aberrations.append(make(
                    'increased_exon_inclusion',
                    'Predicted increased exon inclusion',
                    driving, 'gain_A'))

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

        dp_na = native['dp_na']
        dp_nd = native['dp_nd']

        if activate_acceptor and _cryptic_acceptor_orientation_ok(geo_ag, native, is_reverse=strand == '-'):
            # partial_size = (DP_AG*Strand) - (DP_NA*Strand)
            partial_size = (dp_ag - dp_na) * strand_sign
            driving = ds_ag  # use delta score, not raw ALT score
            if (partial_size > 0 and partial_size < native['native_exon_size']
                    and d_from_exon is not None and d_from_exon <= DISTANCE_EXON_MAX_LOSS):
                aberrations.append(make(
                    'partial_exon_deletion',
                    "5' partial exon deletion (cryptic acceptor)",
                    driving, 'gain_B_acceptor'))
                aberrations[-1]['_partial_size'] = partial_size
            elif (partial_size < 0 and partial_size >= -DISTANCE_PARTIAL_RETENTION_MAX
                  and d_from_exon is not None and d_from_exon <= DISTANCE_PARTIAL_RETENTION_MAX):
                aberrations.append(make(
                    'partial_intron_retention',
                    "5' partial intron retention (cryptic acceptor)",
                    driving, 'gain_B_acceptor'))
                # Store signed partial_size (negative = retention); formatters take abs() for display.
                aberrations[-1]['_partial_size'] = partial_size

        if activate_donor and _cryptic_donor_orientation_ok(geo_dg, native, is_reverse=strand == '-'):
            # partial_size = (DP_ND*Strand) - (DP_DG*Strand)
            partial_size = (dp_nd - dp_dg) * strand_sign
            driving = ds_dg  # use delta score, not raw ALT score
            if (partial_size > 0 and partial_size < native['native_exon_size']
                    and d_from_exon is not None and d_from_exon <= DISTANCE_EXON_MAX_LOSS):
                aberrations.append(make(
                    'partial_exon_deletion',
                    "3' partial exon deletion (cryptic donor)",
                    driving, 'gain_B_donor'))
                aberrations[-1]['_partial_size'] = partial_size
            elif (partial_size < 0 and partial_size >= -DISTANCE_PARTIAL_RETENTION_MAX
                  and d_from_exon is not None and d_from_exon <= DISTANCE_PARTIAL_RETENTION_MAX):
                aberrations.append(make(
                    'partial_intron_retention',
                    "3' partial intron retention (cryptic donor)",
                    driving, 'gain_B_donor'))
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
    """
    if native is None:
        return True
    if is_reverse:
        # Biologically "prev" = genomically "next" for - strand.
        bio_prev_estart = native.get('next_estart')
        if bio_prev_estart is None:
            return True
        return geo_ag - bio_prev_estart < 0
    else:
        prev_eend = native.get('prev_eend')
        if prev_eend is None:
            return True
        return prev_eend - geo_ag < 0


def _cryptic_donor_orientation_ok(geo_dg, native, is_reverse):
    """
    Reference R parser (spliceAI_parser.R:954):
        + strand: GEO_DG - next_eStart < 0    (GEO_DG before next exon's acceptor)
        - strand: next_eEnd - GEO_DG < 0       (in R biological order)

    Same prev/next swap as _cryptic_acceptor_orientation_ok.
    """
    if native is None:
        return True
    if is_reverse:
        # Biologically "next" = genomically "prev" for - strand.
        bio_next_eend = native.get('prev_eend')
        if bio_next_eend is None:
            return True
        return bio_next_eend - geo_dg < 0
    else:
        next_estart = native.get('next_estart')
        if next_estart is None:
            return True
        return geo_dg - next_estart < 0


def sai10k_annotate_frameshift(aberration, cds_start, cds_end, exon_starts, exon_ends):
    """
    Mutate `aberration` in place with size / frameshift / coding fields.

    Adds:
        affects_coding: True | False | None
        frameshift: True | False | None
        size: integer bp when determinable
        size_type: 'exon' | 'intron' | 'pseudoexon' | 'partial'
        frameshift_description: formatted string for display
        exon_numbers: list of biological exon numbers (for exon_skipping)
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

    # Helper: does a closed segment [a, b] overlap the CDS?
    # Returns False for non-coding transcripts (so affects_coding lands at False).
    def in_cds(a, b):
        if not is_coding_transcript:
            return False
        if a is None or b is None:
            return None
        lo = min(a, b)
        hi = max(a, b)
        return not (hi < cds_start or lo > cds_end)

    is_potential = aberration.get('_potential', False)

    if aberration_type == 'exon_skipping':
        strand = aberration.get('_strand', '+')
        result = _find_exons_between_loss_positions(geo_al, geo_dl, strand, exon_starts, exon_ends)
        # Reference-parser fallback: the Canson algorithm (spliceAI_parser.R:239,
        # spliceai_parser.py:312-329) requires BOTH geo_al and geo_dl to sit exactly
        # at exon boundaries. When SpliceAI's DP_AL/DP_DL peak lands interior to the
        # affected exon (e.g. donor-site variants where acceptor-loss peaks inside the
        # same exon), that exact-match lookup returns NA|NA. If the variant is exonic
        # and at least one of geo_al/geo_dl still matches the affected exon's own
        # splice site, treat that single exon as the skipped one.
        if not result['match'] and affected_region and affected_region.get('region_type') == 'exon':
            is_reverse = strand == '-'
            idx = affected_region.get('_genomic_index')
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
            total_size = result['total_size']
            frameshift = total_size % 3 != 0
            segment_in_cds = any(
                in_cds(exon_starts[i], exon_ends[i])
                for i in result['genomic_indices']
            )
            aberration['affects_coding'] = bool(segment_in_cds)
            aberration['frameshift'] = frameshift if segment_in_cds else None
            aberration['size'] = total_size
            aberration['size_type'] = 'exon'
            aberration['exon_numbers'] = result['biological_numbers']
            if is_potential:
                label = 'Potential exon' if len(result['biological_numbers']) == 1 else 'Potential exons'
            else:
                label = 'Exon' if len(result['biological_numbers']) == 1 else 'Exons'
            nums = ', '.join(str(n) for n in result['biological_numbers'])
            aberration['frameshift_description'] = (
                f'{label} {nums} skipping ({total_size}bp) - '
                f'{"frameshift" if frameshift and segment_in_cds else "in-frame" if segment_in_cds else "non-coding"}')
        else:
            aberration['affects_coding'] = None
            aberration['frameshift'] = None
            aberration['frameshift_description'] = 'Exon skipping predicted but skipped exon(s) could not be mapped'
        return

    if aberration_type == 'whole_intron_retention':
        if geo_al is not None and geo_dl is not None:
            intron_size = abs(geo_al - geo_dl) - 1
            # Segment = retained intron: from min(geo_al, geo_dl)+1 to max(...)-1
            seg_lo = min(geo_al, geo_dl) + 1
            seg_hi = max(geo_al, geo_dl) - 1
            segment_in_cds = in_cds(seg_lo, seg_hi)
            frameshift = intron_size % 3 != 0
            aberration['affects_coding'] = bool(segment_in_cds) if segment_in_cds is not None else None
            aberration['frameshift'] = frameshift if segment_in_cds else None
            aberration['size'] = intron_size
            aberration['size_type'] = 'intron'
            label = 'Potential whole intron retention' if is_potential else 'Whole intron retention'
            aberration['frameshift_description'] = (
                f'{label} ({intron_size}bp) - '
                f'{"frameshift" if frameshift and segment_in_cds else "in-frame" if segment_in_cds else "non-coding"}')
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
            segment_in_cds = in_cds(seg_lo, seg_hi)
            frameshift = gex_size % 3 != 0
            aberration['affects_coding'] = bool(segment_in_cds) if segment_in_cds is not None else None
            aberration['frameshift'] = frameshift if segment_in_cds else None
            aberration['size'] = gex_size
            aberration['size_type'] = 'pseudoexon'
            aberration['frameshift_description'] = (
                f'Pseudoexon activation ({gex_size}bp) - '
                f'{"frameshift" if frameshift and segment_in_cds else "in-frame" if segment_in_cds else "non-coding"}')
        else:
            aberration['affects_coding'] = None
            aberration['frameshift'] = None
            aberration['frameshift_description'] = 'Pseudoexon activation - size unclear'
        return

    if aberration_type == 'increased_exon_inclusion':
        native = aberration.get('_native') or {}
        native_exon_size = native.get('native_exon_size')
        if native_exon_size is not None:
            aberration['size'] = native_exon_size
            aberration['size_type'] = 'exon'
            aberration['frameshift'] = None  # same exon, no frame change
            # Check coding by the native exon coords (use the associated exon index,
            # not the intron index, since for intronic variants these can differ).
            idx = native.get('assoc_idx')
            if idx is not None and 0 <= idx < len(exon_starts):
                segment_in_cds = in_cds(exon_starts[idx], exon_ends[idx])
            else:
                segment_in_cds = None
            aberration['affects_coding'] = bool(segment_in_cds) if segment_in_cds is not None else None
            aberration['frameshift_description'] = (
                f'Increased exon inclusion ({native_exon_size}bp) - '
                f'{"coding" if segment_in_cds else "non-coding" if segment_in_cds is False else "unknown"}')
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
        display_size = abs(partial_size)
        frameshift = display_size % 3 != 0
        # Segment = range between native site and predicted cryptic site.
        native = aberration.get('_native') or {}
        branch = aberration.get('_branch', '')
        if branch == 'gain_B_acceptor':
            seg_a = native.get('geo_na')
            seg_b = geo_ag
            shift_type = 'Acceptor'
        else:
            seg_a = native.get('geo_nd')
            seg_b = geo_dg
            shift_type = 'Donor'
        segment_in_cds = in_cds(seg_a, seg_b) if seg_a is not None and seg_b is not None else None
        aberration['affects_coding'] = bool(segment_in_cds) if segment_in_cds is not None else None
        aberration['frameshift'] = frameshift if segment_in_cds else None
        aberration['size'] = display_size
        aberration['size_type'] = 'partial'
        consequence = 'partial exon deletion' if aberration_type == 'partial_exon_deletion' else 'partial intron retention'
        aberration['frameshift_description'] = (
            f'{shift_type} shift ({display_size}bp {consequence}) - '
            f'{"frameshift" if frameshift and segment_in_cds else "in-frame" if segment_in_cds else "non-coding"}')
        return

    # Unknown aberration type — leave blank (should not happen with current paper-only emission).
    aberration['affects_coding'] = None
    aberration['frameshift'] = None
    aberration['frameshift_description'] = f'{aberration_type} - frameshift status unclear'


def sai10k_compute_predictions(transcript_scores, variant_pos):
    """
    Compute splice aberration predictions for a single transcript.

    Args:
        transcript_scores: dict with SpliceAI delta scores, delta positions,
            raw alt-allele site scores, and transcript structure.
        variant_pos: 1-based genomic position of the variant.

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

    ds_ag_alt = transcript_scores.get('DS_AG_ALT', 0)
    ds_al_alt = transcript_scores.get('DS_AL_ALT', 0)
    ds_dg_alt = transcript_scores.get('DS_DG_ALT', 0)
    ds_dl_alt = transcript_scores.get('DS_DL_ALT', 0)

    exon_starts = transcript_scores.get('EXON_STARTS', [])
    exon_ends = transcript_scores.get('EXON_ENDS', [])
    cds_start = transcript_scores.get('CDS_START')
    cds_end = transcript_scores.get('CDS_END')
    strand = transcript_scores.get('STRAND') or transcript_scores.get('t_strand', '+')

    aberrations = sai10k_determine_aberrations(
        ds_ag, ds_al, ds_dg, ds_dl,
        dp_ag, dp_al, dp_dg, dp_dl,
        ds_ag_alt, ds_al_alt, ds_dg_alt, ds_dl_alt,
        variant_pos, exon_starts, exon_ends, strand,
    )

    # Compute the overall max Δ score, the corresponding Δ type, and a
    # confidence label. These are constant per variant (not per aberration) and
    # drive the line-1 display ("Max Δ score: X — donor loss (high confidence)").
    # Reporting the per-aberration "driving" score on its own can mislead readers
    # when, e.g., partial_exon_deletion fires off a 0.20 cryptic gain while the
    # variant also has a 0.92 native donor loss — the user should still see 0.92
    # as the headline confidence (Canson reviewer feedback, 2026-04-17).
    overall_max_score, overall_delta_type = max(
        ((_coerce_float(ds_al), 'acceptor loss'),
         (_coerce_float(ds_dl), 'donor loss'),
         (_coerce_float(ds_ag), 'acceptor gain'),
         (_coerce_float(ds_dg), 'donor gain')),
        key=lambda item: item[0],
    )
    if overall_max_score >= 0.8:
        confidence = 'high'
    elif overall_max_score >= 0.5:
        confidence = 'moderate'
    else:
        confidence = 'low'

    for aberration in aberrations:
        sai10k_annotate_frameshift(
            aberration,
            cds_start, cds_end,
            exon_starts, exon_ends,
        )
        # Overwrite the per-branch driving score with the variant-level overall
        # max so all aberration rows share the same line-1 confidence headline.
        aberration['max_delta_score'] = overall_max_score
        aberration['delta_type'] = overall_delta_type
        aberration['confidence'] = confidence

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


def sai10k_get_transcript_predictions(transcript_scores, variant_pos):
    """Compute predictions for one transcript and attach its metadata."""
    if not transcript_scores:
        return None

    predictions = sai10k_compute_predictions(transcript_scores, variant_pos)
    predictions['transcript_id'] = transcript_scores.get('t_id') or transcript_scores.get('NAME')
    predictions['transcript_priority'] = transcript_scores.get('t_priority', 'N')
    predictions['gene_name'] = transcript_scores.get('g_name')
    return predictions
