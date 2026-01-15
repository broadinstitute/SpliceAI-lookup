"""
SAI-10k predictions module for predicting splice consequences.

This module implements heuristics from SAI-10k-calc to predict
the likely biological consequences of splice-altering variants based
on SpliceAI delta scores and transcript structure.
"""


# Minimum delta score to consider significant
MIN_DELTA_SCORE = 0.2

# Score thresholds for different confidence levels
SCORE_THRESHOLD_HIGH = 0.8
SCORE_THRESHOLD_MEDIUM = 0.5
SCORE_THRESHOLD_LOW = 0.2


def sai10k_find_affected_region(variant_pos, exon_starts, exon_ends):
    """
    Determine which exon or intron contains the variant.

    Args:
        variant_pos: 1-based genomic position of the variant
        exon_starts: List of 1-based exon start positions
        exon_ends: List of 1-based exon end positions

    Returns:
        dict with keys:
            'region_type': 'exon' or 'intron'
            'region_number': 1-based exon or intron number
            'distance_to_boundary': Distance to nearest exon boundary
            'nearest_boundary': 'donor' or 'acceptor'
    """
    if not exon_starts or not exon_ends:
        return None

    num_exons = len(exon_starts)

    # Check each exon
    for i in range(num_exons):
        e_start = exon_starts[i]
        e_end = exon_ends[i]

        if e_start <= variant_pos <= e_end:
            # Variant is in this exon
            dist_to_start = variant_pos - e_start
            dist_to_end = e_end - variant_pos

            return {
                'region_type': 'exon',
                'region_number': i + 1,
                'distance_to_boundary': min(dist_to_start, dist_to_end),
                'nearest_boundary': 'acceptor' if dist_to_start < dist_to_end else 'donor',
            }

    # Check each intron (between exons)
    for i in range(num_exons - 1):
        intron_start = exon_ends[i] + 1
        intron_end = exon_starts[i + 1] - 1

        if intron_start <= variant_pos <= intron_end:
            dist_to_donor = variant_pos - intron_start
            dist_to_acceptor = intron_end - variant_pos

            return {
                'region_type': 'intron',
                'region_number': i + 1,
                'distance_to_boundary': min(dist_to_donor, dist_to_acceptor),
                'nearest_boundary': 'donor' if dist_to_donor < dist_to_acceptor else 'acceptor',
            }

    return None


def sai10k_determine_aberration_type(ds_ag, ds_al, ds_dg, ds_dl, dp_ag, dp_al, dp_dg, dp_dl,
                                     variant_pos, exon_starts, exon_ends, strand):
    """
    Determine the likely splice aberration type based on delta scores.

    Args:
        ds_ag, ds_al, ds_dg, ds_dl: Delta scores for acceptor/donor gain/loss
        dp_ag, dp_al, dp_dg, dp_dl: Delta positions for acceptor/donor gain/loss
        variant_pos: 1-based genomic position
        exon_starts, exon_ends: Lists of exon boundaries
        strand: '+' or '-'

    Returns:
        dict with predicted aberration info
    """
    # Convert scores to float
    ds_ag = float(ds_ag) if ds_ag else 0
    ds_al = float(ds_al) if ds_al else 0
    ds_dg = float(ds_dg) if ds_dg else 0
    ds_dl = float(ds_dl) if ds_dl else 0

    # Find max delta score
    max_score = max(ds_ag, ds_al, ds_dg, ds_dl)
    if max_score < MIN_DELTA_SCORE:
        return {
            'aberration_type': 'none',
            'confidence': 'low',
            'description': 'No significant splice effect predicted',
        }

    # Determine dominant effect
    affected_region = sai10k_find_affected_region(variant_pos, exon_starts, exon_ends)

    # Calculate genomic positions of predicted splice changes
    geo_ag = variant_pos + dp_ag if dp_ag else None
    geo_al = variant_pos + dp_al if dp_al else None
    geo_dg = variant_pos + dp_dg if dp_dg else None
    geo_dl = variant_pos + dp_dl if dp_dl else None

    result = {
        'max_delta_score': max_score,
        'affected_region': affected_region,
        'geo_ag': geo_ag,
        'geo_al': geo_al,
        'geo_dg': geo_dg,
        'geo_dl': geo_dl,
    }

    # Determine confidence level
    if max_score >= SCORE_THRESHOLD_HIGH:
        result['confidence'] = 'high'
    elif max_score >= SCORE_THRESHOLD_MEDIUM:
        result['confidence'] = 'medium'
    else:
        result['confidence'] = 'low'

    # Classify aberration type based on score patterns
    if ds_dl >= MIN_DELTA_SCORE and ds_al >= MIN_DELTA_SCORE:
        # Donor loss + Acceptor loss = likely exon skipping
        result['aberration_type'] = 'exon_skipping'
        result['description'] = 'Predicted exon skipping (donor and acceptor loss)'

    elif ds_dg >= MIN_DELTA_SCORE and ds_ag >= MIN_DELTA_SCORE:
        # Donor gain + Acceptor gain = likely pseudoexon activation
        if affected_region and affected_region['region_type'] == 'intron':
            result['aberration_type'] = 'pseudoexon'
            result['description'] = 'Predicted pseudoexon activation (new splice sites in intron)'
        else:
            result['aberration_type'] = 'cryptic_splice_site'
            result['description'] = 'Predicted cryptic splice site activation'

    elif ds_dl >= MIN_DELTA_SCORE and ds_dg >= MIN_DELTA_SCORE:
        # Donor loss + Donor gain = donor site shift
        result['aberration_type'] = 'donor_shift'
        if geo_dg and geo_dl:
            shift = geo_dg - geo_dl
            if shift > 0:
                result['description'] = f'Predicted donor site shift ({abs(shift)}bp downstream)'
            else:
                result['description'] = f'Predicted donor site shift ({abs(shift)}bp upstream)'
        else:
            result['description'] = 'Predicted donor site shift'

    elif ds_al >= MIN_DELTA_SCORE and ds_ag >= MIN_DELTA_SCORE:
        # Acceptor loss + Acceptor gain = acceptor site shift
        result['aberration_type'] = 'acceptor_shift'
        if geo_ag and geo_al:
            shift = geo_ag - geo_al
            if shift > 0:
                result['description'] = f'Predicted acceptor site shift ({abs(shift)}bp downstream)'
            else:
                result['description'] = f'Predicted acceptor site shift ({abs(shift)}bp upstream)'
        else:
            result['description'] = 'Predicted acceptor site shift'

    elif ds_dl >= MIN_DELTA_SCORE:
        # Donor loss only - could be intron retention or exon skipping
        result['aberration_type'] = 'donor_loss'
        result['description'] = 'Predicted donor site loss (may cause intron retention or exon skipping)'

    elif ds_al >= MIN_DELTA_SCORE:
        # Acceptor loss only
        result['aberration_type'] = 'acceptor_loss'
        result['description'] = 'Predicted acceptor site loss (may cause intron retention or exon skipping)'

    elif ds_dg >= MIN_DELTA_SCORE:
        # Donor gain only
        result['aberration_type'] = 'donor_gain'
        result['description'] = 'Predicted cryptic donor site activation'

    elif ds_ag >= MIN_DELTA_SCORE:
        # Acceptor gain only
        result['aberration_type'] = 'acceptor_gain'
        result['description'] = 'Predicted cryptic acceptor site activation'

    else:
        result['aberration_type'] = 'unknown'
        result['description'] = 'Splice effect predicted but type unclear'

    return result


def sai10k_predict_frameshift(aberration_type, aberration_info, cds_start, cds_end,
                              exon_starts, exon_ends, exon_frames):
    """
    Predict whether the splice aberration causes a frameshift.

    Args:
        aberration_type: Type of aberration from determine_aberration_type()
        aberration_info: Full info dict from determine_aberration_type()
        cds_start, cds_end: CDS boundaries (1-based)
        exon_starts, exon_ends: Exon boundaries (1-based)
        exon_frames: Frame values for each exon (-1 for non-coding)

    Returns:
        dict with frameshift prediction
    """
    if not cds_start or not cds_end:
        return {
            'affects_coding': False,
            'frameshift': None,
            'description': 'Non-coding transcript',
        }

    if not exon_frames:
        return {
            'affects_coding': None,
            'frameshift': None,
            'description': 'Exon frame information not available',
        }

    affected_region = aberration_info.get('affected_region')
    if not affected_region:
        return {
            'affects_coding': None,
            'frameshift': None,
            'description': 'Could not determine affected region',
        }

    # Check if the affected region overlaps with CDS
    region_type = affected_region['region_type']
    region_num = affected_region['region_number']

    if region_type == 'exon':
        exon_idx = region_num - 1
        if exon_idx < len(exon_starts):
            exon_start = exon_starts[exon_idx]
            exon_end = exon_ends[exon_idx]

            # Check if exon overlaps CDS
            if exon_end < cds_start or exon_start > cds_end:
                return {
                    'affects_coding': False,
                    'frameshift': None,
                    'description': 'Affected exon is outside coding region',
                }

    elif region_type == 'intron':
        intron_idx = region_num - 1
        if intron_idx < len(exon_starts) - 1:
            intron_start = exon_ends[intron_idx] + 1
            intron_end = exon_starts[intron_idx + 1] - 1

            # Check if intron is within CDS region
            if intron_end < cds_start or intron_start > cds_end:
                return {
                    'affects_coding': False,
                    'frameshift': None,
                    'description': 'Affected intron is outside coding region',
                }

    # Predict frameshift based on aberration type
    result = {
        'affects_coding': True,
    }

    if aberration_type == 'exon_skipping':
        # Calculate size of potentially skipped exon(s)
        if region_type == 'exon':
            exon_idx = region_num - 1
            if exon_idx < len(exon_starts):
                exon_size = exon_ends[exon_idx] - exon_starts[exon_idx] + 1
                frameshift = exon_size % 3 != 0
                result['frameshift'] = frameshift
                result['exon_size'] = exon_size
                result['description'] = f'Exon {region_num} skipping ({exon_size}bp) - {"frameshift" if frameshift else "in-frame"}'
        else:
            result['frameshift'] = None
            result['description'] = 'Exon skipping predicted but affected exon unclear'

    elif aberration_type == 'donor_loss' or aberration_type == 'acceptor_loss':
        # Could be intron retention
        if region_type == 'intron' or (region_type == 'exon' and aberration_type == 'donor_loss'):
            intron_idx = region_num - 1 if region_type == 'intron' else region_num - 1
            if intron_idx < len(exon_starts) - 1 and intron_idx >= 0:
                intron_size = exon_starts[intron_idx + 1] - exon_ends[intron_idx] - 1
                frameshift = intron_size % 3 != 0
                result['frameshift'] = frameshift
                result['intron_size'] = intron_size
                result['description'] = f'Potential intron {intron_idx + 1} retention ({intron_size}bp) - {"frameshift" if frameshift else "in-frame"}'
            else:
                result['frameshift'] = None
                result['description'] = 'Splice site loss - frameshift status unclear'
        else:
            result['frameshift'] = None
            result['description'] = 'Splice site loss - frameshift status unclear'

    elif aberration_type in ('donor_shift', 'acceptor_shift'):
        # Calculate shift size
        if aberration_type == 'donor_shift':
            geo_new = aberration_info.get('geo_dg')
            geo_old = aberration_info.get('geo_dl')
        else:
            geo_new = aberration_info.get('geo_ag')
            geo_old = aberration_info.get('geo_al')

        if geo_new and geo_old:
            shift_size = abs(geo_new - geo_old)
            frameshift = shift_size % 3 != 0
            result['frameshift'] = frameshift
            result['shift_size'] = shift_size
            result['description'] = f'Splice site shift ({shift_size}bp) - {"frameshift" if frameshift else "in-frame"}'
        else:
            result['frameshift'] = None
            result['description'] = 'Splice site shift - size unclear'

    elif aberration_type == 'pseudoexon':
        # Pseudoexon size from gain positions
        geo_ag = aberration_info.get('geo_ag')
        geo_dg = aberration_info.get('geo_dg')
        if geo_ag and geo_dg:
            pseudo_size = abs(geo_dg - geo_ag)
            frameshift = pseudo_size % 3 != 0
            result['frameshift'] = frameshift
            result['pseudoexon_size'] = pseudo_size
            result['description'] = f'Pseudoexon activation ({pseudo_size}bp) - {"frameshift" if frameshift else "in-frame"}'
        else:
            result['frameshift'] = None
            result['description'] = 'Pseudoexon activation - size unclear'

    else:
        result['frameshift'] = None
        result['description'] = f'{aberration_type} - frameshift status unclear'

    return result


def sai10k_compute_predictions(transcript_scores, variant_pos):
    """
    Compute heuristic predictions for a single transcript.

    This is the main entry point for computing splice heuristics.

    Args:
        transcript_scores: Dict containing SpliceAI scores and transcript info:
            - DS_AG, DS_AL, DS_DG, DS_DL: Delta scores
            - DP_AG, DP_AL, DP_DG, DP_DL: Delta positions
            - TX_START, TX_END: Transcript boundaries
            - CDS_START, CDS_END: CDS boundaries (optional)
            - EXON_STARTS, EXON_ENDS: Exon boundaries
            - EXON_FRAMES: Exon frames (optional)
            - STRAND or t_strand: Strand
        variant_pos: 1-based genomic position of the variant

    Returns:
        dict with heuristic predictions
    """
    # Extract scores
    ds_ag = transcript_scores.get('DS_AG', 0)
    ds_al = transcript_scores.get('DS_AL', 0)
    ds_dg = transcript_scores.get('DS_DG', 0)
    ds_dl = transcript_scores.get('DS_DL', 0)

    dp_ag = transcript_scores.get('DP_AG', 0)
    dp_al = transcript_scores.get('DP_AL', 0)
    dp_dg = transcript_scores.get('DP_DG', 0)
    dp_dl = transcript_scores.get('DP_DL', 0)

    # Extract transcript structure
    exon_starts = transcript_scores.get('EXON_STARTS', [])
    exon_ends = transcript_scores.get('EXON_ENDS', [])
    cds_start = transcript_scores.get('CDS_START')
    cds_end = transcript_scores.get('CDS_END')
    exon_frames = transcript_scores.get('EXON_FRAMES')
    strand = transcript_scores.get('STRAND') or transcript_scores.get('t_strand', '+')

    # Determine aberration type
    aberration_info = sai10k_determine_aberration_type(
        ds_ag, ds_al, ds_dg, ds_dl,
        dp_ag, dp_al, dp_dg, dp_dl,
        variant_pos, exon_starts, exon_ends, strand
    )

    # Predict frameshift
    frameshift_info = sai10k_predict_frameshift(
        aberration_info.get('aberration_type', 'unknown'),
        aberration_info,
        cds_start, cds_end,
        exon_starts, exon_ends, exon_frames
    )

    return {
        'aberration': aberration_info,
        'frameshift': frameshift_info,
        'transcript_info': {
            'strand': strand,
            'tx_start': transcript_scores.get('TX_START'),
            'tx_end': transcript_scores.get('TX_END'),
            'cds_start': cds_start,
            'cds_end': cds_end,
            'num_exons': len(exon_starts) if exon_starts else 0,
            'is_coding': cds_start is not None and cds_end is not None,
        }
    }


def sai10k_select_canonical_transcript(scores_list):
    """
    Select the canonical (highest priority) transcript from a list.

    Selection criteria:
    1. Highest transcript priority (MS > MP > C > N)
    2. Among equal priority, highest max delta score

    Args:
        scores_list: List of transcript score dicts from SpliceAI

    Returns:
        The selected transcript dict, or None if scores_list is empty
    """
    if not scores_list:
        return None

    # Priority order: MS (MANE Select) > MP (MANE Plus Clinical) > C (Canonical) > N (None)
    priority_order = {'MS': 3, 'MP': 2, 'C': 1, 'N': 0}

    best_transcript = None
    best_priority = -1
    best_max_score = 0

    for transcript_scores in scores_list:
        priority = priority_order.get(transcript_scores.get('t_priority', 'N'), 0)

        # Calculate max delta score
        max_score = max(
            abs(float(transcript_scores.get('DS_AG', 0) or 0)),
            abs(float(transcript_scores.get('DS_AL', 0) or 0)),
            abs(float(transcript_scores.get('DS_DG', 0) or 0)),
            abs(float(transcript_scores.get('DS_DL', 0) or 0)),
        )

        # Select highest priority, or highest score if same priority
        if priority > best_priority or (priority == best_priority and max_score > best_max_score):
            best_transcript = transcript_scores
            best_priority = priority
            best_max_score = max_score

    return best_transcript


def sai10k_get_transcript_predictions(transcript_scores, variant_pos):
    """
    Get heuristic predictions for a specific transcript.

    This function computes predictions and adds transcript metadata to the result.

    Args:
        transcript_scores: Dict containing SpliceAI scores and transcript info
        variant_pos: 1-based genomic position of the variant

    Returns:
        dict with predictions including transcript metadata
    """
    if not transcript_scores:
        return None

    predictions = sai10k_compute_predictions(transcript_scores, variant_pos)
    predictions['transcript_id'] = transcript_scores.get('t_id') or transcript_scores.get('NAME')
    predictions['transcript_priority'] = transcript_scores.get('t_priority', 'N')
    predictions['gene_name'] = transcript_scores.get('g_name')

    return predictions


def sai10k_get_canonical_transcript_predictions(scores_list, variant_pos):
    """
    Get heuristic predictions for the canonical (highest priority) transcript.

    This is a convenience function that combines transcript selection and
    prediction computation.

    Args:
        scores_list: List of transcript score dicts from SpliceAI
        variant_pos: 1-based genomic position

    Returns:
        dict with predictions for the canonical transcript, or None if no scores
    """
    best_transcript = sai10k_select_canonical_transcript(scores_list)
    return sai10k_get_transcript_predictions(best_transcript, variant_pos)
