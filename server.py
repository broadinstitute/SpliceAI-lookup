from collections import defaultdict
from datetime import datetime
import json
import markdown2
import os
import pandas as pd
import re
import socket
import subprocess
import traceback
import tempfile
import time

# pangolin imports
from pkg_resources import resource_filename
from pangolin.model import torch, Pangolin, L, W, AR
from pangolin.pangolin import process_variant as process_variant_using_pangolin
import gffutils

# flask imports
from flask import Flask, request, Response, send_from_directory
from flask_cors import CORS
from flask_talisman import Talisman
from intervaltree import IntervalTree, Interval
from spliceai.utils import Annotator, get_delta_scores

# pandas output options
pd.options.display.float_format = "{:,.2f}".format
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('max_colwidth', None)


app = Flask(__name__)

CORS(app)

DEBUG = False if socket.gethostname() == "spliceai-lookup" else True
if not DEBUG:
    Talisman(app)


RATE_LIMIT_WINDOW_SIZE_IN_MINUTES = 1
RATE_LIMIT_REQUESTS_PER_USER_PER_MINUTE = {
    "spliceai:model": 4,
    "spliceai:total": 15,
    "pangolin:model": 4,
    "pangolin:total": 15,
    "liftover:total": 12,
}

RATE_LIMIT_COUNTER_WINDOW_SIZE_IN_DAYS = 3
RATE_LIMIT_OUTLIER_IPS_PATH = os.path.abspath("rate_limit_outlier_ips.txt")

def get_rate_limit_outlier_ips():
    print(f"Reading rate limit outlier IPs: {RATE_LIMIT_OUTLIER_IPS_PATH}")
    if os.path.isfile(RATE_LIMIT_OUTLIER_IPS_PATH):
        with open(RATE_LIMIT_OUTLIER_IPS_PATH, "rt") as f:
            rate_limit_outlier_ips = [l.strip() for l in f]
    else:
        rate_limit_outlier_ips = []

    print(f"Current list of rate limit outlier IPs: {rate_limit_outlier_ips}")
    return rate_limit_outlier_ips


RATE_LIMIT_OUTLIER_IPS = get_rate_limit_outlier_ips()

DISABLE_LOGGING_FOR_IPS = {f"63.143.42.{i}" for i in range(0, 256)}  # ignore uptimerobot.com IPs


HG19_FASTA_PATH = os.path.expanduser("~/hg19.fa")
HG38_FASTA_PATH = os.path.expanduser("~/hg38.fa")
T2T_FASTA_PATH = os.path.expanduser("~/chm13v2.0.fa")

GENCODE_VERSION = "v44"
SPLICEAI_GRCH37_ANNOTATIONS = f"./annotations/gencode.{GENCODE_VERSION}lift37.basic.annotation.txt.gz"
SPLICEAI_GRCH38_ANNOTATIONS = f"./annotations/gencode.{GENCODE_VERSION}.basic.annotation.txt.gz"
PANGOLIN_GRCH37_ANNOTATIONS = f"./annotations/gencode.{GENCODE_VERSION}lift37.basic.annotation.without_chr_prefix.db"
PANGOLIN_GRCH38_ANNOTATIONS = f"./annotations/gencode.{GENCODE_VERSION}.basic.annotation.db"
TRANSCRIPT_GRCH37_ANNOTATIONS = f"./annotations/gencode.{GENCODE_VERSION}lift37.basic.annotation.transcript_annotations.json"
TRANSCRIPT_GRCH38_ANNOTATIONS = f"./annotations/gencode.{GENCODE_VERSION}.basic.annotation.transcript_annotations.json"

UCSC_LIFTOVER_TOOL = "UCSC liftover tool"
BCFTOOLS_LIFTOVER_TOOL = "bcftools liftover plugin"

PANGOLIN_MODELS = []
for i in 0, 2, 4, 6:
    for j in 1, 2, 3:
        model = Pangolin(L, W, AR)
        if torch.cuda.is_available():
            model.cuda()
            weights = torch.load(resource_filename("pangolin", "models/final.%s.%s.3.v2" % (j, i)))
        else:
            weights = torch.load(resource_filename("pangolin", "models/final.%s.%s.3.v2" % (j, i)), map_location=torch.device('cpu'))
        model.load_state_dict(weights)
        model.eval()
        PANGOLIN_MODELS.append(model)


ANNOTATION_INTERVAL_TREES = {
    "37": defaultdict(IntervalTree),
    "38": defaultdict(IntervalTree),
}

for genome_version, annotation_path in ("37", SPLICEAI_GRCH37_ANNOTATIONS), ("38", SPLICEAI_GRCH38_ANNOTATIONS):
    print(f"Loading {annotation_path}", flush=True)
    df = pd.read_table(annotation_path, dtype={"TX_START": int, "TX_END": int})
    for _, row in df.iterrows():
        chrom = row["CHROM"].replace("chr", "")
        ANNOTATION_INTERVAL_TREES[genome_version][chrom].add(Interval(row["TX_START"], row["TX_END"] + 0.1, row["#NAME"]))

SPLICEAI_ANNOTATOR = {
    "37": Annotator(HG19_FASTA_PATH, SPLICEAI_GRCH37_ANNOTATIONS),
    "38": Annotator(HG38_FASTA_PATH, SPLICEAI_GRCH38_ANNOTATIONS),
}

ta37_f = open(TRANSCRIPT_GRCH37_ANNOTATIONS, "rt")
ta38_f = open(TRANSCRIPT_GRCH38_ANNOTATIONS, "rt")
TRANSCRIPT_ANNOTATIONS = {
    "37": json.load(ta37_f),
    "38": json.load(ta38_f),
}
ta37_f.close()
ta38_f.close()

TRANSCRIPT_PRIORITY_ORDER = {
    "MS": 3,  # MANE select transcript
    "MP": 2,  # MANE plus clinical transcript
    "C": 1,   # canonical transcript
    "N": 0
}

# check that json annotations exist for all transcripts in the SpliceAI annotations file
for genome_version in "37", "38":
    json_transcript_ids = set(TRANSCRIPT_ANNOTATIONS[genome_version])
    df = pd.read_table(SPLICEAI_GRCH37_ANNOTATIONS if genome_version == "37" else SPLICEAI_GRCH38_ANNOTATIONS)
    spliceai_annotation_transcript_ids = set(df["#NAME"].apply(lambda t: t.split(".")[0]))
    transcript_ids_without_annotations = spliceai_annotation_transcript_ids - json_transcript_ids
    if len(transcript_ids_without_annotations) > 0:
        raise ValueError(f"Missing {len(transcript_ids_without_annotations)} transcripts in {genome_version} annotations: {transcript_ids_without_annotations}")

SPLICEAI_MAX_DISTANCE_LIMIT = 10000
SPLICEAI_DEFAULT_DISTANCE = 500  # maximum distance between the variant and gained/lost splice site, defaults to 500
SPLICEAI_DEFAULT_MASK = 0        # mask scores representing annotated acceptor/donor gain and unannotated acceptor/donor loss, defaults to 0

SPLICEAI_EXAMPLE = f"/spliceai/?hg=38&distance=500&mask=0&variant=chr8-140300615-C-G"

VARIANT_RE = re.compile(
    "(chr)?(?P<chrom>[0-9XYMTt]{1,2})"
    "[-\s:]+"
    "(?P<pos>[0-9]{1,9})"
    "[-\s:]+"
    "(?P<ref>[ACGT]+)"
    "[-\s:>]+"
    "(?P<alt>[ACGT]+)"
)

USE_REDIS = True

if USE_REDIS:
    import redis
    REDIS = redis.Redis(host='localhost', port=6379, db=0)  # in-memory cache server which may or may not be running
else:
    REDIS = None


def error_response(error_message, source=None):
    response_json = {"error": str(error_message)}
    if source:
        response_json["source"] = source
    return Response(json.dumps(response_json), status=200, mimetype='application/json')


REVERSE_COMPLEMENT_MAP = dict(zip("ACGTN", "TGCAN"))


def reverse_complement(seq):
    return "".join([REVERSE_COMPLEMENT_MAP[n] for n in seq[::-1]])


def parse_variant(variant_str):
    match = VARIANT_RE.match(variant_str)
    if not match:
        raise ValueError(f"Unable to parse variant: {variant_str}")

    return match['chrom'], int(match['pos']), match['ref'], match['alt']


class VariantRecord:
    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alts = [alt]

    def __repr__(self):
        return f"{self.chrom}-{self.pos}-{self.ref}-{self.alts[0]}"


def get_splicing_scores_redis_key(tool_name, variant, genome_version, distance, mask):
    return f"{tool_name}__{variant}__hg{genome_version}__d{distance}__m{mask}"


def get_splicing_scores_from_redis(tool_name, variant, genome_version, distance, mask):
    if REDIS is None:
        return None

    key = get_splicing_scores_redis_key(tool_name, variant, genome_version, distance, mask)
    results = None
    try:
        results_string = REDIS.get(key)
        if results_string:
            results = json.loads(results_string)
            results["source"] += ":redis"
    except Exception as e:
        print(f"Redis error: {e}", flush=True)

    return results


def add_splicing_scores_to_redis(tool_name, variant, genome_version, distance, mask, results):
    if REDIS is None:
        return

    key = get_splicing_scores_redis_key(tool_name, variant, genome_version, distance, mask)
    try:
        results_string = json.dumps(results)
        REDIS.set(key, results_string)
    except Exception as e:
        print(f"Redis error: {e}", flush=True)


def exceeds_rate_limit(user_id, request_type):
    """Checks whether the given address has exceeded rate limits

    Args:
        user_id (str): unique user id
        request_type (str): type of rate limit - can be "spliceai:total", "spliceai:model", or "liftover:total"

    Return str: error message about exceeding the rate limit, or None if the rate limit was not exceeded
    """
    if REDIS is None:
        return False

    if request_type not in RATE_LIMIT_REQUESTS_PER_USER_PER_MINUTE:
        raise ValueError(f"Invalid 'request_type' arg value: {request_type}")

    epoch_time = time.time()  # seconds since 1970

    if epoch_time - int(REDIS.get("rate_limit_outlier_ips_update_time") or 0) > 120:  # time 2 minutes
        REDIS.set("rate_limit_outlier_ips_update_time", int(epoch_time))
        global RATE_LIMIT_OUTLIER_IPS
        RATE_LIMIT_OUTLIER_IPS = get_rate_limit_outlier_ips()

    if user_id in RATE_LIMIT_OUTLIER_IPS:
        print(f"Rate limiting outlier list IP: {user_id}")
        max_requests = 1
    else:
        max_requests_per_minute = RATE_LIMIT_REQUESTS_PER_USER_PER_MINUTE[request_type]
        max_requests = RATE_LIMIT_WINDOW_SIZE_IN_MINUTES * max_requests_per_minute

    try:
        # check number of requests from this user in the last (RATE_LIMIT_WINDOW_SIZE_IN_MINUTES * 60) minutes
        redis_key_prefix = f"request {user_id} {request_type}"
        keys = REDIS.keys(f"{redis_key_prefix}*")
        if len(keys) >= max_requests:
            redis_hit_limit_counter_key = f"request {user_id} rate limit counter"
            redis_hit_limit_counter = REDIS.get(redis_hit_limit_counter_key) or 0
            redis_hit_limit_counter = int(redis_hit_limit_counter) + 1
            REDIS.set(redis_hit_limit_counter_key, redis_hit_limit_counter)
            REDIS.expire(redis_hit_limit_counter_key, RATE_LIMIT_COUNTER_WINDOW_SIZE_IN_DAYS * 24 * 60 * 60)

            if redis_hit_limit_counter > 200:
                error_message = (
                    f"ERROR: You have exceeded the rate limit {redis_hit_limit_counter} times so far "
                    f"over the past few days. To prevent a single user from overwhelming the server and making it "
                    f"unavailable to other users, this tool allows no more than "
                    f"{RATE_LIMIT_REQUESTS_PER_USER_PER_MINUTE[request_type]} computed requests per "
                    f"minute per user. If you continue to exceed this limit, your IP address may be blocked."
                )
            else:
                error_message = (
                    f"ERROR: Rate limit reached. To prevent a user from overwhelming the server and making it "
                    f"unavailable to other users, this tool allows no more than "
                    f"{RATE_LIMIT_REQUESTS_PER_USER_PER_MINUTE[request_type]} computed requests per minute per user."
                )

            return error_message

        # record this request
        REDIS.set(f"{redis_key_prefix}: {epoch_time}", 1)
        REDIS.expire(f"{redis_key_prefix}: {epoch_time}", RATE_LIMIT_WINDOW_SIZE_IN_MINUTES * 60)
    except Exception as e:
        print(f"Redis error: {e}", flush=True)

    return None


def get_spliceai_scores(variant, genome_version, distance_param, mask_param):
    try:
        chrom, pos, ref, alt = parse_variant(variant)
    except ValueError as e:
        return {
            "variant": variant,
            "source": "spliceai",
            "error": f"ERROR: {e}",
        }

    # generate error message if variant falls outside annotated exons or introns
    OTHER_GENOME_VERSION = {"37": "38", "38": "37"}
    chrom_without_chr = chrom.replace("chr", "")
    if not ANNOTATION_INTERVAL_TREES[genome_version][chrom_without_chr].at(pos):
        other_genome_version = OTHER_GENOME_VERSION[genome_version]
        other_genome_overlapping_intervals = ANNOTATION_INTERVAL_TREES[other_genome_version][chrom_without_chr].at(pos)
        if other_genome_overlapping_intervals:
            other_genome_genes = " and ".join(sorted(set([str(i.data) for i in other_genome_overlapping_intervals])))
            return {
                "variant": variant,
                "source": "spliceai",
                "error": f"ERROR: In GRCh{genome_version}, {chrom}-{pos}-{ref}-{alt} falls outside all gencode exons and introns."
                         f"SpliceAI only works for variants within known exons or introns. However, in GRCh{other_genome_version}, "
                         f"{chrom}:{pos} falls within {other_genome_genes}, so perhaps GRCh{genome_version} is not the correct genome version?"
            }
        else:
            return {
                "variant": variant,
                "source": "spliceai",
                "error": f"ERROR: {chrom}-{pos}-{ref}-{alt} falls outside all Gencode exons and introns on "
                f"GRCh{genome_version}. SpliceAI only works for variants that are within known exons or introns.",
            }

            """
            NOTE: The reason SpliceAI currently works only for variants 
            within annotated exons or introns is that, although the SpliceAI neural net takes any 
            arbitrary nucleotide sequence as input, SpliceAI needs 1) the transcript strand 
            to determine whether to reverse-complement the reference genome sequence before passing it 
            to the neural net, and 2) transcript start and end positions to determine where to truncate 
            the reference genome sequence.
            """

    source = None
    scores = []

    # run the SpliceAI model to compute the scores
    if not scores:
        error_message = exceeds_rate_limit(request.remote_addr, request_type="spliceai:model")
        if error_message:
            return {
                "variant": variant,
                "source": "spliceai",
                "error": error_message,
            }

        record = VariantRecord(chrom, pos, ref, alt)
        try:
            scores = get_delta_scores(
                record,
                SPLICEAI_ANNOTATOR[genome_version],
                distance_param,
                mask_param)
            source = "spliceai:model"
        except Exception as e:
            print(f"ERROR while computing SpliceAI scores for {variant}: {e}")
            traceback.print_exc()
            return {
                "variant": variant,
                "source": "spliceai",
                "error": f"ERROR: {type(e)}: {e}",
            }

    if not scores:
        return {
            "variant": variant,
            "source": "spliceai",
            "error": f"ERROR: The SpliceAI model did not return any scores for {variant}. This is typically due to the "
                     f"variant falling outside of all Gencode exons and introns.",
        }

    #scores = [s[s.index("|")+1:] for s in scores]  # drop allele field

    # to reduce the response size, return all non-zero scores only for the canonial transcript (or the 1st transcript)
    all_non_zero_scores = None
    all_non_zero_scores_strand = None
    all_non_zero_scores_transcript_id = None
    all_non_zero_scores_transcript_priority = -1
    for i, transcript_scores in enumerate(scores):
        if "ALL_NON_ZERO_SCORES" not in transcript_scores:
            continue

        transcript_id_without_version = transcript_scores.get("NAME", "").split(".")[0]

        # get json annotations for this transcript
        transcript_annotations = TRANSCRIPT_ANNOTATIONS[genome_version].get(transcript_id_without_version)
        if transcript_annotations is None:
            raise ValueError(f"Missing annotations for {transcript_id_without_version} in {genome_version} annotations")

        # add the extra transcript annotations from the json file to the transcript scores dict
        transcript_scores.update(transcript_annotations)

        # decide whether to use ALL_NON_ZERO_SCORES from this transcript
        current_transcript_priority = TRANSCRIPT_PRIORITY_ORDER[transcript_annotations["t_priority"]]
        if current_transcript_priority > all_non_zero_scores_transcript_priority:
            all_non_zero_scores_transcript_priority = current_transcript_priority
            all_non_zero_scores = transcript_scores["ALL_NON_ZERO_SCORES"]
            all_non_zero_scores_strand = transcript_scores["t_strand"]
            all_non_zero_scores_transcript_id = transcript_scores["t_id"]

        for redundant_key in "ALLELE", "NAME", "STRAND", "ALL_NON_ZERO_SCORES":
            del transcript_scores[redundant_key]

    return {
        "variant": variant,
        "genomeVersion": genome_version,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "distance": distance_param,
        "scores": scores,
        "source": source,
        "allNonZeroScores": all_non_zero_scores,
        "allNonZeroScoresStrand": all_non_zero_scores_strand,
        "allNonZeroScoresTranscriptId": all_non_zero_scores_transcript_id,
    }


def get_pangolin_scores(variant, genome_version, distance_param, mask_param):
    if genome_version not in ("37", "38"):
        raise ValueError(f"Invalid genome_version: {mask_param}")

    if mask_param not in ("True", "False"):
        raise ValueError(f"Invalid mask_param: {mask_param}")

    try:
        chrom, pos, ref, alt = parse_variant(variant)
    except ValueError as e:
        print(f"ERROR while parsing variant {variant}: {e}")
        traceback.print_exc()

        return {
            "variant": variant,
            "source": "pangolin",
            "error": f"ERROR: {e}",
        }

    if len(ref) > 1 and len(alt) > 1:
        return {
            "variant": variant,
            "source": "pangolin",
            "error": f"ERROR: Pangolin does not currently support complex InDels like {chrom}-{pos}-{ref}-{alt}",
        }

    error_message = exceeds_rate_limit(request.remote_addr, request_type="pangolin:model")
    if error_message:
        return {
            "variant": variant,
            "source": "pangolin",
            "error": error_message,
        }

    class PangolinArgs:
        reference_file = HG19_FASTA_PATH if genome_version == "37" else HG38_FASTA_PATH
        distance = distance_param
        mask = mask_param
        score_cutoff = None
        score_exons = "False"

    if genome_version == "37":
        pangolin_gene_db = gffutils.FeatureDB(PANGOLIN_GRCH37_ANNOTATIONS)
    else:
        pangolin_gene_db = gffutils.FeatureDB(PANGOLIN_GRCH38_ANNOTATIONS)

    scores = process_variant_using_pangolin(
        0, chrom, int(pos), ref, alt, pangolin_gene_db, PANGOLIN_MODELS, PangolinArgs)

    if not scores:
        return {
            "variant": variant,
            "source": "pangolin",
            "error": f"ERROR: Pangolin was unable to compute scores for this variant",
        }

    # to reduce the response size, return all non-zero scores only for the canonial transcript (or the 1st transcript)
    all_non_zero_scores = None
    all_non_zero_scores_strand = None
    all_non_zero_scores_transcript_id = None
    max_delta_score_sum = 0
    for i, transcript_scores in enumerate(scores):
        if "ALL_NON_ZERO_SCORES" not in transcript_scores:
            continue

        transcript_id_without_version = transcript_scores.get("NAME", "").split(".")[0]

        # get json annotations for this transcript
        transcript_annotations = TRANSCRIPT_ANNOTATIONS[genome_version].get(transcript_id_without_version)
        if transcript_annotations is None:
            raise ValueError(f"Missing annotations for {transcript_id_without_version} in {genome_version} annotations")

        # add the extra transcript annotations from the json file to the transcript scores dict
        transcript_scores.update(transcript_annotations)

        # decide whether to use ALL_NON_ZERO_SCORES from this gene
        delta_score_sum = sum(abs(float(s.get("SG_ALT", 0)) - float(s.get("SG_REF", 0)))
                              for s in transcript_scores["ALL_NON_ZERO_SCORES"])
        delta_score_sum += sum(abs(float(s.get("SL_ALT", 0)) - float(s.get("SL_REF", 0)))
                               for s in transcript_scores["ALL_NON_ZERO_SCORES"])

        # return all_non_zero_scores for the transcript or gene with the highest delta score sum
        if delta_score_sum > max_delta_score_sum:
            all_non_zero_scores = transcript_scores["ALL_NON_ZERO_SCORES"]
            all_non_zero_scores_strand = transcript_scores["STRAND"]
            all_non_zero_scores_transcript_id = transcript_scores["NAME"]
            max_delta_score_sum = delta_score_sum

        for redundant_key in "NAME", "STRAND", "ALL_NON_ZERO_SCORES":
            del transcript_scores[redundant_key]

    return {
        "variant": variant,
        "genomeVersion": genome_version,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "distance": distance_param,
        "scores": scores,
        "source": "pangolin",
        "allNonZeroScores": all_non_zero_scores,
        "allNonZeroScoresStrand": all_non_zero_scores_strand,
        "allNonZeroScoresTranscriptId": all_non_zero_scores_transcript_id,
    }


@app.route("/spliceai/", methods=['POST', 'GET'])
def run_spliceai():
    return run_splice_prediction_tool(tool_name="spliceai")


@app.route("/pangolin/", methods=['POST', 'GET'])
def run_pangolin():
    return run_splice_prediction_tool(tool_name="pangolin")


def run_splice_prediction_tool(tool_name):
    """Handles API request for splice prediction

    Args:
        tool_name (str): "spliceai" or "pangolin"
    """
    if tool_name not in ("spliceai", "pangolin"):
        raise ValueError(f"Invalid tool_name: {tool_name}")

    start_time = datetime.now()
    logging_prefix = start_time.strftime("%m/%d/%Y %H:%M:%S") + f" t{os.getpid()}"

    # check params
    params = {}
    if request.values:
        params.update(request.values)

    if 'variant' not in params:
        params.update(request.get_json(force=True, silent=True) or {})

    error_message = exceeds_rate_limit(request.remote_addr, request_type=f"{tool_name}:total")
    if error_message:
        print(f"{logging_prefix}: {request.remote_addr}: response: {error_message}", flush=True)
        return error_response(error_message, source=tool_name)

    variant = params.get('variant', '')
    variant = variant.strip().strip("'").strip('"').strip(",")
    if not variant:
        return error_response(f'"variant" not specified. For example: {SPLICEAI_EXAMPLE}\n', source=tool_name)

    if not isinstance(variant, str):
        return error_response(f'"variant" value must be a string rather than a {type(variant)}.\n', source=tool_name)

    genome_version = params.get("hg")
    if not genome_version:
        return error_response(f'"hg" not specified. The URL must include an "hg" arg: hg=37 or hg=38. For example: {SPLICEAI_EXAMPLE}\n', source=tool_name)

    if genome_version not in ("37", "38"):
        return error_response(f'Invalid "hg" value: "{genome_version}". The value must be either "37" or "38". For example: {SPLICEAI_EXAMPLE}\n', source=tool_name)

    distance_param = params.get("distance", SPLICEAI_DEFAULT_DISTANCE)
    try:
        distance_param = int(distance_param)
    except Exception as e:
        return error_response(f'Invalid "distance": "{distance_param}". The value must be an integer.\n', source=tool_name)

    if distance_param > SPLICEAI_MAX_DISTANCE_LIMIT:
        return error_response(f'Invalid "distance": "{distance_param}". The value must be < {SPLICEAI_MAX_DISTANCE_LIMIT}.\n', source=tool_name)

    mask_param = params.get("mask", str(SPLICEAI_DEFAULT_MASK))
    if mask_param not in ("0", "1"):
        return error_response(f'Invalid "mask" value: "{mask_param}". The value must be either "0" or "1". For example: {SPLICEAI_EXAMPLE}\n', source=tool_name)

    if request.remote_addr not in DISABLE_LOGGING_FOR_IPS:
        print(f"{logging_prefix}: {request.remote_addr}: ======================", flush=True)
        print(f"{logging_prefix}: {request.remote_addr}: {variant} processing with hg={genome_version}, "
              f"distance={distance_param}, mask={mask_param}", flush=True)

    # check REDIS cache before processing the variant
    results = get_splicing_scores_from_redis(tool_name, variant, genome_version, distance_param, mask_param)
    if not results:
        try:
            if tool_name == "spliceai":
                results = get_spliceai_scores(variant, genome_version, distance_param, int(mask_param))
            elif tool_name == "pangolin":
                pangolin_mask_param = "True" if mask_param == "1" else "False"
                results = get_pangolin_scores(variant, genome_version, distance_param, pangolin_mask_param)
            else:
                raise ValueError(f"Invalid tool_name: {tool_name}")
        except Exception as e:
            traceback.print_exc()
            return error_response(f"ERROR: {e}", source=tool_name)

        if "error" not in results:
            add_splicing_scores_to_redis(tool_name, variant, genome_version, distance_param, mask_param, results)

    response_json = {}
    response_json.update(params)  # copy input params to output
    response_json.update(results)

    duration = str(datetime.now() - start_time)
    response_json['duration'] = duration

    if request.remote_addr not in DISABLE_LOGGING_FOR_IPS:
        print(f"{logging_prefix}: {request.remote_addr}: {variant} took {duration}", flush=True)

    return Response(json.dumps(response_json), status=200, mimetype='application/json')


LIFTOVER_EXAMPLE = f"/liftover/?hg=hg19-to-hg38&format=interval&chrom=chr8&start=140300615&end=140300620"

CHAIN_FILE_PATHS = {
    "hg19-to-hg38": "hg19ToHg38.over.chain.gz",
    "hg38-to-hg19": "hg38ToHg19.over.chain.gz",
    "hg38-to-t2t": "hg38ToHs1.over.chain.gz", # replaced hg38-chm13v2.over.chain.gz based on advice from Giulio Genovese
    "t2t-to-hg38": "hs1ToHg38.over.chain.gz", # replaced chm13v2-hg38.over.chain.gz based on advice from Giulio Genovese
}

LIFTOVER_REFERENCE_PATHS = {
    "hg19-to-hg38": (HG19_FASTA_PATH, HG38_FASTA_PATH),
    "hg38-to-hg19": (HG38_FASTA_PATH, HG19_FASTA_PATH),
    "hg38-to-t2t": (HG38_FASTA_PATH, T2T_FASTA_PATH),
    "t2t-to-hg38": (T2T_FASTA_PATH, HG38_FASTA_PATH),
}

def run_variant_liftover_tool(hg, chrom, pos, ref, alt, verbose=False):
    if hg not in CHAIN_FILE_PATHS or hg not in LIFTOVER_REFERENCE_PATHS:
        raise ValueError(f"Unexpected hg arg value: {hg}")
    chain_file_path = CHAIN_FILE_PATHS[hg]
    source_fasta_path, destination_fasta_path = LIFTOVER_REFERENCE_PATHS[hg]

    with tempfile.NamedTemporaryFile(suffix=".vcf", mode="wt", encoding="UTF-8") as input_file, \
            tempfile.NamedTemporaryFile(suffix=".vcf", mode="rt", encoding="UTF-8") as output_file:

        #  command syntax: liftOver oldFile map.chain newFile unMapped
        chrom = "chr" + chrom.replace("chr", "")
        input_file.write(f"""##fileformat=VCFv4.2
##contig=<ID={chrom},length=100000000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
{chrom}	{pos}	.	{ref}	{alt}	60	.""")
        input_file.flush()
        command = (
            f"cat {input_file.name} | "
            f"bcftools plugin liftover -- --src-fasta-ref {source_fasta_path} --fasta-ref {destination_fasta_path} --chain {chain_file_path} | "
            f"grep -v ^#  > {output_file.name}"
        )

        try:
            subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT, encoding="UTF-8")
            results = output_file.read()

            if verbose:
                print(f"{BCFTOOLS_LIFTOVER_TOOL} {hg} liftover on {chrom}:{pos} {ref}>{alt} returned: {results}", flush=True)

            # example: chr8	140300616	.	T	G	60	.	.

            result_fields = results.strip().split("\t")
            if len(result_fields) > 5:
                result_fields[1] = int(result_fields[1])

                return {
                    "hg": hg,
                    "chrom": chrom,
                    "start": int(pos) - 1,
                    "end": pos,
                    "output_chrom": result_fields[0],
                    "output_pos": result_fields[1],
                    "output_ref": result_fields[3],
                    "output_alt": result_fields[4],
                    "liftover_tool": BCFTOOLS_LIFTOVER_TOOL,
                    #"output_strand": "-" if "SWAP=-1" in results else "+",
                }

        except Exception as e:
            variant = f"{hg}  {chrom}:{pos} {ref}>{alt}"
            print(f"ERROR during liftover for {variant}: {e}")
            traceback.print_exc()
            raise ValueError(f"liftOver command failed for {variant}: {e}")

        # if bcftools liftover failed, fall back on running UCSC liftover
        result = run_UCSC_liftover_tool(hg, chrom, int(pos)-1, pos, verbose=False)
        result["output_ref"] = ref
        result["output_alt"] = alt
        #if result["output_strand"] == "-":
        #    result["output_ref"] = reverse_complement(result["output_ref"])
        #    result["output_alt"] = reverse_complement(result["output_alt"])
        return result


def run_UCSC_liftover_tool(hg, chrom, start, end, verbose=False):
    if hg not in CHAIN_FILE_PATHS:
        raise ValueError(f"Unexpected hg arg value: {hg}")
    chain_file_path = CHAIN_FILE_PATHS[hg]

    reason_liftover_failed = ""
    with tempfile.NamedTemporaryFile(suffix=".bed", mode="wt", encoding="UTF-8") as input_file, \
        tempfile.NamedTemporaryFile(suffix=".bed", mode="rt", encoding="UTF-8") as output_file, \
        tempfile.NamedTemporaryFile(suffix=".bed", mode="rt", encoding="UTF-8") as unmapped_output_file:

        #  command syntax: liftOver oldFile map.chain newFile unMapped
        chrom = "chr" + chrom.replace("chr", "")
        input_file.write("\t".join(map(str, [chrom, start, end, ".", "0", "+"])) + "\n")
        input_file.flush()
        command = f"liftOver {input_file.name} {chain_file_path} {output_file.name} {unmapped_output_file.name}"

        try:
            subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT, encoding="UTF-8")
            results = output_file.read()
            if verbose:
                print(f"{UCSC_LIFTOVER_TOOL} {hg} liftover on {chrom}:{start}-{end} returned: {results}", flush=True)

            result_fields = results.strip().split("\t")
            if len(result_fields) > 5:
                result_fields[1] = int(result_fields[1])
                result_fields[2] = int(result_fields[2])

                return {
                    "hg": hg,
                    "chrom": chrom,
                    "pos": int(start) + 1,
                    "start": start,
                    "end": end,
                    "output_chrom": result_fields[0],
                    "output_pos":   int(result_fields[1]) + 1,
                    "output_start": result_fields[1],
                    "output_end":    result_fields[2],
                    "output_strand": result_fields[5],
                    "liftover_tool": UCSC_LIFTOVER_TOOL,
                }
            else:
                reason_liftover_failed = unmapped_output_file.readline().replace("#", "").strip()

        except Exception as e:
            variant = f"{hg}  {chrom}:{start}-{end}"
            print(f"ERROR during liftover for {variant}: {e}")
            traceback.print_exc()
            raise ValueError(f"liftOver command failed for {variant}: {e}")

    if reason_liftover_failed:
        raise ValueError(f"{hg} liftover failed for {chrom}:{start}-{end} {reason_liftover_failed}")
    else:
        raise ValueError(f"{hg} liftover failed for {chrom}:{start}-{end} for unknown reasons")


def get_liftover_from_redis(key):
    if REDIS is None:
        return None

    results = None
    try:
        results_string = REDIS.get(key)
        if results_string:
            results = json.loads(results_string)
    except Exception as e:
        print(f"Redis error: {e}", flush=True)

    return results


def add_liftover_to_redis(key, result):
    if REDIS is None:
        return

    try:
        results_string = json.dumps(result)
        REDIS.set(key, results_string)
    except Exception as e:
        print(f"Redis error: {e}", flush=True)


@app.route("/liftover/", methods=['POST', 'GET'])
def run_liftover():
    logging_prefix = datetime.now().strftime("%m/%d/%Y %H:%M:%S") + f" t{os.getpid()}"

    # check params
    params = {}
    if request.values:
        params.update(request.values)

    if "format" not in params:
        params.update(request.get_json(force=True, silent=True) or {})

    error_message = exceeds_rate_limit(request.remote_addr, request_type="liftover:total")
    if error_message:
        print(f"{logging_prefix}: {request.remote_addr}: response: {error_message}", flush=True)
        return error_response(error_message)

    VALID_HG_VALUES = set(CHAIN_FILE_PATHS.keys())
    hg = params.get("hg")
    if not hg or hg not in VALID_HG_VALUES:
        return error_response(f'"hg" param error. It should be set to {" or ".join(VALID_HG_VALUES)}. For example: {LIFTOVER_EXAMPLE}\n')

    VALID_FORMAT_VALUES = ("interval", "variant", "position")
    format = params.get("format", "")
    if not format or format not in VALID_FORMAT_VALUES:
        return error_response(f'"format" param error. It should be set to {" or ".join(VALID_FORMAT_VALUES)}. For example: {LIFTOVER_EXAMPLE}\n')

    chrom = params.get("chrom")
    if not chrom:
        return error_response(f'"chrom" param not specified')

    if format == "interval":
        for key in "start", "end":
            if not params.get(key):
                return error_response(f'"{key}" param not specified')
        start = params.get("start")
        end = params.get("end")
        redis_key = f"{hg}:{chrom}:{start}:{end}"
        variant_log_string = f"{start}-{end}"

    elif format == "position":
        try:
            pos = int(params["pos"])
        except Exception as e:
            return error_response(f'"pos" param error: {e}')

        start = pos - 1
        end = pos
        redis_key = f"{hg}:{chrom}:{pos}"
        variant_log_string = f"{pos} "
    elif format == "variant":
        for key in "pos", "ref", "alt":
            if not params.get(key):
                return error_response(f'"{key}" param not specified')
        pos = params.get("pos")
        ref = params.get("ref")
        alt = params.get("alt")
        redis_key = f"{hg}:{chrom}:{pos}:{ref}:{alt}"
        variant_log_string = f"{pos} {ref}>{alt}"

    verbose = request.remote_addr not in DISABLE_LOGGING_FOR_IPS
    if verbose:
        print(f"{logging_prefix}: {request.remote_addr}: ======================", flush=True)
        print(f"{logging_prefix}: {request.remote_addr}: {hg} liftover {format}: {chrom}:{variant_log_string}", flush=True)

    # check REDIS cache before processing the variant
    result = get_liftover_from_redis(redis_key)
    if result and verbose:
        print(f"{hg} liftover on {variant_log_string} got results from cache: {result}", flush=True)

    if not result:
        try:
            if format == "variant":
                result = run_variant_liftover_tool(hg, chrom, pos, ref, alt, verbose=verbose)
            else:
                result = run_UCSC_liftover_tool(hg, chrom, start, end, verbose=verbose)
        except Exception as e:
            return error_response(str(e))
    
        add_liftover_to_redis(redis_key, result)

    result.update(params)

    return Response(json.dumps(result), mimetype='application/json')

# share static files from the annotations folder to support local installs
@app.route('/annotations/', strict_slashes=False, defaults={'path': ''})
@app.route('/annotations/<path:path>')
def send_annotations(path):
    if os.path.isfile(os.path.join("annotations", path)):
        return send_from_directory('annotations', path)

    # return an html table of available annotation files
    html = "<html><head><title>SpliceAI-lookup: Annotation Files</title></head>"
    html += "<body><table>"
    html += "<tr><th align=left>./annotation files</th><th align=left>last updated</th></tr>"
    for filename in os.listdir("annotations"):
        html += f"<tr><td><a href='/annotations/{filename}'>{filename}</a></td>"
        last_modified = datetime.fromtimestamp(os.path.getmtime(os.path.join('annotations', filename)))
        html += f"<td>{last_modified.strftime('%Y-%m-%d %H:%M:%S')}</td></tr>"
    html += "</table></body></html>"

    return Response(html, mimetype='text/html')



@app.route('/', strict_slashes=False, defaults={'path': ''})
@app.route('/<path:path>/')
def catch_all(path):
    if not path:
        path = "index.html"

    if path in {"index.html", "igv.min.js"}:
        with open(path, "rt") as f:
            html = f.read()
        return Response(html, mimetype='text/html')
    elif path == "favicon.ico":
        return send_from_directory('', 'favicon.ico')
    else:
        with open("README.md") as f:
            return markdown2.markdown(f.read())


print("Initialization completed.", flush=True)

if __name__ == "__main__":
    app.run(debug=DEBUG, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))
