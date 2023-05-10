from collections import defaultdict
from datetime import datetime
import json
import markdown2
import os
import pandas as pd
import pysam
import re
import redis
import socket
import subprocess
import sys
import tempfile
import time

# pangolin imports
from pkg_resources import resource_filename
from pangolin.model import torch, Pangolin, L, W, AR
from pangolin.pangolin import process_variant as process_variant_using_pangolin
import gffutils

# flask imports
from flask import Flask, request, Response
from flask_cors import CORS
from flask_talisman import Talisman
from intervaltree import IntervalTree, Interval
from spliceai.utils import Annotator, get_delta_scores

app = Flask(__name__)

CORS(app)

DEBUG = False
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

SPLICEAI_CACHE_FILES = {}
if socket.gethostname() == "spliceai-lookup":
    for filename in [
        "spliceai_scores.masked.indel.hg19.vcf.gz",
        "spliceai_scores.masked.indel.hg38.vcf.gz",
        "spliceai_scores.masked.snv.hg19.vcf.gz",
        "spliceai_scores.masked.snv.hg38.vcf.gz",
        "spliceai_scores.raw.indel.hg19.vcf.gz",
        "spliceai_scores.raw.indel.hg38.vcf.gz",
        "spliceai_scores.raw.snv.hg19.vcf.gz",
        "spliceai_scores.raw.snv.hg38.vcf.gz",
    ]:
        key = tuple(filename.replace("spliceai_scores.", "").replace(".vcf.gz", "").split("."))
        full_path = os.path.join("/mnt/disks/cache", filename)
        if os.path.isfile(full_path):
            SPLICEAI_CACHE_FILES[key] = pysam.TabixFile(full_path)
else:
    SPLICEAI_CACHE_FILES = {
        ("raw", "indel", "hg38"): pysam.TabixFile("./test_data/spliceai_scores.raw.indel.hg38_subset.vcf.gz"),
        ("raw", "snv", "hg38"): pysam.TabixFile("./test_data/spliceai_scores.raw.snv.hg38_subset.vcf.gz"),
        ("masked", "snv", "hg38"): pysam.TabixFile("./test_data/spliceai_scores.masked.snv.hg38_subset.vcf.gz"),
    }

GRCH37_ANNOTATIONS = "./annotations/gencode.v43lift37.annotation.txt.gz"
GRCH38_ANNOTATIONS = "./annotations/gencode.v43.annotation.txt.gz"
PANGOLIN_GRCH37_ANNOTATIONS = "./annotations/gencode.v43lift37.annotation.db"
PANGOLIN_GRCH38_ANNOTATIONS = "./annotations/gencode.v43.annotation.db"

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

for genome_version, annotation_path in ("37", GRCH37_ANNOTATIONS), ("38", GRCH38_ANNOTATIONS):
    print(f"Loading {annotation_path}", flush=True)
    df = pd.read_table(annotation_path, dtype={"TX_START": int, "TX_END": int})
    for _, row in df.iterrows():
        chrom = row["CHROM"].replace("chr", "")
        ANNOTATION_INTERVAL_TREES[genome_version][chrom].add(Interval(row["TX_START"], row["TX_END"] + 0.1, row["#NAME"]))

SPLICEAI_ANNOTATOR = {
    "37": Annotator(HG19_FASTA_PATH, GRCH37_ANNOTATIONS),
    "38": Annotator(HG38_FASTA_PATH, GRCH38_ANNOTATIONS),
}

SPLICEAI_MAX_DISTANCE_LIMIT = 10000
SPLICEAI_DEFAULT_DISTANCE = 50  # maximum distance between the variant and gained/lost splice site, defaults to 50
SPLICEAI_DEFAULT_MASK = 0       # mask scores representing annotated acceptor/donor gain and unannotated acceptor/donor loss, defaults to 0
USE_PRECOMPUTED_SCORES = 1      # whether to use precomputed scores by default

SPLICEAI_SCORE_FIELDS = "ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL".split("|")

SPLICEAI_EXAMPLE = f"/spliceai/?hg=38&distance=50&mask=0&precomputed=1&variant=chr8-140300615-C-G"

VARIANT_RE = re.compile(
    "(chr)?(?P<chrom>[0-9XYMTt]{1,2})"
    "[-\s:]+"
    "(?P<pos>[0-9]{1,9})"
    "[-\s:]+"
    "(?P<ref>[ACGT]+)"
    "[-\s:>]+"
    "(?P<alt>[ACGT]+)"
)

REDIS = redis.Redis(host='localhost', port=6379, db=0)  # in-memory cache server which may or may not be running


def error_response(error_message):
    return Response(json.dumps({"error": str(error_message)}), status=400, mimetype='application/json')


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


def get_splicing_scores_redis_key(tool_name, variant, genome_version, distance, mask, use_precomputed_scores):
    return f"{tool_name}__{variant}__hg{genome_version}__d{distance}__m{mask}__pre{use_precomputed_scores}"


def get_splicing_scores_from_redis(tool_name, variant, genome_version, distance, mask, use_precomputed_scores):
    key = get_splicing_scores_redis_key(tool_name, variant, genome_version, distance, mask, use_precomputed_scores)
    results = None
    try:
        results_string = REDIS.get(key)
        if results_string:
            results = json.loads(results_string)
            results["source"] += ":redis"
    except Exception as e:
        print(f"Redis error: {e}", flush=True)

    return results


def add_splicing_scores_to_redis(tool_name, variant, genome_version, distance, mask, use_precomputed_scores, results):
    key = get_splicing_scores_redis_key(tool_name, variant, genome_version, distance, mask, use_precomputed_scores)
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


def get_spliceai_scores(variant, genome_version, distance_param, mask_param, use_precomputed_scores):
    try:
        chrom, pos, ref, alt = parse_variant(variant)
    except ValueError as e:
        return {
            "variant": variant,
            "error": f"ERROR: {e}",
        }

    if len(ref) > 1 and len(alt) > 1:
        return {
            "variant": variant,
            "error": f"ERROR: SpliceAI does not currently support complex InDels like {chrom}-{pos}-{ref}-{alt}",
        }

    # generate error message if variant falls outside annotated exons or introns
    OTHER_GENOME_VERSION = {"37": "38", "38": "37"}
    chrom_without_chr = chrom.replace("chr", "")
    if not ANNOTATION_INTERVAL_TREES[genome_version][chrom_without_chr].at(pos):
        other_genome_version = OTHER_GENOME_VERSION[genome_version]
        other_genome_overlapping_intervals = ANNOTATION_INTERVAL_TREES[other_genome_version][chrom_without_chr].at(pos)
        if other_genome_overlapping_intervals:
            other_genome_genes = " and ".join(sorted(set([str(i.data).split("---")[0] for i in other_genome_overlapping_intervals])))
            return {
                "variant": variant,
                "error": f"ERROR: In GRCh{genome_version}, {chrom}-{pos}-{ref}-{alt} falls outside all gencode exons and introns."
                         f"SpliceAI only works for variants within known exons or introns. However, in GRCh{other_genome_version}, "
                         f"{chrom}:{pos} falls within {other_genome_genes}, so perhaps GRCh{genome_version} is not the correct genome version?"
            }
        else:
            return {
                "variant": variant,
                "error": f"ERROR: {chrom}-{pos}-{ref}-{alt} falls outside all Gencode exons and introns on "
                f"GRCh{genome_version}. SpliceAI only works for variants that are within known exons or introns.",
            }

            """
            NOTE: The reason SpliceAI currently works only for variants "
                         f"within annotated exons or introns is that, although the SpliceAI neural net takes any "
                         f"arbitrary nucleotide sequence as input, SpliceAI needs 1) the transcript strand "
                         f"to determine whether to reverse-complement the reference genome sequence before passing it "
                         f"to the neural net, and 2) transcript start and end positions to determine where to truncate "
                         f"the reference genome sequence.
            """

    source = None
    scores = []
    if (len(ref) <= 5 or len(alt) <= 2) and str(distance_param) == str(SPLICEAI_DEFAULT_DISTANCE) and str(use_precomputed_scores) == "1":
        # examples: ("masked", "snv", "hg19")  ("raw", "indel", "hg38")
        key = (
            "masked" if str(mask_param) == "1" else ("raw" if str(mask_param) == "0" else None),
            "snv" if len(ref) == 1 and len(alt) == 1 else "indel",
            "hg19" if genome_version == "37" else ("hg38" if genome_version == "38" else None),
        )
        try:
            results = SPLICEAI_CACHE_FILES[key].fetch(chrom, pos-1, pos+1)
            for line in results:
                # ['1', '739023', '.', 'C', 'CT', '.', '.', 'SpliceAI=CT|AL669831.1|0.00|0.00|0.00|0.00|-1|-37|-48|-37']
                fields = line.split("\t")
                if fields[0] == chrom and int(fields[1]) == pos and fields[3] == ref and fields[4] == alt:
                    scores.append(fields[7])
            if scores:
                source = "splice-ai:lookup"
                #print(f"Fetched: ", scores, flush=True)

        except Exception as e:
            print(f"ERROR: couldn't retrieve scores using tabix: {type(e)}: {e}", flush=True)

    if not scores:
        error_message = exceeds_rate_limit(request.remote_addr, request_type="spliceai:model")
        if error_message:
            return {
                "variant": variant,
                "error": error_message,
            }

        record = VariantRecord(chrom, pos, ref, alt)
        try:
            scores = get_delta_scores(
                record,
                SPLICEAI_ANNOTATOR[genome_version],
                distance_param,
                mask_param)
            source = "splice-ai:model"
            #print(f"Computed: ", scores, flush=True)
        except Exception as e:
            return {
                "variant": variant,
                "error": f"ERROR: {type(e)}: {e}",
            }

    if not scores:
        return {
            "variant": variant,
            "error": f"ERROR: The SpliceAI model did not return any scores for {variant}. This is typically due to the "
                     f"variant falling outside of all Gencode exons and introns.",
        }

    scores = [s[s.index("|")+1:] for s in scores]  # drop allele field

    return {
        "variant": variant,
        "genome_version": genome_version,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "scores": scores,
        "source": source,
    }


def get_pangolin_scores(variant, genome_version, distance_param, mask_param, use_precomputed_scores):
    if genome_version not in ("37", "38"):
        raise ValueError(f"Invalid genome_version: {mask_param}")

    if mask_param not in ("True", "False"):
        raise ValueError(f"Invalid mask_param: {mask_param}")

    try:
        chrom, pos, ref, alt = parse_variant(variant)
    except ValueError as e:
        return {
            "variant": variant,
            "error": f"ERROR: {e}",
        }

    if len(ref) > 1 and len(alt) > 1:
        return {
            "variant": variant,
            "error": f"ERROR: Pangolin does not currently support complex InDels like {chrom}-{pos}-{ref}-{alt}",
        }

    error_message = exceeds_rate_limit(request.remote_addr, request_type="pangolin:model")
    if error_message:
        return {
            "variant": variant,
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

    if scores == -1:
        return {
            "variant": variant,
            "error": f"ERROR: Unable to compute Pangolin scores",
        }

    parsed_scores = []
    for i, scores_for_gene in enumerate(scores.split("ENSG")):
        if i == 0:
            continue

        scores_for_gene = scores_for_gene.split("|")

        gene_id = "ENSG" + scores_for_gene[0]
        splice_gain_pos, splice_gain_score = scores_for_gene[1].split(":")
        splice_loss_pos, splice_loss_score = scores_for_gene[2].split(":")
        warnings = scores_for_gene[3].replace("Warnings:", "").strip()
        parsed_scores.append(
            "|".join([gene_id, splice_gain_score, splice_loss_score, splice_gain_pos, splice_loss_pos])
        )

        if warnings:
            print("Pangolin Warning:", warnings)

    return {
        "variant": variant,
        "genome_version": genome_version,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "scores": parsed_scores,
        "source": "pangolin",
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
        return error_response(error_message)

    variant = params.get('variant', '')
    variant = variant.strip().strip("'").strip('"').strip(",")
    if not variant:
        return error_response(f'"variant" not specified. For example: {SPLICEAI_EXAMPLE}\n')

    if not isinstance(variant, str):
        return error_response(f'"variant" value must be a string rather than a {type(variant)}.\n')

    genome_version = params.get("hg")
    if not genome_version:
        return error_response(f'"hg" not specified. The URL must include an "hg" arg: hg=37 or hg=38. For example: {SPLICEAI_EXAMPLE}\n')

    if genome_version not in ("37", "38"):
        return error_response(f'Invalid "hg" value: "{genome_version}". The value must be either "37" or "38". For example: {SPLICEAI_EXAMPLE}\n')

    distance_param = params.get("distance", SPLICEAI_DEFAULT_DISTANCE)
    try:
        distance_param = int(distance_param)
    except Exception as e:
        return error_response(f'Invalid "distance": "{distance_param}". The value must be an integer.\n')

    if distance_param > SPLICEAI_MAX_DISTANCE_LIMIT:
        return error_response(f'Invalid "distance": "{distance_param}". The value must be < {SPLICEAI_MAX_DISTANCE_LIMIT}.\n')

    mask_param = params.get("mask", str(SPLICEAI_DEFAULT_MASK))
    if mask_param not in ("0", "1"):
        return error_response(f'Invalid "mask" value: "{mask_param}". The value must be either "0" or "1". For example: {SPLICEAI_EXAMPLE}\n')

    use_precomputed_scores = params.get("precomputed", str(USE_PRECOMPUTED_SCORES))
    if use_precomputed_scores not in ("0", "1"):
        return error_response(f'Invalid "precomputed" value: "{use_precomputed_scores}". The value must be either "0" or "1". For example: {SPLICEAI_EXAMPLE}\n')

    use_precomputed_scores = int(use_precomputed_scores)

    if request.remote_addr not in DISABLE_LOGGING_FOR_IPS:
        print(f"{logging_prefix}: {request.remote_addr}: ======================", flush=True)
        print(f"{logging_prefix}: {request.remote_addr}: {variant} processing with hg={genome_version}, "
              f"distance={distance_param}, mask={mask_param}, precomputed={use_precomputed_scores}", flush=True)

    # check REDIS cache before processing the variant
    results = get_splicing_scores_from_redis(tool_name, variant, genome_version, distance_param, mask_param, use_precomputed_scores)
    if not results:
        if tool_name == "spliceai":
            results = get_spliceai_scores(variant, genome_version, distance_param, int(mask_param), use_precomputed_scores)
        elif tool_name == "pangolin":
            pangolin_mask_param = "True" if mask_param == "1" else "False"
            results = get_pangolin_scores(variant, genome_version, distance_param, pangolin_mask_param, use_precomputed_scores)
        else:
            raise ValueError(f"Invalid tool_name: {tool_name}")

        if "error" not in results:
            add_splicing_scores_to_redis(tool_name, variant, genome_version, distance_param, mask_param, use_precomputed_scores, results)

    status = 400 if results.get("error") else 200

    response_json = {}
    response_json.update(params)  # copy input params to output
    response_json.update(results)

    duration = str(datetime.now() - start_time)
    response_json['duration'] = duration

    if request.remote_addr not in DISABLE_LOGGING_FOR_IPS:
        print(f"{logging_prefix}: {request.remote_addr}: {variant} response: {response_json}", flush=True)
        print(f"{logging_prefix}: {request.remote_addr}: {variant} took {duration}", flush=True)

    return Response(json.dumps(response_json), status=status, mimetype='application/json')


LIFTOVER_EXAMPLE = f"/liftover/?hg=hg19-to-hg38&format=interval&chrom=chr8&start=140300615&end=140300620"

CHAIN_FILE_PATHS = {
    "hg19-to-hg38": "hg19ToHg38.over.chain.gz",
    "hg38-to-hg19": "hg38ToHg19.over.chain.gz",
    "hg38-to-t2t": "hg38-chm13v2.over.chain.gz",
    "t2t-to-hg38": "chm13v2-hg38.over.chain.gz",
}


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
                print(f"{hg} liftover on {chrom}:{start}-{end} returned: {results}", flush=True)

            result_fields = results.strip().split("\t")
            if len(result_fields) > 5:
                result_fields[1] = int(result_fields[1])
                result_fields[2] = int(result_fields[2])

                return {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "output_chrom": result_fields[0],
                    "output_start": result_fields[1],
                    "output_end": result_fields[2],
                    "output_strand": result_fields[5],
                }
            else:
                reason_liftover_failed = unmapped_output_file.readline().replace("#", "").strip()
        except Exception as e:
            raise ValueError(f"{hg} liftOver command failed for {chrom}:{start}-{end}: {e}")

    if reason_liftover_failed:
        raise ValueError(f"{hg} liftover failed for {chrom}:{start}-{end} {reason_liftover_failed}")
    else:
        raise ValueError(f"{hg} liftover failed for {chrom}:{start}-{end} for unknown reasons")


def get_liftover_redis_key(genome_version, chrom, start, end):
    return f"liftover_hg{genome_version}__{chrom}_{start}_{end}"


def get_liftover_from_redis(hg, chrom, start, end):
    key = get_liftover_redis_key(hg, chrom, start, end)
    results = None
    try:
        results_string = REDIS.get(key)
        if results_string:
            results = json.loads(results_string)
    except Exception as e:
        print(f"Redis error: {e}", flush=True)

    return results


def add_liftover_to_redis(hg, chrom, start, end, result):
    key = get_liftover_redis_key(hg, chrom, start, end)
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
        start = params.get("start")
        end = params.get("end")
        if not start:
            return error_response(f'"start" param not specified')
        if not end:
            return error_response(f'"end" param not specified')
        variant_log_string = f"{start}-{end}"

    elif format == "position" or format == "variant":
        pos = params.get("pos")
        if not pos:
            return error_response(f'"pos" param not specified')

        pos = int(pos)
        start = pos - 1
        end = pos
        variant_log_string = f"{pos} "
        if params.get('ref') and params.get('alt'):
            variant_log_string += f"{params.get('ref')}>{params.get('alt')}"

    verbose = request.remote_addr not in DISABLE_LOGGING_FOR_IPS
    if verbose:
        print(f"{logging_prefix}: {request.remote_addr}: ======================", flush=True)
        print(f"{logging_prefix}: {request.remote_addr}: {hg} liftover {format}: {chrom}:{variant_log_string}", flush=True)

    # check REDIS cache before processing the variant
    result = get_liftover_from_redis(hg, chrom, start, end)
    if result and verbose:
        print(f"{hg} liftover on {chrom}:{start}-{end} got results from cache: {result}", flush=True)

    if not result:
        try:
            result = run_UCSC_liftover_tool(hg, chrom, start, end, verbose=verbose)
        except Exception as e:
            return error_response(str(e))
    
        add_liftover_to_redis(hg, chrom, start, end, result)

    result.update(params)
    if format == "position" or format == "variant":
        result["pos"] = pos
        result["output_pos"] = result["output_end"]

    if format == "variant":
        result["output_ref"] = result["ref"]
        result["output_alt"] = result["alt"]
        if result["output_strand"] == "-":
            result["output_ref"] = reverse_complement(result["output_ref"])
            result["output_alt"] = reverse_complement(result["output_alt"])

    return Response(json.dumps(result), mimetype='application/json')


@app.route('/', defaults={'path': ''})
@app.route('/<path:path>/')
def catch_all(path):
    with open("README.md") as f:
        return markdown2.markdown(f.read())


print("Initialization completed.", flush=True)

if __name__ == "__main__":
    app.run(debug=DEBUG, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))
