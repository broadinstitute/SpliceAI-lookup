from collections import defaultdict
from datetime import datetime
import json
import os
import pandas as pd
from pprint import pprint
import pyfastx
import psycopg2
import re
import traceback
import time

# flask imports
from flask import Flask, request, Response, send_from_directory
from flask_cors import CORS
from flask_talisman import Talisman

app = Flask(__name__)

CORS(app)


DEBUG = True # if socket.gethostname() == "spliceai-lookup" else True
if not DEBUG:
    Talisman(app)


DEFAULT_DISTANCE = 500  # maximum distance between the variant and gained/lost splice site, defaults to 500
MAX_DISTANCE_LIMIT = 10000
DEFAULT_MASK = 0        # mask scores representing annotated acceptor/donor gain and unannotated acceptor/donor loss, defaults to 0

SPLICEAI_EXAMPLE_URL = f"/spliceai/?hg=38&distance=500&mask=0&variant=chr8-140300615-C-G"
PANGOLIN_EXAMPLE_URL = f"/pangolin/?hg=38&distance=500&mask=0&variant=chr8-140300615-C-G"


VARIANT_RE = re.compile(
    "(chr)?(?P<chrom>[0-9XYMTt]{1,2})"
    "[-\s:]+"
    "(?P<pos>[0-9]{1,9})"
    "[-\s:]+"
    "(?P<ref>[ACGT]+)"
    "[-\s:>]+"
    "(?P<alt>[ACGT]+)"
)

FASTA_PATH = {
    "37": "/hg19.fa.gz",
    "38": "/hg38.fa.gz",
}

PYFASTX_REF = {}

GENCODE_VERSION = "v46"

SHARED_TRANSCRIPT_ANNOTATIONS = {}
SHARED_TRANSCRIPT_ANNOTATION_PATHS = {
    "37": f"/gencode.{GENCODE_VERSION}lift37.basic.annotation.transcript_annotations.json",
    "38": f"/gencode.{GENCODE_VERSION}.basic.annotation.transcript_annotations.json",
}

TRANSCRIPT_PRIORITY_ORDER = {
    "MS": 3,  # MANE select transcript
    "MP": 2,  # MANE plus clinical transcript
    "C": 1,   # canonical transcript
    "N": 0
}

TOOL = os.environ.get("TOOL")
GENOME_VERSION = os.environ.get("GENOME_VERSION")
if GENOME_VERSION not in ("37", "38"):
    raise ValueError(f'Environment variable "GENOME_VERSION" should be set to either "37" or "38" instead of: "{os.environ.get("GENOME_VERSION")}"')

if TOOL == "spliceai":
    from spliceai.utils import Annotator, get_delta_scores

    class VariantRecord:
        def __init__(self, chrom, pos, ref, alt):
            self.chrom = chrom
            self.pos = pos
            self.ref = ref
            self.alts = [alt]

        def __repr__(self):
            return f"{self.chrom}-{self.pos}-{self.ref}-{self.alts[0]}"

    SPLICEAI_ANNOTATOR = {}
    SPLICEAI_ANNOTATION_PATHS = {
        "37": f"/gencode.{GENCODE_VERSION}lift37.basic.annotation.txt.gz",
        "38": f"/gencode.{GENCODE_VERSION}.basic.annotation.txt.gz",
    }

elif TOOL == "pangolin":
    from pkg_resources import resource_filename
    from pangolin.pangolin import process_variant as process_variant_using_pangolin
    from pangolin.model import torch, Pangolin, L, W, AR
    import gffutils

    PANGOLIN_MODELS = []

    PANGOLIN_ANNOTATION_PATHS = {
        "37": f"/gencode.{GENCODE_VERSION}lift37.basic.annotation.without_chr_prefix.db",
        "38": f"/gencode.{GENCODE_VERSION}.basic.annotation.db",
    }
else:
    raise ValueError(f'Environment variable "TOOL" should be set to either "spliceai" or "pangolin" instead of: "{os.environ.get("TOOL")}"')


def init_reference(genome_version):
    if genome_version not in PYFASTX_REF:
        PYFASTX_REF[genome_version] = pyfastx.Fasta(FASTA_PATH[genome_version])


def init_spliceai(genome_version):
    
    if genome_version not in SPLICEAI_ANNOTATOR:
        SPLICEAI_ANNOTATOR[genome_version] = Annotator(FASTA_PATH[genome_version], SPLICEAI_ANNOTATION_PATHS[genome_version])

def init_pangolin(genome_version):

    if not PANGOLIN_MODELS:
        torch.set_num_threads(os.cpu_count()*2)

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

def init_transcript_annotations(genome_version):
    if genome_version in SHARED_TRANSCRIPT_ANNOTATIONS:
        return

    # init shared transcript annotations
    with open(SHARED_TRANSCRIPT_ANNOTATION_PATHS[genome_version], "rt") as ta_f:
        SHARED_TRANSCRIPT_ANNOTATIONS[genome_version] = json.load(ta_f)


def error_response(error_message, source=None):
    response_json = {"error": str(error_message)}
    if source:
        response_json["source"] = source
    return Response(json.dumps(response_json), status=200, mimetype='application/json')


def parse_variant(variant_str):
    match = VARIANT_RE.match(variant_str)
    if not match:
        raise ValueError(f"Unable to parse variant: {variant_str}")

    return match['chrom'], int(match['pos']), match['ref'], match['alt']


DATABASE_CONNECTION = None
def init_database_connection():
    global DATABASE_CONNECTION
    if DATABASE_CONNECTION is not None:
        try:
            DATABASE_CONNECTION.isolation_level
            return  # connection is already initialized
        except psycopg2.OperationalError as oe:
            print("WARNING: Connection has closed. Reopenning...")

    print("Environment:")
    pprint(dict(os.environ))
    try:
        DATABASE_CONNECTION = psycopg2.connect(
            host="/cloudsql/spliceai-lookup-412920:us-central1:spliceai-lookup-db",
            database="spliceai-lookup-db",
            user="postgres",
            password=os.environ.get("DB_PASSWORD"))
        DATABASE_CONNECTION.autocommit = True
    except Exception as e:
        print(f"ERROR: Unable to connect to the database: {e}")
        traceback.print_exc()
        return

    def does_table_exist(table_name):
        results = run_sql(f"SELECT EXISTS (SELECT 1 AS result FROM pg_tables WHERE tablename=%s)", (table_name,))
        does_table_already_exist = results[0][0]
        return does_table_already_exist

    if not does_table_exist("cache"):
        print("Creating cache table")
        run_sql("""CREATE TABLE cache (key TEXT UNIQUE, value TEXT, counter INT, accessed TIMESTAMP DEFAULT now())""")
        run_sql("""CREATE INDEX cache_index ON cache (key)""")

    if not does_table_exist("log"):
        print("Creating event_log table")
        run_sql("""CREATE TABLE log (event_name TEXT, ip TEXT, logtime TIMESTAMP DEFAULT now(), duration REAL, variant TEXT, genome VARCHAR(10), distance INT, mask INT4, details TEXT)""")
        run_sql("""CREATE INDEX log_index1 ON log (event_name)""")
        run_sql("""CREATE INDEX log_index2 ON log (ip)""")


def run_sql(sql_query, *args):
    with DATABASE_CONNECTION:
        c = DATABASE_CONNECTION.cursor()
        c.execute(sql_query, *args)
        try:
            results = c.fetchall()
        except:
            results = []
        c.close()
    return results


def exceeds_rate_limit(user_ip):
    """Place holder for implementing rate limits based on user ip address"""
    if DATABASE_CONNECTION is None:
        return False

    """
    SELECT * FROM log WHERE event_name like '%computed' AND duration > 2 AND ip='61.68.128.103' AND logtime >= NOW() - INTERVAL '5 minutes' ;
    SELECT ip, count(*) FROM log WHERE event_name like '%computed' AND duration > 2 AND logtime >= NOW() - INTERVAL '20 minutes' GROUP BY ip ORDER BY count DESC;
    """
    return False

    #rows = run_sql("SELECT COUNT(*) FROM log WHERE event_name like '%computed' AND duration > 2 AND ip=%s AND logtime >= NOW() - INTERVAL '5 minutes'", user_ip)
    #request_count = int(rows[0][0])
    #if request_count > 50: return True

    # check rate limit
    #try:
    #    rows = run_sql(f"SELECT count(*) FROM log WHERE ip=%s AND logtime > now() - interval '1 second' AND event_name=%s", (ip, request_type))
    #    if rows and rows[0][0] > 4:
    #        return f"Rate limit exceeded for {request_type} requests"
    #except Exception as e:
    #    print(f"Rate limit error: {e}", flush=True)


def get_splicing_scores_cache_key(tool_name, variant, genome_version, distance, mask):
    return f"{tool_name}__{variant}__hg{genome_version}__d{distance}__m{mask}"


def get_splicing_scores_from_cache(tool_name, variant, genome_version, distance, mask):
    if DATABASE_CONNECTION is None:
        return None

    key = get_splicing_scores_cache_key(tool_name, variant, genome_version, distance, mask)
    results = None
    try:
        rows = run_sql(f"SELECT value FROM cache WHERE key=%s", (key,))
        if rows:
            results = json.loads(rows[0][0])
            results["source"] += ":cache"
    except Exception as e:
        print(f"Cache error: {e}", flush=True)

    return results


def add_splicing_scores_to_cache(tool_name, variant, genome_version, distance, mask, results):
    if DATABASE_CONNECTION is None:
        return

    key = get_splicing_scores_cache_key(tool_name, variant, genome_version, distance, mask)
    try:
        results_string = json.dumps(results)

        run_sql(r"""INSERT INTO cache (key, value, counter, accessed) VALUES (%s, %s, 1, now()) """ +
                r"""ON CONFLICT (key) DO """ +
                r"""UPDATE SET key=%s, value=%s, counter=cache.counter+1, accessed=now()""", (key, results_string, key, results_string))
    except Exception as e:
        print(f"Cache error: {e}", flush=True)


def check_reference_allele(genome_version, chrom, pos, ref, alt):
    # check that variant matches the reference allele
    if genome_version not in PYFASTX_REF:
        raise ValueError(f"Invalid genome_version: {genome_version}")

    chrom = chrom.replace("chr", "")
    if genome_version == "37":
        if chrom.upper() in ("M", "MT"):
            chrom = "MT"
    else:
        if chrom.upper() in ("M", "MT"):
            chrom = "M"
        chrom = "chr" + chrom


    variant = f"{chrom}-{pos}-{ref}-{alt}"
    try:
        ref_sequence = PYFASTX_REF[genome_version][chrom][pos-1:pos+len(ref)-1].seq
        if ref_sequence.upper() != ref.upper():
            return {
                "variant": variant,
                "source": "spliceai",
                "error": f"Unexpected reference allele in {chrom}-{pos}-{ref}-{alt}. The reference allele should be: {ref_sequence}",
            }
    except Exception as e:
        print(f"ERROR while checking the reference allele for {variant}: {e}")

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

    error_resonse = check_reference_allele(genome_version, chrom, pos, ref, alt)
    if error_resonse:
        return error_resonse

    # generate error message if variant falls outside annotated exons or introns
    OTHER_GENOME_VERSION = {"37": "38", "38": "37"}
    chrom_without_chr = chrom.replace("chr", "")
    source = None
    scores = []

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
            "error": f"ERROR: The SpliceAI model did not return any scores for {variant}",
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
        transcript_annotations = SHARED_TRANSCRIPT_ANNOTATIONS[genome_version].get(transcript_id_without_version)
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

    error_resonse = check_reference_allele(genome_version, chrom, pos, ref, alt)
    if error_resonse:
        return error_resonse

    if len(ref) > 1 and len(alt) > 1:
        return {
            "variant": variant,
            "source": "pangolin",
            "error": f"ERROR: Pangolin does not currently support complex InDels like {chrom}-{pos}-{ref}-{alt}",
        }

    class PangolinArgs:
        reference_file = FASTA_PATH[genome_version]
        distance = distance_param
        mask = mask_param
        score_cutoff = None
        score_exons = "False"

    features_db = gffutils.FeatureDB(PANGOLIN_ANNOTATION_PATHS[GENOME_VERSION])
    scores = process_variant_using_pangolin(
        0, chrom, int(pos), ref, alt, features_db, PANGOLIN_MODELS, PangolinArgs)

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
        transcript_annotations = SHARED_TRANSCRIPT_ANNOTATIONS[genome_version].get(transcript_id_without_version)
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

    if tool_name != TOOL:
        return error_response(f"ERROR: This server is configured to run {TOOL} rather than {tool_name}.\n", source=tool_name)

    start_time = datetime.now()
    logging_prefix = start_time.strftime("%m/%d/%Y %H:%M:%S") + f" t{os.getpid()}"
    example_url = SPLICEAI_EXAMPLE_URL if tool_name == "spliceai" else PANGOLIN_EXAMPLE_URL

    # check params
    params = {}
    if request.values:
        params.update(request.values)

    if 'variant' not in params:
        params.update(request.get_json(force=True, silent=True) or {})

    variant = params.get('variant', '')
    variant = variant.strip().strip("'").strip('"').strip(",")
    if not variant:
        return error_response(f'"variant" not specified.\n', source=tool_name)

    if not isinstance(variant, str):
        return error_response(f'"variant" value must be a string rather than a {type(variant)}.\n', source=tool_name)

    genome_version = params.get("hg")
    if not genome_version:
        return error_response(f'"hg" not specified. The URL must include an "hg" arg: hg=37 or hg=38. For example: {example_url}\n', source=tool_name)

    if genome_version not in ("37", "38"):
        return error_response(f'Invalid "hg" value: "{genome_version}". The value must be either "37" or "38". For example: {example_url}\n', source=tool_name)

    distance_param = params.get("distance", DEFAULT_DISTANCE)
    try:
        distance_param = int(distance_param)
    except Exception as e:
        return error_response(f'Invalid "distance": "{distance_param}". The value must be an integer.\n', source=tool_name)

    if distance_param > MAX_DISTANCE_LIMIT:
        return error_response(f'Invalid "distance": "{distance_param}". The value must be < {MAX_DISTANCE_LIMIT}.\n', source=tool_name)

    mask_param = params.get("mask", str(DEFAULT_MASK))
    if mask_param not in ("0", "1"):
        return error_response(f'Invalid "mask" value: "{mask_param}". The value must be either "0" or "1". For example: {example_url}\n', source=tool_name)

    print(f"{logging_prefix}: ======================", flush=True)
    print(f"{logging_prefix}: {variant} processing with hg={genome_version}, "
          f"distance={distance_param}, mask={mask_param}", flush=True)

    
    if tool_name == "spliceai":
        init_spliceai(genome_version)
        example_url = SPLICEAI_EXAMPLE_URL
    elif tool_name == "pangolin":
        init_pangolin(genome_version)
        example_url = PANGOLIN_EXAMPLE_URL
    else:
        raise ValueError(f"Invalid tool_name: {tool_name}")

    init_reference(genome_version)
    init_transcript_annotations(genome_version)
    init_database_connection()

    user_ip = get_user_ip(request)

    # check cache before processing the variant
    results = get_splicing_scores_from_cache(tool_name, variant, genome_version, distance_param, mask_param)
    duration = (datetime.now() - start_time).total_seconds()
    if results:
        log(f"{tool_name}:from-cache", ip=user_ip, variant=variant, genome=genome_version, distance=distance_param, mask=mask_param)
    else:
        error_message = exceeds_rate_limit(user_ip)
        if error_message:
            print(f"{logging_prefix}: {user_ip}: response: {error_message}", flush=True)
            return error_response(error_message, source=tool_name)

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

        duration = (datetime.now() - start_time).total_seconds()
        log(f"{tool_name}:computed", ip=user_ip, duration=duration, variant=variant, genome=genome_version, distance=distance_param, mask=mask_param)

        if "error" not in results:
            add_splicing_scores_to_cache(tool_name, variant, genome_version, distance_param, mask_param, results)

    if "error" in results:
        log(f"{tool_name}:error", ip=user_ip, variant=variant, genome=genome_version, distance=distance_param, mask=mask_param, details=results["error"])

    response_json = {}
    response_json.update(params)  # copy input params to output
    response_json.update(results)

    print(f"{logging_prefix}: {variant} took {str(datetime.now() - start_time)}", flush=True)

    return Response(json.dumps(response_json), status=200, mimetype='application/json', headers=[
        ('Access-Control-Allow-Origin', '*'),
    ])


def log(event_name, ip=None, duration=None, variant=None, genome=None, distance=None, mask=None, details=None):
    """Utility method for logging an event"""

    try:
        if duration is not None: duration = float(duration)
        if distance is not None: distance = int(distance)
        if mask is not None: mask = int(mask)
    except Exception as e:
        print(f"Error parsing log params: {e}", flush=True)
        return

    try:
        run_sql(r"INSERT INTO log (event_name, ip, duration, variant, genome, distance, mask, details) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)",
                (event_name, ip, duration, variant, genome, distance, mask, details))
    except Exception as e:
        print(f"Log error: {e}", flush=True)

def get_user_ip(request):
    return request.environ.get("HTTP_X_FORWARDED_FOR")

@app.route('/log/<string:name>/', strict_slashes=False)
def log_event(name):
    init_database_connection()

    if DATABASE_CONNECTION is None:
        message = f"Log error: Database not available"
        print(message, flush=True)
        return error_response(f"ERROR: {message}")

    if name != "show_igv":
        message = f"Log error: invalid event name: {name}"
        print(message, flush=True)
        return error_response(f"ERROR: {message}")

    # check params
    params = {}
    if request.values:
        params.update(request.values)
    if not params:
        params.update(request.get_json(force=True, silent=True) or {})

    variant = params.get("variant")
    genome_version = params.get("hg")
    distance_param = params.get("distance")
    mask_param = params.get("mask")
    details = params.get("details")
    if details:
        details = str(details)
        details = details[:2000]

    user_ip = get_user_ip(request)
    logging_prefix = datetime.now().strftime("%m/%d/%Y %H:%M:%S") + f" {user_ip} t{os.getpid()}"
    print(f"{logging_prefix}: ======================", flush=True)
    print(f"{logging_prefix}: {variant} show igv with hg={genome_version}, "
          f"distance={distance_param}, mask={mask_param}", flush=True)

    log(name,
        ip=user_ip,
        variant=variant,
        genome=genome_version,
        distance=distance_param,
        mask=mask_param,
        details=details)

    return Response(json.dumps({"status": "Done"}), status=200, mimetype='application/json', headers=[
        ('Access-Control-Allow-Origin', '*'),
    ])

@app.route('/', strict_slashes=False, defaults={'path': ''})
@app.route('/<path:path>/')
def catch_all(path):
    return f"SpliceAI-lookup APIs: invalid endpoint {path}"

app.run(debug=DEBUG, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))
