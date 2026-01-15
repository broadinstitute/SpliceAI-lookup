from datetime import datetime
import gzip
import json
import logging
import os
import psycopg2
import re
import time
import traceback


# used for DB connection pooling
from psycopg2.pool import SimpleConnectionPool
from contextlib import contextmanager

# flask imports
from flask import Flask, g, request, Response, send_from_directory
from flask_cors import CORS
from flask_talisman import Talisman

# SAI-10k style predictions for splice consequences
from sai10k_predictions import sai10k_get_transcript_predictions

app = Flask(__name__)

CORS(app)


DEBUG = True # if socket.gethostname() == "spliceai-lookup" else True
if not DEBUG:
    Talisman(app)

logging.getLogger('werkzeug').disabled = True

DEFAULT_DISTANCE = 500  # maximum distance between the variant and gained/lost splice site, defaults to 500
MAX_DISTANCE_LIMIT = 10000
DEFAULT_MASK = 0        # mask scores representing annotated acceptor/donor gain and unannotated acceptor/donor loss, defaults to 0

SPLICEAI_EXAMPLE_URL = f"/spliceai/?hg=38&distance=500&mask=0&variant=chr8-140300615-C-G&bc=basic"
PANGOLIN_EXAMPLE_URL = f"/pangolin/?hg=38&distance=500&mask=0&variant=chr8-140300615-C-G&bc=basic"


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

GENCODE_VERSION = "v49"

SHARED_TRANSCRIPT_ANNOTATIONS = {}
SHARED_TRANSCRIPT_ANNOTATION_PATHS = {
    ("37", "basic"): f"/gencode.{GENCODE_VERSION}lift37.basic.annotation.transcript_annotations.json.gz",
    ("38", "basic"): f"/gencode.{GENCODE_VERSION}.basic.annotation.transcript_annotations.json.gz",
    ("37", "comprehensive"): f"/gencode.{GENCODE_VERSION}lift37.annotation.transcript_annotations.json.gz",
    ("38", "comprehensive"): f"/gencode.{GENCODE_VERSION}.annotation.transcript_annotations.json.gz",
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
        ("37", "basic"): f"/gencode.{GENCODE_VERSION}lift37.basic.annotation.txt.gz",
        ("38", "basic"): f"/gencode.{GENCODE_VERSION}.basic.annotation.txt.gz",
        ("37", "comprehensive"): f"/gencode.{GENCODE_VERSION}lift37.annotation.txt.gz",
        ("38", "comprehensive"): f"/gencode.{GENCODE_VERSION}.annotation.txt.gz",
    }

elif TOOL == "pangolin":
    from pkg_resources import resource_filename
    from pangolin.pangolin import process_variant as process_variant_using_pangolin
    from pangolin.model import torch, Pangolin, L, W, AR
    import gffutils

    PANGOLIN_ANNOTATION_PATHS = {
        ("37", "basic"): f"/gencode.{GENCODE_VERSION}lift37.basic.annotation.without_chr_prefix.db",
        ("38", "basic"): f"/gencode.{GENCODE_VERSION}.basic.annotation.db",
        ("37", "comprehensive"): f"/gencode.{GENCODE_VERSION}lift37.annotation.without_chr_prefix.db",
        ("38", "comprehensive"): f"/gencode.{GENCODE_VERSION}.annotation.db",
    }
else:
    raise ValueError(f'Environment variable "TOOL" should be set to either "spliceai" or "pangolin" instead of: "{os.environ.get("TOOL")}"')


RATE_LIMIT_ERROR_MESSAGE = (
    f"Rate limit exceeded. This server only supports interactive use. To process large numbers of variants programmatically, "
    f"please install a local instance of the API server, or just run the prediction models directly. Attempts to query large "
    f"numbers of variants programmatically will result in loss of access to this API for an extended period of time. Contact "
    f"us at https://github.com/broadinstitute/SpliceAI-lookup/issues if you have any questions."
)


def init_spliceai(genome_version, basic_or_comprehensive):
    
    if (genome_version, basic_or_comprehensive) not in SPLICEAI_ANNOTATOR:
        SPLICEAI_ANNOTATOR[(genome_version, basic_or_comprehensive)] = Annotator(
            FASTA_PATH[genome_version],
            SPLICEAI_ANNOTATION_PATHS[(genome_version, basic_or_comprehensive)]
        )


def init_transcript_annotations(genome_version, basic_or_comprehensive):
    if (genome_version, basic_or_comprehensive) in SHARED_TRANSCRIPT_ANNOTATIONS:
        return

    # init shared transcript annotations
    with gzip.open(SHARED_TRANSCRIPT_ANNOTATION_PATHS[(genome_version, basic_or_comprehensive)], "rt") as ta_f:
        SHARED_TRANSCRIPT_ANNOTATIONS[(genome_version, basic_or_comprehensive)] = json.load(ta_f)


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


#while True:
#    # https://groups.google.com/g/google-cloud-sql-discuss/c/mxsaf-YDrbA?pli=1
#    # https://cloud.google.com/sql/docs/postgres/flags#gcloud
#
#    error_count = 0
#    try:
#        DATABASE_CONNECTION_POOL = SimpleConnectionPool(
#            minconn=1,
#            maxconn=5,
#            dbname="spliceai-lookup-db",
#            user="postgres",
#            password=os.environ.get("DB_PASSWORD"),
#            host="/cloudsql/spliceai-lookup-412920:us-central1:spliceai-lookup-db",
#            port="5432",
#            connect_timeout=5,
#        )
#        print(f"Successfully connected to database", flush=True)
#        break
#    except psycopg2.Error as e:
#        error_count += 1
#        time.sleep(2)
#        print(f"Error connecting to database: {e}", flush=True)
#        traceback.print_exc()
#        if error_count > 5:
#            print(f"Error connecting to database. Exiting...", flush=True)
#            sys.exit(1)


@contextmanager
def get_db_connection():
    """Get a database connection"""
    #conn = DATABASE_CONNECTION_POOL.getconn()
    try:
        conn = psycopg2.connect(
            dbname="spliceai-lookup-db",
            user="postgres",
            password=os.environ.get("DB_PASSWORD"),
            host="/cloudsql/spliceai-lookup-412920:us-central1:spliceai-lookup-db",
            port="5432",
            connect_timeout=5,
        )
    except Exception as e:
        print(f"ERROR: Unable to connect to SQL database: {e}")
        conn = None

    try:
        yield conn
    finally:
        if conn is not None:
            conn.close()
            #DATABASE_CONNECTION_POOL.putconn(conn)

@contextmanager
def get_db_cursor(conn):
    """Get a database cursor"""
    if conn is None:
        return

    cursor = conn.cursor()
    try:
        yield cursor
        conn.commit()
    finally:
        cursor.close()


def run_sql(conn, sql_query, *params):
    if conn is None:
        return

    with get_db_cursor(conn) as cursor:
        cursor.execute(sql_query, *params)
        try:
            results = cursor.fetchall()
        except:
            results = []
    return results


def get_transcript_structure(conn, transcript_id, genome_version):
    """Retrieve transcript structure from the database.

    Args:
        conn: Database connection
        transcript_id: Transcript ID (with or without version, e.g., "ENST00000123456" or "ENST00000123456.5")
        genome_version: "37" or "38"

    Returns:
        dict with transcript structure fields, or None if not found:
            - EXON_STARTS: list of 1-based exon start positions
            - EXON_ENDS: list of 1-based exon end positions
            - CDS_START: 1-based CDS start position (or None if non-coding)
            - CDS_END: 1-based CDS end position (or None if non-coding)
            - EXON_FRAMES: list of frame values for each exon (or None)
            - STRAND: '+' or '-'
    """
    if conn is None:
        return None

    # Remove version suffix if present
    transcript_id_without_version = transcript_id.split(".")[0]

    table_name = f"transcripts_hg{genome_version}"
    rows = run_sql(
        conn,
        f"""SELECT strand, cds_start, cds_end, exon_starts, exon_ends, exon_frames
           FROM {table_name} WHERE transcript_id = %s""",
        (transcript_id_without_version,)
    )

    if not rows:
        return None

    strand, cds_start, cds_end, exon_starts_str, exon_ends_str, exon_frames_str = rows[0]

    # Parse comma-separated values and convert from 0-based to 1-based coordinates
    # genePred format uses 0-based start and 0-based half-open end
    exon_starts_0based = [int(s) for s in exon_starts_str.rstrip(",").split(",") if s]
    exon_ends_0based = [int(s) for s in exon_ends_str.rstrip(",").split(",") if s]

    # Convert to 1-based coordinates (start +1, end stays same since half-open becomes closed)
    exon_starts_1based = [s + 1 for s in exon_starts_0based]
    exon_ends_1based = exon_ends_0based  # Already correct for 1-based closed interval

    # Convert CDS coordinates (0-based half-open to 1-based closed)
    cds_start_1based = cds_start + 1 if cds_start is not None else None
    cds_end_1based = cds_end if cds_end is not None else None

    # Parse exon frames
    exon_frames = None
    if exon_frames_str:
        exon_frames = [int(f) for f in exon_frames_str.rstrip(",").split(",") if f]

    return {
        "EXON_STARTS": exon_starts_1based,
        "EXON_ENDS": exon_ends_1based,
        "CDS_START": cds_start_1based,
        "CDS_END": cds_end_1based,
        "EXON_FRAMES": exon_frames,
        "STRAND": strand,
    }


#def does_table_exist(table_name):
#    results = run_sql(f"SELECT EXISTS (SELECT 1 AS result FROM pg_tables WHERE tablename=%s)", (table_name,))
#    does_table_already_exist = results[0][0]
#    return does_table_already_exist

#if not does_table_exist("cache"):
#    print("Creating cache table")
#    run_sql("""CREATE TABLE cache (key TEXT UNIQUE, value TEXT, counter INT, accessed TIMESTAMP DEFAULT now())""")
#    run_sql("""CREATE INDEX cache_index ON cache (key)""")

#if not does_table_exist("log"):
#    print("Creating event_log table")
#    run_sql("""CREATE TABLE log (event_name TEXT, ip TEXT, logtime TIMESTAMP DEFAULT now(), duration REAL, variant TEXT, genome VARCHAR(10), bc VARCHAR(20), distance INT, mask INT4, details TEXT, variant_consequence TEXT)""")
#    run_sql("""CREATE INDEX idx_log_ip_logtime ON log USING btree (ip, logtime DESC)""")
#    run_sql("""CREATE INDEX idx_log_event_name ON log USING btree (event_name)""")

#if not does_table_exist("restricted_ips"):
#    print("Creating restricted_ips table")
#    run_sql("""CREATE TABLE restricted_ips (ip TEXT UNIQUE, created TIMESTAMP DEFAULT now())""")
#    run_sql("""CREATE INDEX idx_restricted_ips_created ON restricted_ips USING btree (created)""")

# Query to add ip to the restricted_ips table
#run_sql("""INSERT INTO restricted_ips (ip) VALUES ('210.3.222.157')""")

def is_user_on_whitelist(conn, user_ip):
    """Check if the user is on the whitelist"""
    if conn is None or not user_ip:
        return False

    if not re.match(r"^\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3}$", user_ip):
        return False

    rows = run_sql(conn, "SELECT COUNT(ip) FROM whitelist_ips WHERE ip=%s", (user_ip,))
    return rows and int(rows[0][0]) > 0

def exceeds_rate_limit(conn, user_ip, params):
    """Rate limit requests based on user ip address"""

    #"""
    #SELECT * FROM log WHERE event_name like '%computed' AND duration > 2 AND ip='210.3.222.157' AND logtime >= NOW() - INTERVAL '5 minutes' ;
    #SELECT ip, count(*) FROM log WHERE event_name like '%computed' AND duration > 2 AND logtime >= NOW() - INTERVAL '20 minutes' GROUP BY ip ORDER BY count DESC;
    #"""

    try:
        if conn is None:
            return False

        if is_user_on_whitelist(conn, params.get("ip")):
            return False

        # check if the user has exceeded the rate limit or is on the list of restricted IPs
        rows = run_sql(conn, "SELECT COUNT(ip) FROM restricted_ips WHERE ip=%s AND created >= NOW() - INTERVAL '1 weeks'", (user_ip,))
        is_user_currently_blocked = rows and int(rows[0][0]) > 0
        if is_user_currently_blocked:
            return RATE_LIMIT_ERROR_MESSAGE

        rows = run_sql(conn, "SELECT COUNT(ip) FROM log WHERE event_name LIKE %s AND ip=%s AND logtime >= NOW() - INTERVAL '7 minutes'", ("%computed%", user_ip))
        did_user_exceed_rate_limit = rows and int(rows[0][0]) >= 50
        if did_user_exceed_rate_limit and not is_user_on_whitelist(conn, user_ip):
            # the user has exceeded the rate limit: computing scores for 50 or more variants in the last 7 minutes
            rows = run_sql(conn, "SELECT COUNT(ip) FROM log WHERE event_name='rate_limit_exceeded' AND ip=%s AND logtime >= NOW() - INTERVAL '5 minutes'", (user_ip,))
            user_hit_rate_limit_exceeded_recently = rows and int(rows[0][0]) > 0
            if not user_hit_rate_limit_exceeded_recently:
                # the user will receive at most one "rate_limit_exceeded" event every 5 minutes
                log(conn, f"rate_limit_exceeded", ip=user_ip)
                rows = run_sql(conn, "SELECT COUNT(ip) FROM log WHERE event_name='rate_limit_exceeded' AND ip=%s AND logtime >= NOW() - INTERVAL '1 days'", (user_ip,))
                user_triggered_too_many_rate_limit_exceeded_errors_today = rows and int(rows[0][0]) >= 5
                if user_triggered_too_many_rate_limit_exceeded_errors_today:
                    # the user has hit the limit of 5 or more "rate_limit_exceeded" events during the last 24 hours
                    rows = run_sql(conn, "SELECT COUNT(ip) FROM restricted_ips WHERE ip=%s", (user_ip,))
                    need_to_delete_previous_restricted_ip_record = rows and int(rows[0][0]) > 0
                    if need_to_delete_previous_restricted_ip_record:
                        # delete the previous record
                        run_sql(conn, "DELETE FROM restricted_ips WHERE ip=%s", (user_ip,))

                    # block the user's IP for 1 week
                    run_sql(conn, "INSERT INTO restricted_ips (ip) VALUES (%s)", (user_ip,))

            return RATE_LIMIT_ERROR_MESSAGE

    except Exception as e:
        print(f"Error while checking rate limit: {e}", flush=True)
        # print traceback
        traceback.print_exc()
        return False


def get_splicing_scores_cache_key(tool_name, variant, genome_version, distance, mask, basic_or_comprehensive="basic"):
    return f"{tool_name}__{variant}__hg{genome_version}__d{distance}__m{mask}__{basic_or_comprehensive}"


def get_splicing_scores_from_cache(conn, tool_name, variant, genome_version, distance, mask, basic_or_comprehensive="basic"):
    results = {}
    key = get_splicing_scores_cache_key(tool_name, variant, genome_version, distance, mask, basic_or_comprehensive)
    try:
        rows = run_sql(conn, f"SELECT value FROM cache WHERE key=%s", (key,))
        if rows:
            results = json.loads(rows[0][0])
            results["source"] += ":cache"
    except Exception as e:
        print(f"Cache error: {e}", flush=True)

    return results


def add_splicing_scores_to_cache(conn, tool_name, variant, genome_version, distance, mask, basic_or_comprehensive, results):
    key = get_splicing_scores_cache_key(tool_name, variant, genome_version, distance, mask, basic_or_comprehensive)
    try:
        results_string = json.dumps(results)

        run_sql(conn,
                r"""INSERT INTO cache (key, value, counter, accessed) VALUES (%s, %s, 1, now()) """ +
                r"""ON CONFLICT (key) DO """ +
                r"""UPDATE SET key=%s, value=%s, counter=cache.counter+1, accessed=now()""", (key, results_string, key, results_string))
    except Exception as e:
        print(f"Cache error: {e}", flush=True)


def get_spliceai_scores(conn, variant, genome_version, distance_param, mask_param, basic_or_comprehensive_param):
    try:
        chrom, pos, ref, alt = parse_variant(variant)
    except ValueError as e:
        return {
            "variant": variant,
            "source": "spliceai",
            "error": str(e),
        }

    # generate error message if variant falls outside annotated exons or introns
    record = VariantRecord(chrom, pos, ref, alt)
    try:
        scores = get_delta_scores(
            record,
            SPLICEAI_ANNOTATOR[(genome_version, basic_or_comprehensive_param)],
            distance_param,
            mask_param)
    except Exception as e:
        print(f"ERROR while computing SpliceAI scores for {variant}: {e}")
        traceback.print_exc()
        return {
            "variant": variant,
            "source": "spliceai",
            "error": f"{type(e)}: {e}",
        }

    if not scores:
        return {
            "variant": variant,
            "source": "spliceai",
            "error": f"The SpliceAI model did not return any scores for {variant}. This may be because the variant does "
                     f"not overlap any exons or introns defined by the GENCODE '{basic_or_comprehensive_param}' annotation.",
        }

    #scores = [s[s.index("|")+1:] for s in scores]  # drop allele field

    # to reduce the response size, return all non-zero scores only for the canonical transcript (or the 1st transcript)
    all_non_zero_scores = None
    all_non_zero_scores_strand = None
    all_non_zero_scores_transcript_id = None
    all_non_zero_scores_transcript_priority = -1
    max_delta_score_sum = 0
    selected_transcript = None
    for i, transcript_scores in enumerate(scores):
        if "ALL_NON_ZERO_SCORES" not in transcript_scores:
            continue

        transcript_id_without_version = transcript_scores.get("NAME", "").split(".")[0]

        # get json annotations for this transcript
        transcript_annotations = SHARED_TRANSCRIPT_ANNOTATIONS[(genome_version, basic_or_comprehensive_param)].get(transcript_id_without_version)
        if transcript_annotations is None:
            raise ValueError(f"Missing annotations for {transcript_id_without_version} in {genome_version} annotations")

        # add the extra transcript annotations from the json file to the transcript scores dict
        transcript_scores.update(transcript_annotations)

        # add transcript structure from the database for heuristic predictions
        transcript_structure = get_transcript_structure(conn, transcript_id_without_version, genome_version)
        if transcript_structure:
            transcript_scores.update(transcript_structure)

        # decide whether to use ALL_NON_ZERO_SCORES from this transcript
        current_transcript_priority = TRANSCRIPT_PRIORITY_ORDER[transcript_annotations["t_priority"]]
        current_delta_score_sum = sum(abs(float(transcript_scores[key])) for key in ("DS_AG", "DS_AL", "DS_DG", "DS_DL"))
        if current_transcript_priority > all_non_zero_scores_transcript_priority:
            max_delta_score_sum = current_delta_score_sum
            all_non_zero_scores_transcript_priority = current_transcript_priority
            all_non_zero_scores = transcript_scores["ALL_NON_ZERO_SCORES"]
            all_non_zero_scores_strand = transcript_scores["t_strand"]
            all_non_zero_scores_transcript_id = transcript_scores["t_id"]
            selected_transcript = transcript_scores

        elif current_transcript_priority == all_non_zero_scores_transcript_priority and current_delta_score_sum > max_delta_score_sum:
            # select the one with the highest delta score sum
            max_delta_score_sum = current_delta_score_sum
            all_non_zero_scores = transcript_scores["ALL_NON_ZERO_SCORES"]
            all_non_zero_scores_strand = transcript_scores["t_strand"]
            all_non_zero_scores_transcript_id = transcript_scores["t_id"]
            selected_transcript = transcript_scores

        for redundant_key in "ALLELE", "NAME", "STRAND", "ALL_NON_ZERO_SCORES":
            del transcript_scores[redundant_key]

    # Compute heuristic predictions for the selected transcript
    heuristic_predictions = None
    if selected_transcript:
        try:
            heuristic_predictions = sai10k_get_transcript_predictions(selected_transcript, pos)
        except Exception as e:
            print(f"WARNING: Error computing heuristic predictions for {variant}: {e}")
            traceback.print_exc()

    return {
        "variant": variant,
        "genomeVersion": genome_version,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "distance": distance_param,
        "mask": mask_param,
        "scores": scores,
        "source": "spliceai:model",
        "allNonZeroScores": all_non_zero_scores,
        "allNonZeroScoresStrand": all_non_zero_scores_strand,
        "allNonZeroScoresTranscriptId": all_non_zero_scores_transcript_id,
        "heuristicPredictions": heuristic_predictions,
    }


def get_pangolin_scores(variant, genome_version, distance_param, mask_param, basic_or_comprehensive_param):
    if genome_version not in ("37", "38"):
        raise ValueError(f"Invalid genome_version: {mask_param}")

    if mask_param not in ("True", "False"):
        raise ValueError(f"Invalid mask_param: {mask_param}")

    if basic_or_comprehensive_param not in ("basic", "comprehensive"):
        raise ValueError(f"Invalid basic_or_comprehensive_param: {basic_or_comprehensive_param}")

    try:
        chrom, pos, ref, alt = parse_variant(variant)
    except ValueError as e:
        print(f"ERROR while parsing variant {variant}: {e}")
        traceback.print_exc()

        return {
            "variant": variant,
            "source": "pangolin",
            "error": str(e),
        }

    if len(ref) > 1 and len(alt) > 1:
        return {
            "variant": variant,
            "source": "pangolin",
            "error": f"Pangolin does not currently support complex InDels like {chrom}-{pos}-{ref}-{alt}",
        }

    class PangolinArgs:
        reference_file = FASTA_PATH[genome_version]
        distance = distance_param
        mask = mask_param
        score_cutoff = None
        score_exons = "False"

    pangolin_models = []

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
            pangolin_models.append(model)

    features_db = gffutils.FeatureDB(PANGOLIN_ANNOTATION_PATHS[(GENOME_VERSION, basic_or_comprehensive_param)])
    scores = process_variant_using_pangolin(
        0, chrom, int(pos), ref, alt, features_db, pangolin_models, PangolinArgs)

    if not scores:
        return {
            "variant": variant,
            "source": "pangolin",
            "error": f"Pangolin was unable to compute scores for this variant",
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
        transcript_annotations = SHARED_TRANSCRIPT_ANNOTATIONS[(genome_version, basic_or_comprehensive_param)].get(transcript_id_without_version)
        if transcript_annotations is None:
            raise ValueError(f"Missing annotations for {transcript_id_without_version} in {genome_version} annotations")

        # add the extra transcript annotations from the json file to the transcript scores dict
        transcript_scores.update(transcript_annotations)

        # decide whether to use ALL_NON_ZERO_SCORES from this gene
        current_delta_score_sum = sum(abs(float(s.get("SG_ALT", 0)) - float(s.get("SG_REF", 0)))
                              for s in transcript_scores["ALL_NON_ZERO_SCORES"])
        current_delta_score_sum += sum(abs(float(s.get("SL_ALT", 0)) - float(s.get("SL_REF", 0)))
                               for s in transcript_scores["ALL_NON_ZERO_SCORES"])

        # return all_non_zero_scores for the transcript or gene with the highest delta score sum
        if current_delta_score_sum > max_delta_score_sum:
            all_non_zero_scores = transcript_scores["ALL_NON_ZERO_SCORES"]
            all_non_zero_scores_strand = transcript_scores["STRAND"]
            all_non_zero_scores_transcript_id = transcript_scores["NAME"]
            max_delta_score_sum = current_delta_score_sum

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
        "mask": mask_param,
        "scores": scores,
        "source": "pangolin:model",
        "allNonZeroScores": all_non_zero_scores,
        "allNonZeroScoresStrand": all_non_zero_scores_strand,
        "allNonZeroScoresTranscriptId": all_non_zero_scores_transcript_id,
    }


@app.route("/spliceai/", methods=['POST', 'GET'])
def run_spliceai():
    with get_db_connection() as conn:
        return run_splice_prediction_tool(conn, tool_name="spliceai")


@app.route("/pangolin/", methods=['POST', 'GET'])
def run_pangolin():
    with get_db_connection() as conn:
        return run_splice_prediction_tool(conn, tool_name="pangolin")


def run_splice_prediction_tool(conn, tool_name):
    """Handles API request for splice prediction

    Args:
        conn (psycopg2.connection): Database connection
        tool_name (str): "spliceai" or "pangolin"
    """

    if tool_name != TOOL:
        return error_response(f"ERROR: This server is configured to run {TOOL} rather than {tool_name}.\n", source=tool_name)

    user_ip = get_user_ip(request)


    start_time = datetime.now()
    #logging_prefix = start_time.strftime("%m/%d/%Y %H:%M:%S") + f" t{os.getpid()} ip:{user_ip}"
    logging_prefix = f"t{os.getpid()} ip:{user_ip}"
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

    basic_or_comprehensive_param = params.get("bc", "basic")
    if basic_or_comprehensive_param not in ("basic", "comprehensive"):
        return error_response(f'Invalid "bc" value: "{basic_or_comprehensive_param}". The value must be either "basic" or "comprehensive". For example: {example_url}\n', source=tool_name)

    variant_consequence = params.get("variant_consequence")

    force = params.get("force")  # ie. don't use cache

    print(f"{logging_prefix}: ======================", flush=True)
    print(f"{logging_prefix}: {variant} tool={tool_name} hg={genome_version}, distance={distance_param}, mask={mask_param}, bc={basic_or_comprehensive_param}", flush=True)

    if tool_name == "spliceai":
        init_spliceai(genome_version, basic_or_comprehensive_param)

    init_transcript_annotations(genome_version, basic_or_comprehensive_param)

    # check cache before processing the variant
    results = {}
    if not force:
        results = get_splicing_scores_from_cache(conn, tool_name, variant, genome_version, distance_param, mask_param, basic_or_comprehensive_param)

    duration = (datetime.now() - start_time).total_seconds()
    if results:
        log(conn, f"{tool_name}:from-cache", ip=user_ip, variant=variant, genome=genome_version, distance=distance_param, mask=mask_param, bc=basic_or_comprehensive_param, variant_consequence=variant_consequence)
    else:
        error_message = exceeds_rate_limit(conn, user_ip, params)
        if error_message:
            print(f"{logging_prefix}: {user_ip}: response: {error_message}", flush=True)
            return error_response(error_message, source=tool_name)

        try:
            if tool_name == "spliceai":
                results = get_spliceai_scores(conn, variant, genome_version, distance_param, int(mask_param), basic_or_comprehensive_param)
            elif tool_name == "pangolin":
                pangolin_mask_param = "True" if mask_param == "1" else "False"
                results = get_pangolin_scores(variant, genome_version, distance_param, pangolin_mask_param, basic_or_comprehensive_param)
            else:
                raise ValueError(f"Invalid tool_name: {tool_name}")
        except Exception as e:
            traceback.print_exc()
            return error_response(f"ERROR: {e}", source=tool_name)

        duration = (datetime.now() - start_time).total_seconds()
        log(conn, f"{tool_name}:computed", ip=user_ip, duration=duration, variant=variant, genome=genome_version, distance=distance_param, mask=mask_param, bc=basic_or_comprehensive_param, variant_consequence=variant_consequence)

        if "error" not in results:
            add_splicing_scores_to_cache(conn, tool_name, variant, genome_version, distance_param, mask_param, basic_or_comprehensive_param, results)

    if "error" in results:
        log(conn, f"{tool_name}:error", ip=user_ip, variant=variant, genome=genome_version, distance=distance_param, mask=mask_param, details=results["error"], bc=basic_or_comprehensive_param, variant_consequence=variant_consequence)

    response_json = {}
    response_json.update(params)  # copy input params to output
    response_json.update(results)

    response_log_string = ", ".join([f"{k}: {v}" for k, v in response_json.items() if not k.startswith("allNonZeroScores")])
    print(f"{logging_prefix}: {variant} response took {str(datetime.now() - start_time)}: {response_log_string}", flush=True)

    return Response(json.dumps(response_json), status=200, mimetype='application/json', headers=[
        ('Access-Control-Allow-Origin', '*'),
    ])


def log(conn, event_name, ip=None, duration=None, variant=None, genome=None, distance=None, mask=None, bc=None, details=None, variant_consequence=None):
    """Utility method for logging an event"""

    try:
        if duration is not None: duration = float(duration)
        if distance is not None: distance = int(distance)
        if mask is not None: mask = int(mask)
    except Exception as e:
        print(f"Error parsing log params: {e}", flush=True)
        return

    try:
        run_sql(conn,
                r"INSERT INTO log (event_name, ip, duration, variant, genome, distance, mask, bc, details, variant_consequence) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
                (event_name, ip, duration, variant, genome, distance, mask, bc, details, variant_consequence))
    except Exception as e:
        print(f"Log error: {e}", flush=True)


def get_user_ip(request):
    return request.environ.get("HTTP_X_FORWARDED_FOR")


@app.route('/log/<string:name>/', strict_slashes=False)
def log_event(name):

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
    basic_or_comprehensive_param = params.get("bc")
    details = params.get("details")
    variant_consequence = params.get("variant_consequence")
    if details:
        details = str(details)
        details = details[:2000]

    user_ip = get_user_ip(request)
    logging_prefix = datetime.now().strftime("%m/%d/%Y %H:%M:%S") + f" {user_ip} t{os.getpid()}"
    print(f"{logging_prefix}: ======================", flush=True)
    print(f"{logging_prefix}: {variant} show igv with hg={genome_version}, distance={distance_param}, mask={mask_param}", flush=True)

    with get_db_connection() as conn:
        log(conn,
            name,
            ip=user_ip,
            variant=variant,
            genome=genome_version,
            distance=distance_param,
            mask=mask_param,
            bc=basic_or_comprehensive_param,
            details=details,
            variant_consequence=variant_consequence)

    return Response(json.dumps({"status": "Done"}), status=200, mimetype='application/json', headers=[
        ('Access-Control-Allow-Origin', '*'),
    ])


@app.route('/', strict_slashes=False, defaults={'path': ''})
@app.route('/<path:path>/')
def catch_all(path):
    return f"SpliceAI-lookup APIs: invalid endpoint {path}"


if '__main__' == __name__ or os.environ.get('RUNNING_ON_GOOGLE_CLOUD_RUN'):
    app.run(debug=DEBUG, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))
