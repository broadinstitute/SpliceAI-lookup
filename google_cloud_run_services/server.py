import time
_PROCESS_START_TIME = time.time()
from datetime import datetime
import gzip
import json
import logging
import os
import psycopg2
import re
import threading
import traceback


# used for DB connection pooling
from psycopg2.pool import ThreadedConnectionPool
from contextlib import contextmanager

# flask imports
from flask import Flask, g, request, Response
from flask_cors import CORS
from flask_talisman import Talisman

# SAI-10k-calc predictions for splice consequences
from sai10k_predictions import sai10k_get_transcript_predictions, sai10k_select_transcript, TRANSCRIPT_PRIORITY_ORDER

app = Flask(__name__)

# Intentional: this is a read-only public API with no cookie/session auth.
# Wildcard CORS lets notebook and cross-origin research tools call the API
# directly. Restricting to a single frontend origin would break those users
# without adding real security (no auth state to protect).
CORS(app)


# On Cloud Run, disable Werkzeug's debug PIN / interactive traceback; keep it
# on for local development.
DEBUG = not os.environ.get('RUNNING_ON_GOOGLE_CLOUD_RUN')

# Security headers: HSTS, CSP, X-Frame-Options, X-Content-Type-Options, etc.
# force_https=False because Cloud Run's load balancer terminates TLS and
# forwards plain HTTP to the container with X-Forwarded-Proto: https — the
# LB already enforces HTTPS at the edge, so an app-level redirect would loop
# with the LB.
Talisman(app, force_https=False)

logging.getLogger('werkzeug').disabled = True

DEFAULT_DISTANCE = 500  # maximum distance between the variant and gained/lost splice site, defaults to 500
MAX_DISTANCE_LIMIT = 10000
DEFAULT_MASK = 0        # mask scores representing annotated acceptor/donor gain and unannotated acceptor/donor loss, defaults to 0

SPLICEAI_EXAMPLE_URL = f"/spliceai/?hg=38&distance=500&mask=0&variant=chr8-140300615-C-G&bc=basic"
PANGOLIN_EXAMPLE_URL = f"/pangolin/?hg=38&distance=500&mask=0&variant=chr8-140300615-C-G&bc=basic"


VARIANT_RE = re.compile(
    r"(chr)?(?P<chrom>[0-9XYMTt]{1,2})"
    r"[-\s:]+"
    r"(?P<pos>[0-9]{1,9})"
    r"[-\s:]+"
    r"(?P<ref>[ACGT]+)"
    r"[-\s:>]+"
    r"(?P<alt>[ACGT]+)"
)

FASTA_PATH = {
    "37": "/hg19.fa.gz",
    "38": "/hg38.fa.gz",
}

# Lazy pyfastx Fasta singletons keyed by genome_version, used by SAI-10k-calc's
# premature-stop detection. Mirrors the SPLICEAI_ANNOTATOR cache pattern below.
# pyfastx (already in the container's spliceai/requirements.txt) handles
# bgzipped .fa.gz natively. Init failures are tolerated: detection silently
# falls back to None on every aberration, leaving the rest of the SAI-10k
# response intact.
SAI10K_FASTA = {}
_SAI10K_FASTA_LOCK = threading.Lock()


def _get_sai10k_fasta(genome_version):
    if genome_version in SAI10K_FASTA:
        return SAI10K_FASTA[genome_version]
    with _SAI10K_FASTA_LOCK:
        if genome_version in SAI10K_FASTA:
            return SAI10K_FASTA[genome_version]
        try:
            import pyfastx
            SAI10K_FASTA[genome_version] = pyfastx.Fasta(FASTA_PATH[genome_version])
        except Exception as e:
            print(f"WARNING: Failed to open FASTA for hg{genome_version} "
                  f"(SAI-10k premature-stop detection disabled): "
                  f"{type(e).__name__}: {e}")
            SAI10K_FASTA[genome_version] = None
    return SAI10K_FASTA[genome_version]

GENCODE_VERSION = "v49"

SHARED_TRANSCRIPT_ANNOTATIONS = {}
SHARED_TRANSCRIPT_ANNOTATION_PATHS = {
    ("37", "basic"): f"/gencode.{GENCODE_VERSION}lift37.basic.annotation.transcript_annotations.json.gz",
    ("38", "basic"): f"/gencode.{GENCODE_VERSION}.basic.annotation.transcript_annotations.json.gz",
    ("37", "comprehensive"): f"/gencode.{GENCODE_VERSION}lift37.annotation.transcript_annotations.json.gz",
    ("38", "comprehensive"): f"/gencode.{GENCODE_VERSION}.annotation.transcript_annotations.json.gz",
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
    # genome_version -> bare mito sequence name as it appears in the FASTA
    # (with any leading "chr" stripped, since spliceai.get_delta_scores calls
    # normalise_chrom() which re-adds it based on the fasta's first key).
    # Populated lazily by init_spliceai. hg19's fasta uses "MT", hg38's uses
    # "chrM"; without this remap, user-submitted "M"/"chrM" 500s on hg19.
    MITO_CHROM_NAME = {}
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
        t0 = time.time()
        print(f"[startup pid={os.getpid()}] init_spliceai({genome_version}, {basic_or_comprehensive}) START "
              f"+{t0 - _PROCESS_START_TIME:.2f}s", flush=True)
        SPLICEAI_ANNOTATOR[(genome_version, basic_or_comprehensive)] = Annotator(
            FASTA_PATH[genome_version],
            SPLICEAI_ANNOTATION_PATHS[(genome_version, basic_or_comprehensive)]
        )
        print(f"[startup pid={os.getpid()}] init_spliceai({genome_version}, {basic_or_comprehensive}) "
              f"Annotator() ready in {time.time() - t0:.2f}s", flush=True)
        if genome_version not in MITO_CHROM_NAME:
            keys = set(SPLICEAI_ANNOTATOR[(genome_version, basic_or_comprehensive)].ref_fasta.keys())
            for candidate in ('chrM', 'chrMT', 'MT', 'M'):
                if candidate in keys:
                    MITO_CHROM_NAME[genome_version] = candidate[3:] if candidate.startswith('chr') else candidate
                    break


def init_transcript_annotations(genome_version, basic_or_comprehensive):
    if (genome_version, basic_or_comprehensive) in SHARED_TRANSCRIPT_ANNOTATIONS:
        return

    t0 = time.time()
    with gzip.open(SHARED_TRANSCRIPT_ANNOTATION_PATHS[(genome_version, basic_or_comprehensive)], "rt") as ta_f:
        SHARED_TRANSCRIPT_ANNOTATIONS[(genome_version, basic_or_comprehensive)] = json.load(ta_f)
    print(f"[startup pid={os.getpid()}] init_transcript_annotations({genome_version}, {basic_or_comprehensive}) "
          f"loaded in {time.time() - t0:.2f}s (+{time.time() - _PROCESS_START_TIME:.2f}s since start)", flush=True)


def error_response(error_message, source=None, status=400):
    # Default to HTTP 400 (Bad Request) so clients and monitoring can
    # distinguish failures from successful responses. Callers that wrap an
    # internal exception pass status=500 explicitly. The previous default of
    # 200 made it impossible for the frontend's `xhr.status < 300` check —
    # and any generic API client — to detect the error.
    response_json = {"error": str(error_message)}
    if source:
        response_json["source"] = source
    return Response(json.dumps(response_json), status=status, mimetype='application/json')


def parse_variant(variant_str):
    match = VARIANT_RE.fullmatch(variant_str)
    if not match:
        raise ValueError(f"Unable to parse variant: {variant_str}")

    return match['chrom'], int(match['pos']), match['ref'], match['alt']


DB_CONNECT_KWARGS = dict(
    dbname="spliceai-lookup-db",
    user="postgres",
    password=os.environ.get("DB_PASSWORD"),
    host="/cloudsql/spliceai-lookup-412920:us-central1:spliceai-lookup-db",
    port="5432",
    connect_timeout=5,
)

# Module-level connection pool for Cloud SQL. Flask under Cloud Run typically
# serves multiple concurrent requests per instance via threaded workers, so use
# ThreadedConnectionPool (thread-safe) rather than SimpleConnectionPool.
# maxconn is sized to match Cloud Run's default per-instance request concurrency
# (80) so the 81st cache-miss request never blocks on getconn(). Each Cloud Run
# instance is one Python process; with ~10 instances the absolute peak is ~800
# Cloud SQL connections, well under the per-tier limit.
# If initialisation fails (e.g. transient Cloud SQL hiccup at startup, or
# DB_PASSWORD not set in a local dev environment), DATABASE_CONNECTION_POOL
# stays None and get_db_connection falls back to opening a connection per
# request. _try_init_database_pool retries the init lazily (throttled to once per minute)
# so a momentary outage at startup doesn't permanently disable pooling for the
# container's lifetime.
DATABASE_CONNECTION_POOL = None
_DATABASE_POOL_INIT_RETRY_SECONDS = 60
_database_pool_init_last_attempt = 0.0
_database_pool_init_lock = threading.Lock()


def _try_init_database_pool():
    """Attempt to (re)initialise DATABASE_CONNECTION_POOL, throttled to once per _DATABASE_POOL_INIT_RETRY_SECONDS."""
    global DATABASE_CONNECTION_POOL, _database_pool_init_last_attempt
    with _database_pool_init_lock:
        if DATABASE_CONNECTION_POOL is not None:
            return
        now = time.monotonic()
        if _database_pool_init_last_attempt and now - _database_pool_init_last_attempt < _DATABASE_POOL_INIT_RETRY_SECONDS:
            return
        _database_pool_init_last_attempt = now
        try:
            DATABASE_CONNECTION_POOL = ThreadedConnectionPool(minconn=1, maxconn=80, **DB_CONNECT_KWARGS)
            print("Successfully initialised DB connection pool", flush=True)
        except Exception as e:
            print(f"WARNING: DB connection pool init failed; falling back to per-request connections: {e}", flush=True)


_try_init_database_pool()


@contextmanager
def get_db_connection():
    """Get a database connection from the pool (or open a per-request connection if the pool is unavailable).

    Standard transaction discipline: commit when the with-block exits cleanly,
    rollback when an exception escapes. Previously this committed inside the
    cursor scope and then unconditionally rolled back on connection exit — an
    extra round-trip per request that cost real Cloud SQL latency under load.

    Broken connections (conn.closed != 0) are discarded instead of recycled.

    If the pool has not yet been initialised (e.g. Cloud SQL was unavailable at
    container startup), retry pool init lazily so a transient outage doesn't
    permanently disable pooling for this instance.
    """
    if DATABASE_CONNECTION_POOL is None:
        _try_init_database_pool()

    conn = None
    from_pool = False
    if DATABASE_CONNECTION_POOL is not None:
        try:
            conn = DATABASE_CONNECTION_POOL.getconn()
            from_pool = True
        except Exception as e:
            print(f"ERROR: Unable to get DB connection from pool: {e}")
            conn = None
    else:
        try:
            conn = psycopg2.connect(**DB_CONNECT_KWARGS)
        except Exception as e:
            print(f"ERROR: Unable to connect to SQL database: {e}")
            conn = None

    raised = False
    try:
        yield conn
    except Exception:
        raised = True
        raise
    finally:
        if conn is not None:
            try:
                if not conn.closed:
                    if raised:
                        conn.rollback()
                    else:
                        conn.commit()
            except Exception as txn_err:
                print(f"ERROR finalising DB transaction: {txn_err}", flush=True)
            if from_pool:
                try:
                    DATABASE_CONNECTION_POOL.putconn(conn, close=bool(conn.closed))
                except Exception as put_err:
                    print(f"ERROR returning connection to pool: {put_err}")
                    try:
                        conn.close()
                    except Exception:
                        pass
            else:
                try:
                    conn.close()
                except Exception:
                    pass

@contextmanager
def get_db_cursor(conn):
    """Yield a database cursor, or None when the connection is None.

    Transaction commit/rollback is handled at the connection scope (see
    get_db_connection), not here — yielding the cursor and closing it on exit
    is all this needs to do.

    A bare `return` before `yield` in a @contextmanager generator raises
    `RuntimeError("generator didn't yield")` when used via `with`, so yield
    None explicitly and let callers guard on the result.
    """
    if conn is None:
        yield None
        return

    cursor = conn.cursor()
    try:
        yield cursor
    finally:
        cursor.close()


def run_sql(conn, sql_query, *params):
    if conn is None:
        return []

    try:
        with get_db_cursor(conn) as cursor:
            cursor.execute(sql_query, *params)
            try:
                results = cursor.fetchall()
            except psycopg2.ProgrammingError:
                # No result set (e.g. from DELETE/INSERT/UPDATE); caller just needs [].
                results = []
    except psycopg2.Error:
        # Commit/rollback now happens at the connection scope (one commit per
        # `with get_db_connection()` block). A failed query leaves the conn in
        # an aborted-transaction state where any further query raises
        # InFailedSqlTransaction. Rollback here so the conn is usable for the
        # next query in the same scope, then re-raise so the caller can decide.
        try:
            if conn and not conn.closed:
                conn.rollback()
        except Exception:
            pass
        raise
    return results


def get_transcript_structures(conn, transcript_ids, genome_version):
    """Batch-fetch transcript structure for many transcripts in one round trip.

    Args:
        conn: Database connection.
        transcript_ids: iterable of transcript IDs WITHOUT version suffix
            (e.g. "ENST00000123456"). Caller is responsible for stripping ".N".
        genome_version: "37" or "38".

    Returns:
        dict mapping transcript_id -> structure dict (same fields the prior
        per-row helper produced):
            - EXON_STARTS: list of 1-based exon start positions
            - EXON_ENDS: list of 1-based exon end positions
            - CDS_START: 1-based CDS start position (or None if non-coding)
            - CDS_END: 1-based CDS end position (or None if non-coding)
            - STRAND: '+' or '-'
        Transcripts absent from the DB are simply missing from the result dict
        (caller distinguishes via `id in result`). Returns {} when conn is
        None or the input is empty. Returns None when the query raises so the
        caller can distinguish "DB unreachable mid-query" from "query
        succeeded with zero matches" and suppress per-row "not found"
        warnings.
    """
    if conn is None or not transcript_ids:
        return {}

    transcript_ids = list(transcript_ids)
    table_name = f"transcripts_hg{genome_version}"
    try:
        rows = run_sql(
            conn,
            f"""SELECT transcript_id, strand, cds_start, cds_end, exon_starts, exon_ends
               FROM {table_name} WHERE transcript_id = ANY(%s)""",
            (transcript_ids,)
        )
    except psycopg2.Error as e:
        # A transient OperationalError (e.g. broken connection mid-query)
        # returns None so SAI-10k falls back to annotation-based defaults
        # without spamming N per-transcript "not found" warnings.
        print(f"DB error fetching transcript structures for hg{genome_version}: {e}", flush=True)
        return None

    result = {}
    for transcript_id, strand, cds_start, cds_end, exon_starts_str, exon_ends_str in rows:
        # genePred uses 0-based half-open coordinates. Convert to 1-based closed.
        exon_starts_1based = [int(s) + 1 for s in exon_starts_str.rstrip(",").split(",") if s]
        exon_ends_1based = [int(s) for s in exon_ends_str.rstrip(",").split(",") if s]
        result[transcript_id] = {
            "EXON_STARTS": exon_starts_1based,
            "EXON_ENDS": exon_ends_1based,
            "CDS_START": cds_start + 1 if cds_start is not None else None,
            "CDS_END": cds_end if cds_end is not None else None,
            "STRAND": strand,
        }
    return result


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

    try:
        rows = run_sql(conn, "SELECT COUNT(ip) FROM whitelist_ips WHERE ip=%s", (user_ip,))
    except psycopg2.Error as e:
        # Fail closed (treat as not on whitelist) so a DB blip doesn't
        # accidentally skip the rate-limit check. exceeds_rate_limit's own
        # try/except already fails open if the rate-limit query itself fails.
        print(f"DB error checking whitelist for {user_ip}: {e}", flush=True)
        return False
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
        # Fail open so a transient DB hiccup doesn't lock everyone out, but log
        # loudly so repeated failures are visible in the Cloud Run logs — a
        # silent-always-allow would let an attacker DoS the DB to bypass the
        # rate limiter.
        print(f"SECURITY: rate-limit check failed (failing open): {e}", flush=True)
        traceback.print_exc()
        return False


# Bump SAI10K_VERSION whenever sai10k_predictions.py changes its classification
# logic or output shape, so cached responses from older algorithm versions are
# invalidated and recomputed.
SAI10K_VERSION = "v17"


def get_splicing_scores_cache_key(tool_name, variant, genome_version, distance, mask, basic_or_comprehensive="basic"):
    suffix = f"__sai10k-{SAI10K_VERSION}" if tool_name == "spliceai" else ""
    return f"{tool_name}__{variant}__hg{genome_version}__d{distance}__m{mask}__{basic_or_comprehensive}{suffix}"


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


def get_spliceai_scores(variant, genome_version, distance_param, mask_param, basic_or_comprehensive_param):
    try:
        chrom, pos, ref, alt = parse_variant(variant)
    except ValueError as e:
        return {
            "variant": variant,
            "source": "spliceai",
            "error": str(e),
        }

    # spliceai's normalise_chrom() handles "chr" prefix mismatches but not the
    # M↔MT alias, so a user submitting M/chrM against hg19 (which uses "MT")
    # would otherwise hit KeyError. Remap to whichever name the fasta uses.
    if chrom.upper() in {"M", "MT"} and genome_version in MITO_CHROM_NAME:
        chrom = MITO_CHROM_NAME[genome_version]

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

    # Enrich each transcript_scores with in-memory annotations (no DB).
    candidate_transcripts = []
    for transcript_scores in scores:
        if "ALL_NON_ZERO_SCORES" not in transcript_scores:
            continue

        transcript_id_without_version = transcript_scores.get("NAME", "").split(".")[0]

        transcript_annotations = SHARED_TRANSCRIPT_ANNOTATIONS[(genome_version, basic_or_comprehensive_param)].get(transcript_id_without_version)
        if transcript_annotations is None:
            raise ValueError(f"Missing annotations for {transcript_id_without_version} in {genome_version} annotations")
        transcript_scores.update(transcript_annotations)

        candidate_transcripts.append(transcript_scores)

    # Brief DB scope: enrich every candidate with transcript structure (used by
    # SAI-10k-calc and returned in the response) via a single batched SELECT.
    # Held only for that one query — not across model inference, which the
    # caller already ran outside any pooled connection.
    skip_cache = False
    db_enrich_t0 = time.perf_counter()
    structures = {}
    db_unavailable = False
    if candidate_transcripts:
        candidate_ids = [
            transcript_scores.get("NAME", "").split(".")[0]
            for transcript_scores in candidate_transcripts
        ]
        with get_db_connection() as conn:
            if conn is None:
                # DB unavailable: structures dict stays empty, SAI-10k-calc
                # silently falls back to its annotation defaults instead of
                # EXON_STARTS/EXON_ENDS/CDS_*/STRAND, and the degraded result
                # would be cached under the same key. Log loudly and skip the
                # cache write so the next request retries.
                print(f"WARNING: DB unavailable for transcript-structure lookup for {variant}; "
                      f"SAI-10k will use annotation defaults and the result will not be cached.", flush=True)
                skip_cache = True
                db_unavailable = True
            else:
                structures = get_transcript_structures(conn, candidate_ids, genome_version)
                if structures is None:
                    # Query failed mid-flight. Treat the same as a missing
                    # connection: log once, skip the cache, and suppress the
                    # per-row "not found" warnings below.
                    print(f"WARNING: DB query for transcript-structure lookup failed for {variant}; "
                          f"SAI-10k will use annotation defaults and the result will not be cached.", flush=True)
                    structures = {}
                    skip_cache = True
                    db_unavailable = True

        for transcript_scores, transcript_id_without_version in zip(candidate_transcripts, candidate_ids):
            transcript_structure = structures.get(transcript_id_without_version)
            if transcript_structure:
                transcript_scores.update(transcript_structure)
            elif not db_unavailable:
                # DB was reachable but this row is missing. Without skip_cache,
                # the degraded result (SAI-10k falling back to annotation
                # defaults) would be cached permanently and re-served forever.
                print(f"WARNING: transcript {transcript_id_without_version} not found in "
                      f"transcripts_hg{genome_version} for {variant}; "
                      f"SAI-10k will use annotation defaults and the result will not be cached.", flush=True)
                skip_cache = True
    db_enrich_ms = (time.perf_counter() - db_enrich_t0) * 1000

    # Single source of truth for canonical-transcript selection (priority, then
    # sum of |DS_*|). Used for both (a) which transcript's ALL_NON_ZERO_SCORES
    # to return to the client and (b) which transcript to feed into SAI-10k-calc.
    sai10k_t0 = time.perf_counter()
    selected_transcript = sai10k_select_transcript(candidate_transcripts)
    all_non_zero_scores = selected_transcript["ALL_NON_ZERO_SCORES"] if selected_transcript else None
    # Prefer STRAND (from the SpliceAI annotator, structurally guaranteed) and
    # fall back to t_strand from the external transcript-annotations JSON. This
    # matches sai10k_predictions.py:1150 so the strand reported in the JSON
    # response matches the strand used to compute the SAI-10k aberrations.
    all_non_zero_scores_strand = (selected_transcript.get("STRAND") or selected_transcript.get("t_strand")) if selected_transcript else None
    all_non_zero_scores_transcript_id = selected_transcript["t_id"] if selected_transcript else None

    # Compute SAI-10k-calc predictions for the selected transcript (before we
    # delete ALL_NON_ZERO_SCORES / STRAND from the scores dicts below).
    sai10k_predictions = None
    sai10k_predictions_error = None
    fasta_open_ms = 0.0
    if selected_transcript:
        fasta_t0 = time.perf_counter()
        sai10k_fasta = _get_sai10k_fasta(genome_version)
        fasta_open_ms = (time.perf_counter() - fasta_t0) * 1000
        # Premature-stop detection requires the FASTA. When it failed to open,
        # the resulting predictions have null stop_codon_introduced / aa_change
        # fields — don't cache that degraded response, otherwise it would
        # persist past FASTA recovery until SAI10K_VERSION is bumped.
        if sai10k_fasta is None:
            skip_cache = True
        try:
            sai10k_predictions = sai10k_get_transcript_predictions(
                selected_transcript, pos,
                chrom=chrom, ref=ref, alt=alt,
                fasta=sai10k_fasta,
            )
        except Exception as e:
            # Log the full exception server-side; return a generic message to the
            # client so internal details (file paths, transcript IDs, dict-key
            # names from KeyErrors, etc.) aren't echoed back through the JSON
            # response that index.html renders.
            print(f"WARNING: Error computing SAI-10k predictions for {variant}: {type(e).__name__}: {e}")
            traceback.print_exc()
            sai10k_predictions_error = "Internal error computing SAI-10k predictions."
            # Don't cache responses where SAI-10k bailed mid-request — the
            # exception may have been a transient DB / FASTA / parser hiccup,
            # and the client would otherwise see the error message forever.
            skip_cache = True

    sai10k_total_ms = (time.perf_counter() - sai10k_t0) * 1000

    # Strip the internal timing dict before serializing predictions to the
    # client; emit a single SAI10K_TIMING log line so processing times can be
    # derived from the Cloud Run logs (filter prefix: "SAI10K_TIMING ").
    inner_timing = sai10k_predictions.pop('_timing_ms', None) if sai10k_predictions else None
    n_exons = len(selected_transcript.get('EXON_STARTS', []) or []) if selected_transcript else 0
    selected_t_id = selected_transcript.get('t_id') if selected_transcript else None
    breakdown = ''
    if inner_timing:
        breakdown = (
            f" determine={inner_timing['determine']:.1f}ms"
            f" annotate={inner_timing['annotate']:.1f}ms"
            f" premature_stop={inner_timing['premature_stop']:.1f}ms"
            f" n_aberrations={inner_timing['n_aberrations']}"
            f" n_premature_stop_calls={inner_timing['n_premature_stop_calls']}"
        )
    print(
        f"SAI10K_TIMING variant={variant} hg{genome_version} "
        f"total={sai10k_total_ms:.1f}ms db_enrich={db_enrich_ms:.1f}ms "
        f"fasta_open={fasta_open_ms:.1f}ms{breakdown} "
        f"n_candidates={len(candidate_transcripts)} n_exons={n_exons} "
        f"selected_transcript={selected_t_id} error={bool(sai10k_predictions_error)}",
        flush=True,
    )

    for transcript_scores in candidate_transcripts:
        for redundant_key in ("ALLELE", "NAME", "STRAND", "ALL_NON_ZERO_SCORES"):
            transcript_scores.pop(redundant_key, None)

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
        "sai10kPredictions": sai10k_predictions,
        "sai10kPredictionsError": sai10k_predictions_error,
        # Internal sentinel: stripped by run_splice_prediction_tool before the
        # response is returned to the client. True when any of the following
        # produced a degraded response that should not be cached past the
        # underlying recovery:
        #   - per-request DB connection couldn't be acquired (transcript-
        #     structure enrichment skipped — see the get_db_connection block).
        #   - SAI-10k FASTA failed to open (premature-stop detection skipped).
        #   - SAI-10k computation raised an exception.
        "_skip_cache": skip_cache,
    }


def get_pangolin_scores(variant, genome_version, distance_param, mask_param, basic_or_comprehensive_param):
    if genome_version not in ("37", "38"):
        raise ValueError(f"Invalid genome_version: {genome_version}")

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

    features_db = gffutils.FeatureDB(PANGOLIN_ANNOTATION_PATHS[(genome_version, basic_or_comprehensive_param)])
    scores = process_variant_using_pangolin(
        0, chrom, int(pos), ref, alt, features_db, pangolin_models, PangolinArgs)

    if not scores:
        return {
            "variant": variant,
            "source": "pangolin",
            "error": f"Pangolin was unable to compute scores for this variant",
        }

    # Enrich each transcript with annotations, then select the best one for visualization
    # using the same priority-based logic as SpliceAI (MS > MP > C > N, tie-break on score sum).
    candidate_transcripts = []
    for transcript_scores in scores:
        if "ALL_NON_ZERO_SCORES" not in transcript_scores:
            continue

        transcript_id_without_version = transcript_scores.get("NAME", "").split(".")[0]

        transcript_annotations = SHARED_TRANSCRIPT_ANNOTATIONS[(genome_version, basic_or_comprehensive_param)].get(transcript_id_without_version)
        if transcript_annotations is None:
            raise ValueError(f"Missing annotations for {transcript_id_without_version} in {genome_version} annotations")

        transcript_scores.update(transcript_annotations)
        candidate_transcripts.append(transcript_scores)

    # Select transcript: highest priority, then highest sum of |DS_SL| + |DS_SG|
    selected_transcript = None
    best_priority = -1
    best_score_sum = -1.0
    for transcript_scores in candidate_transcripts:
        priority = TRANSCRIPT_PRIORITY_ORDER.get(transcript_scores.get('t_priority', 'N'), 0)
        score_sum = abs(float(transcript_scores.get('DS_SL', 0))) + abs(float(transcript_scores.get('DS_SG', 0)))
        if priority > best_priority or (priority == best_priority and score_sum > best_score_sum):
            selected_transcript = transcript_scores
            best_priority = priority
            best_score_sum = score_sum

    all_non_zero_scores = selected_transcript["ALL_NON_ZERO_SCORES"] if selected_transcript else None
    all_non_zero_scores_strand = selected_transcript["STRAND"] if selected_transcript else None
    all_non_zero_scores_transcript_id = selected_transcript["NAME"] if selected_transcript else None

    for transcript_scores in candidate_transcripts:
        for redundant_key in ("NAME", "STRAND", "ALL_NON_ZERO_SCORES"):
            transcript_scores.pop(redundant_key, None)

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
    return run_splice_prediction_tool(tool_name="spliceai")


@app.route("/pangolin/", methods=['POST', 'GET'])
def run_pangolin():
    return run_splice_prediction_tool(tool_name="pangolin")


_FIRST_REQUEST_LOGGED = False


def run_splice_prediction_tool(tool_name):
    """Handles API request for splice prediction.

    DB connections are taken from the pool only for short bursts (cache lookup,
    rate-limit check, log writes, cache writes, and per-call transcript-structure
    SELECTs inside get_spliceai_scores). The model inference runs without holding
    any pooled connection, so a slow inference can't starve the pool.

    Args:
        tool_name (str): "spliceai" or "pangolin"
    """

    global _FIRST_REQUEST_LOGGED
    if not _FIRST_REQUEST_LOGGED:
        _FIRST_REQUEST_LOGGED = True
        print(f"[startup pid={os.getpid()}] first request received "
              f"+{time.time() - _PROCESS_START_TIME:.2f}s after process start", flush=True)

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
    # Type-check before .strip() — a non-string payload (e.g. {"variant": 123}
    # from request.get_json) would otherwise raise AttributeError and 500
    # instead of producing a clean 400.
    if not isinstance(variant, str):
        return error_response(f'"variant" value must be a string rather than a {type(variant)}.\n', source=tool_name)

    variant = variant.strip().strip("'").strip('"').strip(",")
    if not variant:
        return error_response(f'"variant" not specified.\n', source=tool_name)

    genome_version = params.get("hg")
    if not genome_version:
        return error_response(f'"hg" not specified. The URL must include an "hg" arg: hg=37 or hg=38. For example: {example_url}\n', source=tool_name)

    if genome_version not in ("37", "38"):
        return error_response(f'Invalid "hg" value: "{genome_version}". The value must be either "37" or "38". For example: {example_url}\n', source=tool_name)

    if genome_version != GENOME_VERSION:
        return error_response(f'This service only handles hg{GENOME_VERSION} requests, but received hg={genome_version}. Route hg{genome_version} requests to the matching per-genome service.\n', source=tool_name)

    distance_param = params.get("distance", DEFAULT_DISTANCE)
    try:
        distance_param = int(distance_param)
    except Exception as e:
        return error_response(f'Invalid "distance": "{distance_param}". The value must be an integer.\n', source=tool_name)

    if distance_param < 0:
        return error_response(f'Invalid "distance": "{distance_param}". The value must be non-negative.\n', source=tool_name)

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

    # check cache before processing the variant (short DB scope)
    results = {}
    if not force:
        with get_db_connection() as conn:
            results = get_splicing_scores_from_cache(conn, tool_name, variant, genome_version, distance_param, mask_param, basic_or_comprehensive_param)

    if results:
        # Cache hit: brief log scope, then fall through to response building.
        with get_db_connection() as conn:
            log(conn, f"{tool_name}:from-cache", ip=user_ip, variant=variant, genome=genome_version, distance=distance_param, mask=mask_param, bc=basic_or_comprehensive_param, variant_consequence=variant_consequence)
            if "error" in results:
                log(conn, f"{tool_name}:error", ip=user_ip, variant=variant, genome=genome_version, distance=distance_param, mask=mask_param, details=results["error"], bc=basic_or_comprehensive_param, variant_consequence=variant_consequence)
    else:
        # Rate-limit check (short DB scope).
        with get_db_connection() as conn:
            error_message = exceeds_rate_limit(conn, user_ip, params)
        if error_message:
            print(f"{logging_prefix}: {user_ip}: response: {error_message}", flush=True)
            return error_response(error_message, source=tool_name, status=429)

        # Model inference runs without a pooled DB connection. get_spliceai_scores
        # acquires its own short-lived connection internally for transcript-structure
        # SELECTs after inference completes; get_pangolin_scores does no DB work.
        try:
            if tool_name == "spliceai":
                results = get_spliceai_scores(variant, genome_version, distance_param, int(mask_param), basic_or_comprehensive_param)
            elif tool_name == "pangolin":
                pangolin_mask_param = "True" if mask_param == "1" else "False"
                results = get_pangolin_scores(variant, genome_version, distance_param, pangolin_mask_param, basic_or_comprehensive_param)
            else:
                raise ValueError(f"Invalid tool_name: {tool_name}")
        except Exception as e:
            # Internal exceptions can carry implementation detail (file paths,
            # KeyError on internal dict shapes, missing-annotation messages
            # naming specific transcript IDs). Don't echo them to clients —
            # log server-side and return a generic message. User-input
            # validation errors above are passed straight to error_response and
            # are unaffected by this wrapping.
            print(f"{logging_prefix}: 500 in {tool_name} for variant {variant}: {type(e).__name__}: {e}", flush=True)
            traceback.print_exc()
            return error_response(
                "Internal server error while computing predictions. "
                "If this persists, please file an issue at "
                "https://github.com/broadinstitute/SpliceAI-lookup/issues.",
                source=tool_name,
                status=500,
            )

        # Strip the internal sentinel before anything downstream sees `results`.
        # Set by get_spliceai_scores when the per-request DB connection couldn't
        # be acquired, so the transcript-structure enrichment was skipped and
        # the result reflects degraded inputs — don't cache it.
        skip_cache = results.pop("_skip_cache", False)

        # Post-inference: log + cache write + (if error) error log, all in one short DB scope.
        duration = (datetime.now() - start_time).total_seconds()
        with get_db_connection() as conn:
            log(conn, f"{tool_name}:computed", ip=user_ip, duration=duration, variant=variant, genome=genome_version, distance=distance_param, mask=mask_param, bc=basic_or_comprehensive_param, variant_consequence=variant_consequence)
            if "error" not in results and not skip_cache:
                add_splicing_scores_to_cache(conn, tool_name, variant, genome_version, distance_param, mask_param, basic_or_comprehensive_param, results)
            elif "error" in results:
                log(conn, f"{tool_name}:error", ip=user_ip, variant=variant, genome=genome_version, distance=distance_param, mask=mask_param, details=results["error"], bc=basic_or_comprehensive_param, variant_consequence=variant_consequence)

    # Echo only a whitelist of input params back to the client. A prior version
    # used `response_json.update(params)`, which reflected every query-string
    # field unchanged — combined with the `.html()` sinks in index.html and
    # the URL-hash auto-submit on page load, that let a crafted link execute
    # arbitrary HTML in visitors' browsers.
    ECHO_PARAM_KEYS = (
        "variant", "hg", "bc", "distance", "mask", "raw", "variant_consequence",
    )
    response_json = {k: params[k] for k in ECHO_PARAM_KEYS if k in params}
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
    # On Cloud Run the X-Forwarded-For header is "<client-supplied>..., <verified-client>",
    # where the final entry is appended by GCP's load balancer and is the only value
    # the client cannot forge. Using the whole header (or the first entry) lets an
    # attacker spoof a different IP to bypass per-IP rate limits or frame a victim
    # IP into the 1-week restricted_ips ban list.
    xff = request.environ.get("HTTP_X_FORWARDED_FOR", "")
    if not xff:
        return None
    return xff.rsplit(",", 1)[-1].strip() or None


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
    # Serve as text/plain so a `path` containing HTML/JS can't be rendered as
    # markup by browsers — Flask's default mimetype is text/html, which would
    # make this a reflected-XSS sink. Talisman's default CSP blocks inline
    # scripts but does not strip the HTML content type itself.
    return Response(
        f"SpliceAI-lookup APIs: invalid endpoint {path}",
        status=404,
        mimetype='text/plain',
    )


print(f"[startup pid={os.getpid()}] server.py module loaded in "
      f"{time.time() - _PROCESS_START_TIME:.2f}s (tool={TOOL}, genome={GENOME_VERSION})", flush=True)


if '__main__' == __name__ or os.environ.get('RUNNING_ON_GOOGLE_CLOUD_RUN'):
    app.run(debug=DEBUG, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))
