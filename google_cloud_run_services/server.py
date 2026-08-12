import time
_PROCESS_START_TIME = time.time()
from datetime import datetime
import gc
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

# Matches a bare chrom+pos with no ref/alt (e.g. "chr8-140300615", "chr8:140300615"),
# used for REF-only score requests where no ALT allele is available. fullmatch()
# against this is tried only after VARIANT_RE's fullmatch fails, so a full
# variant string is never misparsed as a position (the trailing "-C-G" etc.
# keeps it from matching this pattern at all).
POSITION_RE = re.compile(
    r"(chr)?(?P<chrom>[0-9XYMTt]{1,2})"
    r"[-\s:]+"
    r"(?P<pos>[0-9]{1,9})"
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

# The single Gencode gene set this service answers for, pinned the same way GENOME_VERSION is.
# Each container therefore loads exactly one annotator instead of accumulating one per gene set
# its workers happened to be asked for, which is what pushed instances past their memory limit.
# Defaults to "basic" so the long-standing service URLs keep their existing behaviour; the
# comprehensive services set GENE_SET=comprehensive explicitly at deploy time.
GENE_SET = os.environ.get("GENE_SET", "basic")
if GENE_SET not in ("basic", "comprehensive"):
    raise ValueError(f'Environment variable "GENE_SET" should be set to either "basic" or "comprehensive" instead of: "{GENE_SET}"')

# "prod" or "dev". Mixed into the cache key (see get_splicing_scores_cache_key) so a dev
# revision experimenting against the shared database can neither read nor overwrite the entries
# production serves. Defaults to "prod" and contributes nothing to the key in that case, so
# every existing production cache entry stays valid.
DEPLOYMENT = os.environ.get("DEPLOYMENT", "prod")

if TOOL == "spliceai":
    from spliceai.utils import Annotator, get_delta_scores, get_reference_scores, MIN_SCORE_THRESHOLD

    # The per-position score fields of an ALL_NON_ZERO_SCORES row from a delta-score response,
    # and the threshold a row must clear to count as reportable. Both are tool-specific;
    # MIN_SCORE_THRESHOLD is imported from the model package so it can't drift from the value the
    # model itself filtered on. REF-only responses carry a different, smaller set of row fields,
    # which is why count_scores_above_threshold is only called on the delta-score path.
    PER_POSITION_SCORE_FIELDS = ("RA", "AA", "RD", "AD")

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
    from pangolin.pangolin import process_position as process_position_using_pangolin
    from pangolin.pangolin import MIN_SCORE_THRESHOLD

    # see the comment on the spliceai branch above
    PER_POSITION_SCORE_FIELDS = ("SL_REF", "SL_ALT", "SG_REF", "SG_ALT")
    from pangolin.model import torch, Pangolin, L, W, AR
    import gffutils

    PANGOLIN_ANNOTATION_PATHS = {
        ("37", "basic"): f"/gencode.{GENCODE_VERSION}lift37.basic.annotation.without_chr_prefix.db",
        ("38", "basic"): f"/gencode.{GENCODE_VERSION}.basic.annotation.db",
        ("37", "comprehensive"): f"/gencode.{GENCODE_VERSION}lift37.annotation.without_chr_prefix.db",
        ("38", "comprehensive"): f"/gencode.{GENCODE_VERSION}.annotation.db",
    }

    # The 12 Pangolin models (4 splice-score types x 3 replicates each), populated by
    # init_pangolin(). Their weights ship with the pangolin package and depend on neither the
    # genome version nor the basic/comprehensive gene set, so one flat list serves every request.
    #
    # DANGER, and the reason init_pangolin() builds a local list and publishes it in a single
    # assignment: pangolin.compute_score() slices this list positionally, `for model in
    # models[3*j: 3*j+3]` with j in range(4). Python slicing returns a short (or empty) slice
    # instead of raising, so a list observed while it is still being filled yields *silently
    # wrong* scores. That is exactly the April-November 2024 bug, when the cache was a dict of
    # lists appended to in place behind a `if not PANGOLIN_MODELS[key]` guard: a second
    # concurrent request saw the guard satisfied after the first append and scored against a
    # 1-of-12-model list. Never append to the published list; rebuild and reassign.
    PANGOLIN_MODELS = None
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


def init_pangolin():
    """Populate the module-level PANGOLIN_MODELS cache with all 12 Pangolin models.

    Idempotent, and safe to call from concurrent requests: the models are built into a local
    list and published in one assignment, so PANGOLIN_MODELS is only ever None or a complete
    12-model list. Two callers racing here both build a full list and one overwrites the other
    -- wasteful but correct. See the PANGOLIN_MODELS comment for why a partially published
    list would silently corrupt scores.
    """
    global PANGOLIN_MODELS
    if PANGOLIN_MODELS is not None:
        return

    t0 = time.time()
    models = []
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
            models.append(model)

    PANGOLIN_MODELS = models
    print(f"[startup pid={os.getpid()}] init_pangolin() loaded {len(models)} models in {time.time() - t0:.2f}s "
          f"(+{time.time() - _PROCESS_START_TIME:.2f}s since start)", flush=True)


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


def parse_position(position_str):
    match = POSITION_RE.fullmatch(position_str)
    if not match:
        raise ValueError(f"Unable to parse position: {position_str}")

    return match['chrom'], int(match['pos'])


def _env_flag(name, default=False):
    """Parse a boolean environment variable.

    Treats only 1/true/yes/on (case-insensitive) as true and 0/false/no/off/""
    as false, so e.g. NAME=0 reads as false instead of the way
    bool(os.environ.get(name)) would make any non-empty string truthy. Returns
    `default` when the variable is unset.
    """
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() in ("1", "true", "yes", "on")


# Connection parameters. Each defaults to the production Cloud Run / Cloud SQL
# value so the deployed service is unchanged, but every field can be overridden
# via env vars to point a local instance at your own PostgreSQL (e.g.
# DB_HOST=localhost DB_PORT=5432). A DB_HOST starting with "/" is treated by
# psycopg2 as a Unix-socket directory (the Cloud SQL default); a regular
# hostname uses TCP.
DB_CONNECT_KWARGS = dict(
    dbname=os.environ.get("DB_NAME", "spliceai-lookup-db"),
    user=os.environ.get("DB_USER", "postgres"),
    password=os.environ.get("DB_PASSWORD"),
    host=os.environ.get("DB_HOST", "/cloudsql/spliceai-lookup-412920:us-central1:spliceai-lookup-db"),
    port=os.environ.get("DB_PORT", "5432"),
    connect_timeout=5,
)

# Whether to use a database at all. DB-backed features (response caching, per-IP
# rate limiting, SAI-10k transcript-structure enrichment) are only active when
# this is true; otherwise get_db_connection yields None and every caller
# degrades gracefully. Defaults to on whenever a DB_PASSWORD is present (Cloud
# Run injects it as a secret) and off otherwise (a local `docker run` with no
# env), preserving prior behavior with no deploy change. Set DATABASE_ENABLED=1
# explicitly to attach a local PostgreSQL that uses passwordless (trust) auth,
# or DATABASE_ENABLED=0 to force the no-DB path even when a password is set.
DATABASE_ENABLED = _env_flag("DATABASE_ENABLED", default=bool(os.environ.get("DB_PASSWORD")))

# Explicitly disable per-IP rate limiting (independent of the database), e.g. for
# a private instance that DOES have a database but should not throttle. Rate
# limiting is also off whenever DATABASE_ENABLED is False, since it is entirely
# DB-backed.
DISABLE_RATE_LIMIT = _env_flag("DISABLE_RATE_LIMIT")

# Comma-separated list of blocked IP addresses, rejected at the door (see
# block_ips) before any routing, DB query, or model inference.
BLOCKED_IPS = frozenset(ip.strip() for ip in os.environ.get("BLOCKED_IPS", "").split(",") if ip.strip())
if not BLOCKED_IPS:
    print("WARNING: BLOCKED_IPS env var is unset/empty; no IPs will be blocked at the door", flush=True)

# Module-level connection pool for Cloud SQL. Flask under Cloud Run typically
# serves multiple concurrent requests per instance via threaded workers, so use
# ThreadedConnectionPool (thread-safe) rather than SimpleConnectionPool.
# Each gunicorn worker is a separate process (forked under --preload) with its
# own pool and runs single-threaded (--threads 1), so it serves one request at a
# time and needs only one DB connection. maxconn=2 keeps a one-slot safety margin
# while bounding the total: workers x instances x services x maxconn must stay
# under the server's max_connections (75). The previous maxconn=80 assumed Cloud
# Run's default concurrency of 80, but the deploy pins --workers/--concurrency to
# 6, so 80 per worker let a single instance open far more connections than the
# tier allows -- exhausting max_connections and leaving dozens of idle backends.
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

# DDL for the small operational tables the server reads/writes for caching and
# rate limiting. All idempotent (CREATE ... IF NOT EXISTS) so running them
# against a database where they already exist (e.g. production Cloud SQL) is a
# no-op. The large transcripts_hg37/38 enrichment tables are intentionally NOT
# created here — they require a separate bulk data load from genePred files (see
# build_and_deploy.py's update_transcript_tables command), and SAI-10k degrades
# gracefully to the bundled annotations when they are absent.
_SCHEMA_DDL_STATEMENTS = (
    "CREATE TABLE IF NOT EXISTS cache (key TEXT UNIQUE, value TEXT, counter INT, accessed TIMESTAMP DEFAULT now())",
    "CREATE INDEX IF NOT EXISTS cache_index ON cache (key)",
    "CREATE TABLE IF NOT EXISTS log (event_name TEXT, ip TEXT, logtime TIMESTAMP DEFAULT now(), duration REAL, variant TEXT, genome VARCHAR(10), bc VARCHAR(20), distance INT, mask INT4, details TEXT, variant_consequence TEXT)",
    "CREATE INDEX IF NOT EXISTS idx_log_ip_logtime ON log USING btree (ip, logtime DESC)",
    "CREATE INDEX IF NOT EXISTS idx_log_event_name ON log USING btree (event_name)",
    "CREATE TABLE IF NOT EXISTS restricted_ips (ip TEXT UNIQUE, created TIMESTAMP DEFAULT now())",
    "CREATE INDEX IF NOT EXISTS idx_restricted_ips_created ON restricted_ips USING btree (created)",
    "CREATE TABLE IF NOT EXISTS whitelist_ips (ip TEXT UNIQUE, created TIMESTAMP DEFAULT now())",
)
_database_schema_init_attempted = False
_database_schema_lock = threading.Lock()


def _init_database_schema(conn):
    """Create the cache / log / restricted_ips / whitelist_ips tables if missing.

    Runs at most once per process (guarded by _database_schema_init_attempted)
    the first time a usable connection is obtained, so a freshly-pointed local
    PostgreSQL works without any manual bootstrap. Idempotent and safe on
    production where the tables already exist. The attempt flag is set even when
    the DDL fails (e.g. the DB user lacks CREATE rights) so a permanent failure
    isn't retried — re-running the 8 statements under the lock on every request
    would serialize traffic. The failure is logged and non-fatal: the server
    still runs via the fail-open paths.
    """
    global _database_schema_init_attempted
    if _database_schema_init_attempted or conn is None:
        return
    with _database_schema_lock:
        if _database_schema_init_attempted:
            return
        try:
            with conn.cursor() as cursor:
                for ddl in _SCHEMA_DDL_STATEMENTS:
                    cursor.execute(ddl)
            conn.commit()
            print("Ensured DB schema exists (cache, log, restricted_ips, whitelist_ips)", flush=True)
        except Exception as e:
            try:
                conn.rollback()
            except Exception:
                pass
            print(f"WARNING: DB schema init failed (continuing without it): {e}", flush=True)
        finally:
            _database_schema_init_attempted = True


def _try_init_database_pool():
    """Attempt to (re)initialise DATABASE_CONNECTION_POOL, throttled to once per _DATABASE_POOL_INIT_RETRY_SECONDS."""
    global DATABASE_CONNECTION_POOL, _database_pool_init_last_attempt
    if not DATABASE_ENABLED:
        return
    with _database_pool_init_lock:
        if DATABASE_CONNECTION_POOL is not None:
            return
        now = time.monotonic()
        if _database_pool_init_last_attempt and now - _database_pool_init_last_attempt < _DATABASE_POOL_INIT_RETRY_SECONDS:
            return
        _database_pool_init_last_attempt = now
        try:
            DATABASE_CONNECTION_POOL = ThreadedConnectionPool(minconn=1, maxconn=2, **DB_CONNECT_KWARGS)
            print("Successfully initialised DB connection pool", flush=True)
        except Exception as e:
            print(f"WARNING: DB connection pool init failed; falling back to per-request connections: {e}", flush=True)


# The pool is initialised lazily on the first request (see get_db_connection),
# not at import time. Under gunicorn --preload the module is imported once in the
# arbiter before workers fork; initialising the pool here would share a single
# connection's socket across all forked workers and corrupt the protocol. Lazy
# init means each worker opens its own pool after the fork.


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
    if not DATABASE_ENABLED:
        # No database configured (e.g. a local instance) — yield None so callers
        # skip caching, rate limiting, and DB lookups without attempting (and
        # logging a failed) connection on every request.
        yield None
        return

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

    # Ensure the operational tables exist the first time we get a usable
    # connection (no-op after the first success, and on production where they
    # already exist).
    if conn is not None:
        _init_database_schema(conn)

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


# Operational tables (cache, log, restricted_ips, whitelist_ips) are created
# automatically by _init_database_schema (see _SCHEMA_DDL_STATEMENTS above).

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
        if DISABLE_RATE_LIMIT or conn is None:
            return False

        if is_user_on_whitelist(conn, params.get("ip")):
            return False

        # check if the user has exceeded the rate limit or is on the list of restricted IPs
        rows = run_sql(conn, "SELECT COUNT(ip) FROM restricted_ips WHERE ip=%s AND created >= NOW() - INTERVAL '1 weeks'", (user_ip,))
        is_user_currently_blocked = rows and int(rows[0][0]) > 0
        if is_user_currently_blocked:
            return RATE_LIMIT_ERROR_MESSAGE

        rows = run_sql(conn, "SELECT COUNT(ip) FROM log WHERE event_name LIKE %s AND ip=%s AND logtime >= NOW() - INTERVAL '7 minutes'", ("%computed%", user_ip))
        did_user_exceed_rate_limit = rows and int(rows[0][0]) >= 150
        if did_user_exceed_rate_limit and not is_user_on_whitelist(conn, user_ip):
            # the user has exceeded the rate limit: computing scores for 150 or more variants in the last 7 minutes
            rows = run_sql(conn, "SELECT COUNT(ip) FROM log WHERE event_name='rate_limit_exceeded' AND ip=%s AND logtime >= NOW() - INTERVAL '5 minutes'", (user_ip,))
            user_hit_rate_limit_exceeded_recently = rows and int(rows[0][0]) > 0
            if not user_hit_rate_limit_exceeded_recently:
                # the user will receive at most one "rate_limit_exceeded" event every 5 minutes
                log(conn, f"rate_limit_exceeded", ip=user_ip)
                rows = run_sql(conn, "SELECT COUNT(ip) FROM log WHERE event_name='rate_limit_exceeded' AND ip=%s AND logtime >= NOW() - INTERVAL '1 days'", (user_ip,))
                user_triggered_too_many_rate_limit_exceeded_errors_today = rows and int(rows[0][0]) >= 15
                if user_triggered_too_many_rate_limit_exceeded_errors_today:
                    # the user has hit the limit of 15 or more "rate_limit_exceeded" events during the last 24 hours
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
# invalidated and recomputed. v22 discards the v21 entries, whose protein windows are
# null for variants whose footprint straddles a shifted splice boundary (issue #134) --
# those entries predate the substitution-clipping fix and would otherwise keep serving
# a response with no protein sequence for exactly the variants the fix targets.
SAI10K_VERSION = "v22"

# Bump CACHE_VERSION whenever the shape of the cached response changes for either tool.
# v2 added the per-position rows (ALL_NON_ZERO_SCORES) and nNonZeroScores to the cached copy
# so the /scores endpoints can serve any transcript without re-running the model. v3 discards
# the v2 entries, which were written by a revision whose model packages predated the per-position
# ref/alt bases, and so hold rows the current code cannot render.
CACHE_VERSION = "v3"


def is_score_above_threshold(row):
    """Whether one ALL_NON_ZERO_SCORES row cleared the model's reporting threshold.

    The model packages also report the variant's own position and the delta-score maxima, whose
    scores can be below their threshold. Judging each row by its own scores keeps that
    distinction here, where the consumers live, instead of asking the model packages to mark
    rows on their behalf.

    The scores are strings the model already rounded to 3 decimals, so a true score just under
    the threshold (0.0095 renders as "0.010") still passes. That residual is the safe direction
    for the count below, and for the visualization it is a rare single-position overshoot.

    Args:
        row (dict): one ALL_NON_ZERO_SCORES entry from a delta-score response. REF-only
            responses carry a different, smaller set of row fields and are not handled here.

    Returns:
        bool
    """
    return max(float(row[score_field]) for score_field in PER_POSITION_SCORE_FIELDS) >= MIN_SCORE_THRESHOLD


def count_scores_above_threshold(transcript_scores):
    """Count the per-position rows that cleared the model's reporting threshold.

    Zero is what tells the UI to disable that transcript's table icon, so the positions the
    model adds below its threshold must not be counted here.

    Args:
        transcript_scores (dict): one transcript's entry, carrying ALL_NON_ZERO_SCORES.

    Returns:
        int
    """
    return sum(1 for row in (transcript_scores.get("ALL_NON_ZERO_SCORES") or []) if is_score_above_threshold(row))


def get_splicing_scores_cache_key(tool_name, variant, genome_version, distance, mask, basic_or_comprehensive="basic", is_position_only=False):
    # REF-only (position-only) spliceai responses never run SAI-10k-calc (it requires
    # an ALT allele), so bumping SAI10K_VERSION should not invalidate their cache entries.
    suffix = f"__sai10k-{SAI10K_VERSION}" if tool_name == "spliceai" and not is_position_only else ""
    # Dev revisions share production's database, so without this a dev revision computing a
    # result would write it onto the key production reads -- publishing unreviewed output to
    # every user. Production adds nothing to the key, keeping existing entries valid.
    if DEPLOYMENT != "prod":
        suffix += f"__{DEPLOYMENT}"
    return f"{tool_name}__{variant}__hg{genome_version}__d{distance}__m{mask}__{basic_or_comprehensive}__{CACHE_VERSION}{suffix}"


def get_splicing_scores_from_cache(conn, tool_name, variant, genome_version, distance, mask, basic_or_comprehensive="basic", is_position_only=False):
    results = {}
    key = get_splicing_scores_cache_key(tool_name, variant, genome_version, distance, mask, basic_or_comprehensive, is_position_only)
    try:
        rows = run_sql(conn, f"SELECT value FROM cache WHERE key=%s", (key,))
        if rows:
            results = json.loads(rows[0][0])
            results["source"] += ":cache"
    except Exception as e:
        print(f"Cache error: {e}", flush=True)

    return results


def add_splicing_scores_to_cache(conn, tool_name, variant, genome_version, distance, mask, basic_or_comprehensive, results, is_position_only=False):
    key = get_splicing_scores_cache_key(tool_name, variant, genome_version, distance, mask, basic_or_comprehensive, is_position_only)
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
    # The visualization tracks re-filter these rows with their own cutoffs, which are looser
    # than the model's, so hand them only the rows that cleared the model threshold. Otherwise
    # the positions the model adds below threshold would be drawn as marks that were never
    # drawn before. The /scores endpoints serve the full set from the cached copy.
    all_non_zero_scores = [
        row for row in selected_transcript["ALL_NON_ZERO_SCORES"] if is_score_above_threshold(row)
    ] if selected_transcript else None
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

    # ALL_NON_ZERO_SCORES stays in the dict so it reaches the cache, which lets the
    # /spliceai/scores endpoint serve any transcript's per-position rows without re-running
    # the model. run_splice_prediction_tool drops it from the HTTP response and keeps only
    # nNonZeroScores, which is all the results table needs in order to decide whether the
    # per-position table is worth offering for a given transcript.
    for transcript_scores in candidate_transcripts:
        transcript_scores["nNonZeroScores"] = count_scores_above_threshold(transcript_scores)
        for redundant_key in ("ALLELE", "NAME", "STRAND"):
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


def get_spliceai_reference_scores(variant, genome_version, distance_param, basic_or_comprehensive_param):
    """REF-only counterpart of get_spliceai_scores: variant is a bare chrom-pos
    position (no ref/alt). No SAI-10k-calc predictions or DB transcript-
    structure enrichment -- those require an ALT allele.
    """
    chrom, pos = parse_position(variant)

    # spliceai's normalise_chrom() handles "chr" prefix mismatches but not the
    # M↔MT alias -- see the matching comment in get_spliceai_scores.
    if chrom.upper() in {"M", "MT"} and genome_version in MITO_CHROM_NAME:
        chrom = MITO_CHROM_NAME[genome_version]

    try:
        scores = get_reference_scores(
            chrom, pos,
            SPLICEAI_ANNOTATOR[(genome_version, basic_or_comprehensive_param)],
            distance_param)
    except Exception as e:
        print(f"ERROR while computing SpliceAI REF scores for {variant}: {e}")
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
            "error": f"The SpliceAI model did not return any scores for {variant}. This may be because the position does "
                     f"not overlap any exons or introns defined by the GENCODE '{basic_or_comprehensive_param}' annotation.",
        }

    # Enrich each transcript_scores with in-memory annotations (no DB), same as get_spliceai_scores.
    candidate_transcripts = []
    for transcript_scores in scores:
        transcript_id_without_version = transcript_scores.get("NAME", "").split(".")[0]
        transcript_annotations = SHARED_TRANSCRIPT_ANNOTATIONS[(genome_version, basic_or_comprehensive_param)].get(transcript_id_without_version)
        if transcript_annotations is None:
            raise ValueError(f"Missing annotations for {transcript_id_without_version} in {genome_version} annotations")
        transcript_scores.update(transcript_annotations)
        candidate_transcripts.append(transcript_scores)

    # Select the transcript to visualize: highest priority (MS > MP > C > N),
    # then highest combined |RA_MAX| + |RD_MAX|. Mirrors the selection logic in
    # get_pangolin_scores below.
    selected_transcript = None
    best_priority = -1
    best_score_sum = -1.0
    for transcript_scores in candidate_transcripts:
        priority = TRANSCRIPT_PRIORITY_ORDER.get(transcript_scores.get('t_priority', 'N'), 0)
        score_sum = abs(float(transcript_scores['RA_MAX'])) + abs(float(transcript_scores['RD_MAX']))
        if priority > best_priority or (priority == best_priority and score_sum > best_score_sum):
            selected_transcript = transcript_scores
            best_priority = priority
            best_score_sum = score_sum

    all_non_zero_scores = selected_transcript["ALL_NON_ZERO_SCORES"] if selected_transcript else None
    all_non_zero_scores_strand = (selected_transcript.get("STRAND") or selected_transcript.get("t_strand")) if selected_transcript else None
    all_non_zero_scores_transcript_id = selected_transcript["t_id"] if selected_transcript else None

    for transcript_scores in candidate_transcripts:
        for redundant_key in ("NAME", "STRAND", "ALL_NON_ZERO_SCORES"):
            transcript_scores.pop(redundant_key, None)

    return {
        "variant": variant,
        "genomeVersion": genome_version,
        "chrom": chrom,
        "pos": pos,
        "distance": distance_param,
        "scores": scores,
        "source": "spliceai:model",
        "isPositionOnly": True,
        "allNonZeroScores": all_non_zero_scores,
        "allNonZeroScoresStrand": all_non_zero_scores_strand,
        "allNonZeroScoresTranscriptId": all_non_zero_scores_transcript_id,
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

    init_pangolin()

    # FeatureDB stays per-request: it wraps a sqlite3 connection, which is bound to the thread
    # that opened it, and opening one is a few milliseconds against an on-disk index (unlike the
    # models, it does not read the database into memory). Caching it would trade nothing for a
    # cross-thread sharing hazard.
    features_db = gffutils.FeatureDB(PANGOLIN_ANNOTATION_PATHS[(genome_version, basic_or_comprehensive_param)])
    scores = process_variant_using_pangolin(
        0, chrom, int(pos), ref, alt, features_db, PANGOLIN_MODELS, PangolinArgs)

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

    # Brief DB scope: add EXON_STARTS/EXON_ENDS so the per-position table can label annotated
    # acceptors and donors for Pangolin the same way it does for SpliceAI. As on the SpliceAI
    # path, a failed lookup means the response is degraded (here: a blank Notes column), so
    # don't cache it -- otherwise it would be re-served long after the DB recovered.
    skip_cache = False
    if candidate_transcripts:
        candidate_ids = [transcript_scores.get("NAME", "").split(".")[0] for transcript_scores in candidate_transcripts]
        with get_db_connection() as conn:
            structures = get_transcript_structures(conn, candidate_ids, genome_version) if conn is not None else None

        if structures is None:
            # The DB was unreachable or the query failed. That is usually transient, so don't
            # cache a response whose Notes column would stay blank long after it recovered.
            print(f"WARNING: transcript-structure lookup unavailable for {variant}; the Pangolin "
                  f"per-position table will have no exon annotations and the result will not be cached.", flush=True)
            skip_cache = True
        else:
            for transcript_scores, transcript_id_without_version in zip(candidate_transcripts, candidate_ids):
                transcript_structure = structures.get(transcript_id_without_version)
                if transcript_structure:
                    # Copy only the exon coordinates. get_transcript_structures also returns
                    # CDS_START/CDS_END/STRAND, and merging those wholesale would replace
                    # Pangolin's own STRAND, which is read below into allNonZeroScoresStrand.
                    for exon_key in ("EXON_STARTS", "EXON_ENDS"):
                        transcript_scores[exon_key] = transcript_structure[exon_key]
                else:
                    # The query worked, this transcript simply has no row. That won't change on a
                    # retry, so cache the response anyway rather than recomputing it forever --
                    # only the Notes column is affected, and CACHE_VERSION can flush it later.
                    print(f"WARNING: transcript {transcript_id_without_version} not found in "
                          f"transcripts_hg{genome_version} for {variant}; the Pangolin per-position "
                          f"table will have no exon annotations for it.", flush=True)

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

    # The visualization tracks re-filter these rows with their own cutoffs, which are looser
    # than the model's, so hand them only the rows that cleared the model threshold. Otherwise
    # the positions the model adds below threshold would be drawn as marks that were never
    # drawn before. The /scores endpoints serve the full set from the cached copy.
    all_non_zero_scores = [
        row for row in selected_transcript["ALL_NON_ZERO_SCORES"] if is_score_above_threshold(row)
    ] if selected_transcript else None
    all_non_zero_scores_strand = selected_transcript["STRAND"] if selected_transcript else None
    all_non_zero_scores_transcript_id = selected_transcript["NAME"] if selected_transcript else None

    # see the matching comment in get_spliceai_scores: the per-position rows stay in the dict
    # so they reach the cache for /pangolin/scores, and are stripped from the HTTP response
    for transcript_scores in candidate_transcripts:
        transcript_scores["nNonZeroScores"] = count_scores_above_threshold(transcript_scores)
        for redundant_key in ("NAME", "STRAND"):
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
        # Internal sentinel, stripped by run_splice_prediction_tool before the response is
        # returned. True when the transcript-structure lookup failed, so the exon annotations
        # are missing and the response should not be cached past the underlying recovery.
        "_skip_cache": skip_cache,
    }


def get_pangolin_reference_scores(variant, genome_version, distance_param, basic_or_comprehensive_param):
    """REF-only counterpart of get_pangolin_scores: variant is a bare chrom-pos
    position (no ref/alt).
    """
    if genome_version not in ("37", "38"):
        raise ValueError(f"Invalid genome_version: {genome_version}")

    if basic_or_comprehensive_param not in ("basic", "comprehensive"):
        raise ValueError(f"Invalid basic_or_comprehensive_param: {basic_or_comprehensive_param}")

    chrom, pos = parse_position(variant)

    class PangolinArgs:
        reference_file = FASTA_PATH[genome_version]
        distance = distance_param

    init_pangolin()

    # Per-request FeatureDB: see the note in get_pangolin_scores.
    features_db = gffutils.FeatureDB(PANGOLIN_ANNOTATION_PATHS[(genome_version, basic_or_comprehensive_param)])
    scores = process_position_using_pangolin(0, chrom, int(pos), features_db, PANGOLIN_MODELS, PangolinArgs)

    if not scores:
        return {
            "variant": variant,
            "source": "pangolin",
            "error": f"Pangolin was unable to compute scores for this position",
        }

    candidate_transcripts = []
    for transcript_scores in scores:
        transcript_id_without_version = transcript_scores.get("NAME", "").split(".")[0]

        transcript_annotations = SHARED_TRANSCRIPT_ANNOTATIONS[(genome_version, basic_or_comprehensive_param)].get(transcript_id_without_version)
        if transcript_annotations is None:
            raise ValueError(f"Missing annotations for {transcript_id_without_version} in {genome_version} annotations")

        transcript_scores.update(transcript_annotations)
        candidate_transcripts.append(transcript_scores)

    # Select transcript: highest priority, then highest |S_REF|
    selected_transcript = None
    best_priority = -1
    best_score = -1.0
    for transcript_scores in candidate_transcripts:
        priority = TRANSCRIPT_PRIORITY_ORDER.get(transcript_scores.get('t_priority', 'N'), 0)
        score = abs(float(transcript_scores['S_REF']))
        if priority > best_priority or (priority == best_priority and score > best_score):
            selected_transcript = transcript_scores
            best_priority = priority
            best_score = score

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
        "distance": distance_param,
        "scores": scores,
        "source": "pangolin:model",
        "isPositionOnly": True,
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


@app.route("/spliceai/scores/", strict_slashes=False, methods=['POST', 'GET'])
def run_spliceai_scores():
    return run_splice_prediction_tool(tool_name="spliceai", scores_for_one_transcript=True)


@app.route("/pangolin/scores/", strict_slashes=False, methods=['POST', 'GET'])
def run_pangolin_scores():
    return run_splice_prediction_tool(tool_name="pangolin", scores_for_one_transcript=True)


def per_transcript_scores_response(results, transcript_id, tool_name):
    """Build the /scores response: the per-position rows for a single transcript.

    The rows are read out of the same cached result the main endpoint produced, so opening the
    per-position table for a transcript never re-runs the model.

    Args:
        results (dict): a full result dict, from the cache or freshly computed.
        transcript_id (str): the t_id to return rows for, e.g. "ENST00000234420.11".
        tool_name (str): "spliceai" or "pangolin".

    Returns:
        A flask Response.
    """
    if "error" in results:
        return error_response(results["error"], source=tool_name)

    if not transcript_id:
        return error_response('"transcript" not specified.\n', source=tool_name)

    if results.get("isPositionOnly"):
        return error_response("Per-position scores are not available for position-only queries.\n", source=tool_name)

    for transcript_scores in results.get("scores") or []:
        if transcript_scores.get("t_id") == transcript_id:
            return Response(json.dumps({
                "transcript": transcript_id,
                "strand": transcript_scores.get("t_strand"),
                "rows": transcript_scores.get("ALL_NON_ZERO_SCORES") or [],
            }), status=200, mimetype='application/json', headers=[
                ('Access-Control-Allow-Origin', '*'),
            ])

    return error_response(f'Transcript "{transcript_id}" not found in the results for this variant.\n', source=tool_name)


_FIRST_REQUEST_LOGGED = False


def run_splice_prediction_tool(tool_name, scores_for_one_transcript=False):
    """Handles API request for splice prediction.

    DB connections are taken from the pool only for short bursts (cache lookup,
    rate-limit check, log writes, cache writes, and per-call transcript-structure
    SELECTs inside get_spliceai_scores and get_pangolin_scores). The model inference
    runs without holding any pooled connection, so a slow inference can't starve the pool.

    Args:
        tool_name (str): "spliceai" or "pangolin"
        scores_for_one_transcript (bool): when True, respond with just the per-position rows
            for the transcript named by the "transcript" param (the /scores endpoints) instead
            of the full results. Everything up to that point -- param validation, rate limiting,
            the cache lookup and the cache write -- is shared with the main endpoint, so a
            /scores call for an already-computed variant is a cache read.
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

    # Each service pins one gene set (see GENE_SET) and loads only that annotator, so a request
    # for the other one has to go to the service that has it. The message names that service
    # rather than just refusing, since API clients that predate the split will land here.
    if basic_or_comprehensive_param != GENE_SET:
        other_service = f"{tool_name}-{GENOME_VERSION}" + ("" if basic_or_comprehensive_param == "basic" else "-comprehensive")
        return error_response(
            f'This service only handles bc={GENE_SET} requests, but received bc={basic_or_comprehensive_param}. '
            f'Send bc={basic_or_comprehensive_param} requests to the "{other_service}" service instead '
            f'(https://{other_service}-xwkwwwxdwq-uc.a.run.app/{tool_name}/).\n',
            source=tool_name)

    # A bare chrom-pos position (no ref/alt) requests REF-only scores instead
    # of the usual REF-vs-ALT delta scores. Try the full variant format first
    # so a well-formed variant is never misparsed as a position.
    try:
        parse_variant(variant)
        is_position_only = False
    except ValueError:
        try:
            parse_position(variant)
            is_position_only = True
        except ValueError:
            return error_response(
                f'Unable to parse "variant": "{variant}". Expected either a chrom-pos-ref-alt variant '
                f'(e.g. chr8-140300615-C-G) or a chrom-pos position (e.g. chr8-140300615) for REF-only scores.\n',
                source=tool_name)

    variant_consequence = params.get("variant_consequence")

    force = params.get("force")  # ie. don't use cache

    print(f"{logging_prefix}: ======================", flush=True)
    print(f"{logging_prefix}: {variant} tool={tool_name} hg={genome_version}, distance={distance_param}, mask={mask_param}, bc={basic_or_comprehensive_param}", flush=True)

    # Everything this service can serve was loaded at startup by preload_models(), so there is
    # no init_* call on the request path any more -- the first request is as fast as the
    # thousandth. That also means a request for the other gene set has nothing loaded to answer
    # it, hence the check above.

    # check cache before processing the variant (short DB scope)
    results = {}
    if not force:
        with get_db_connection() as conn:
            results = get_splicing_scores_from_cache(conn, tool_name, variant, genome_version, distance_param, mask_param, basic_or_comprehensive_param, is_position_only)

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
        # and get_pangolin_scores each acquire their own short-lived connection
        # internally for transcript-structure SELECTs after inference completes.
        try:
            if is_position_only:
                if tool_name == "spliceai":
                    results = get_spliceai_reference_scores(variant, genome_version, distance_param, basic_or_comprehensive_param)
                elif tool_name == "pangolin":
                    results = get_pangolin_reference_scores(variant, genome_version, distance_param, basic_or_comprehensive_param)
                else:
                    raise ValueError(f"Invalid tool_name: {tool_name}")
            elif tool_name == "spliceai":
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
        # Set by get_spliceai_scores or get_pangolin_scores when the transcript-structure
        # lookup couldn't run (no DB connection, or the query failed), so that enrichment
        # was skipped and the result reflects degraded inputs — don't cache it.
        skip_cache = results.pop("_skip_cache", False)

        # Post-inference: log + cache write + (if error) error log, all in one short DB scope.
        duration = (datetime.now() - start_time).total_seconds()
        with get_db_connection() as conn:
            log(conn, f"{tool_name}:computed", ip=user_ip, duration=duration, variant=variant, genome=genome_version, distance=distance_param, mask=mask_param, bc=basic_or_comprehensive_param, variant_consequence=variant_consequence)
            if "error" not in results and not skip_cache:
                add_splicing_scores_to_cache(conn, tool_name, variant, genome_version, distance_param, mask_param, basic_or_comprehensive_param, results, is_position_only)
            elif "error" in results:
                log(conn, f"{tool_name}:error", ip=user_ip, variant=variant, genome=genome_version, distance=distance_param, mask=mask_param, details=results["error"], bc=basic_or_comprehensive_param, variant_consequence=variant_consequence)

    if scores_for_one_transcript:
        return per_transcript_scores_response(results, params.get("transcript"), tool_name)

    # The per-position rows are kept in the cached copy (see get_spliceai_scores) but not sent
    # with the main response: the results table only needs nNonZeroScores, and the /scores
    # endpoints serve the rows for whichever transcript the user opens.
    for transcript_scores in results.get("scores") or []:
        transcript_scores.pop("ALL_NON_ZERO_SCORES", None)

    # Echo only a whitelist of input params back to the client. A prior version
    # used `response_json.update(params)`, which reflected every query-string
    # field unchanged — combined with the `.html()` sinks in index.html and
    # the URL-hash auto-submit on page load, that let a crafted link execute
    # arbitrary HTML in visitors' browsers.
    ECHO_PARAM_KEYS = (
        "variant", "hg", "bc", "distance", "mask", "raw", "variant_consequence",
    )
    response_json = {k: params[k] for k in ECHO_PARAM_KEYS if k in params}
    # REF-only (position-only) results never set a "mask" key (mask has no meaning
    # without an ALT allele) -- drop a client-supplied mask echo too, so the response
    # never implies a mask value was applied when it wasn't.
    if is_position_only:
        response_json.pop("mask", None)
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


@app.before_request
def block_ips():
    """Reject hard-blocked IPs before routing, DB access, or model inference.

    Returns a 429 at the earliest point in the request lifecycle so a flood from
    a blocked IP costs no DB connection and no compute — unlike the
    DB-backed restricted_ips check, which runs later and only on cache-miss.
    Returning None (for every other IP) lets the request proceed normally. The
    `BLOCKED_IPS and` guard skips the X-Forwarded-For parse when no IPs are
    blocked, so the hook is free in the common case.
    """
    if BLOCKED_IPS and get_user_ip(request) in BLOCKED_IPS:
        return error_response(RATE_LIMIT_ERROR_MESSAGE, status=429)


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


@app.route('/', strict_slashes=False)
def index():
    # Bare-root probes (uptime checkers, scanners polling `*.run.app/`) would
    # otherwise generate 404 noise that drowns out real client errors in the
    # monitoring breakdown. Return a tiny 200 so they look like normal traffic.
    return Response("OK\n", status=200, mimetype='text/plain')


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


def preload_models():
    """Load everything this service can serve, before it accepts any request.

    A service answers for exactly one (GENOME_VERSION, GENE_SET) pair, so the full set is known
    at startup and there is nothing left to load lazily. That is the point: the request path no
    longer calls init_*, so no user's request pays a model load, and a container that cannot
    load its models fails at startup instead of 500ing on the unlucky first request.

    This must run in the WORKER, not in the gunicorn arbiter, which is why the Dockerfiles no
    longer pass `--preload`: without it gunicorn imports this module separately in each worker
    after forking, so each worker builds its own copies. That is deliberate. Annotator.__init__
    opens the reference FASTA with pyfastx (an open file descriptor plus a sqlite .fxi index)
    and loads 5 Keras models, and neither a sqlite connection nor the TensorFlow runtime
    survives fork(): sharing them across forked workers hung every SpliceAI inference until
    gunicorn's --timeout killed the worker, turning each request into a 503. The DB pool comment
    above describes the same hazard for psycopg2, and takes the same way out.

    The cost of not forking shared copies is that each worker holds its own annotator. Pinning
    one gene set per service is what makes that affordable: a worker can no longer accumulate
    both, which is what exhausted the instance's memory before the split.
    """
    t0 = time.time()
    if TOOL == "spliceai":
        init_spliceai(GENOME_VERSION, GENE_SET)
    elif TOOL == "pangolin":
        init_pangolin()
    init_transcript_annotations(GENOME_VERSION, GENE_SET)
    print(f"[startup pid={os.getpid()}] preload_models() ready for hg{GENOME_VERSION} {GENE_SET} "
          f"in {time.time() - t0:.2f}s", flush=True)


preload_models()

# Move everything loaded above into the GC's permanent generation, which the collector never
# walks again. The transcript annotations alone are hundreds of thousands of nested dicts that
# live for the life of the process, and every full collection would otherwise re-traverse all
# of them to prove what is already known: none of it is garbage.
gc.freeze()

print(f"[startup pid={os.getpid()}] server.py module loaded in "
      f"{time.time() - _PROCESS_START_TIME:.2f}s (tool={TOOL}, genome={GENOME_VERSION})", flush=True)


# Start the Werkzeug dev server only when this file is run directly
# (python server.py). Under gunicorn the module is imported (so __name__ is
# "server", not "__main__") and gunicorn serves `app` itself. The previous
# `or os.environ.get('RUNNING_ON_GOOGLE_CLOUD_RUN')` clause fired app.run() at
# import time under gunicorn --preload, blocking the arbiter before it could fork
# workers — silently defeating the gunicorn worker recycling (--timeout 120) this
# image relies on to recover from stuck inferences.
if __name__ == '__main__':
    app.run(debug=DEBUG, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))
