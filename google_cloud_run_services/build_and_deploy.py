import argparse
import logging
import os
import platform
import subprocess
import time

from dotenv import load_dotenv
import pandas as pd
import psycopg2
from psycopg2.extras import execute_values
import re
from tqdm import tqdm

load_dotenv()

logging.basicConfig(level=logging.INFO, format="%(asctime)s: %(message)s")

VALID_COMMANDS = {
    "update_annotations", "update_transcript_tables", "build", "deploy", "test", "test2", "run",
}

GCLOUD_PROJECT = "spliceai-lookup-412920"
DOCKERHUB_REPO = "docker.io/weisburd"

DB_HOST = os.environ.get("SPLICEAI_LOOKUP_DB_HOST")
if not DB_HOST:
    raise ValueError("SPLICEAI_LOOKUP_DB_HOST not set. Please add it to .env file.")
DB_NAME = "spliceai-lookup-db"
DB_USER = "postgres"


def get_db_connection():
    """Get a database connection using password from .pgpass file."""
    pgpass_path = os.path.join(os.path.dirname(__file__), ".pgpass")
    if not os.path.exists(pgpass_path):
        raise FileNotFoundError(f"Database password file not found: {pgpass_path}")

    with open(pgpass_path, "r") as f:
        password = f.read().strip()

    return psycopg2.connect(
        host=DB_HOST,
        dbname=DB_NAME,
        user=DB_USER,
        password=password,
    )


def update_transcript_tables(genome_versions, gencode_version):
    """Populate the transcripts table in the database from genePred files.

    Args:
        genome_versions: List of genome versions to process (e.g., ["37", "38"])
        gencode_version: The gencode version string (e.g., "v49")
    """
    conn = get_db_connection()
    cursor = conn.cursor()

    # Process each genome version
    total_records = 0
    for genome_version in genome_versions:
        table_name = f"transcripts_hg{genome_version}"
        temp_table_name = f"{table_name}_reloading"

        # Create the temporary table for this genome version
        logging.info(f"Creating {temp_table_name} table...")
        cursor.execute(f"DROP TABLE IF EXISTS {temp_table_name}")
        cursor.execute(f"""
            CREATE TABLE {temp_table_name} (
                transcript_id VARCHAR(50) PRIMARY KEY,
                chrom VARCHAR(25) NOT NULL,
                strand VARCHAR(1) NOT NULL,
                tx_start INTEGER NOT NULL,
                tx_end INTEGER NOT NULL,
                cds_start INTEGER,
                cds_end INTEGER,
                exon_count INTEGER NOT NULL,
                exon_starts TEXT NOT NULL,
                exon_ends TEXT NOT NULL
            )
        """)
        conn.commit()
        # Look for genePred files in the annotations directory
        gene_pred_path = f"./docker/ref/GRCh{genome_version}/gencode.{gencode_version}.GRCh{genome_version}.sorted.txt.gz"

        if not os.path.exists(gene_pred_path):
            # Try alternate path patterns
            alt_path = f"gencode.{gencode_version}.GRCh{genome_version}.sorted.txt.gz"
            if os.path.exists(alt_path):
                gene_pred_path = alt_path
            else:
                logging.warning(f"GenePred file not found: {gene_pred_path}")
                continue

        logging.info(f"Reading genePred file: {gene_pred_path}")

        # genePred extended format columns
        column_names = [
            "name",        # transcript ID
            "chrom",
            "strand",
            "txStart",     # 0-based
            "txEnd",
            "cdsStart",    # 0-based
            "cdsEnd",
            "exonCount",
            "exonStarts",  # comma-separated, 0-based
            "exonEnds",    # comma-separated
            "score",
            "name2",       # gene name
            "cdsStartStat",
            "cdsEndStat",
            "exonFrames",  # comma-separated
        ]

        # The sorted file has an additional index column at the start
        df = pd.read_table(gene_pred_path, names=["i"] + column_names)

        # Prepare all rows for bulk insert
        rows = []
        for _, row in tqdm(df.iterrows(), total=len(df), desc=f"GRCh{genome_version}"):
            transcript_id = row["name"].split(".")[0]  # Remove version suffix
            cds_start = int(row["cdsStart"]) if pd.notna(row["cdsStart"]) else None
            cds_end = int(row["cdsEnd"]) if pd.notna(row["cdsEnd"]) else None

            # Check if CDS start equals CDS end (non-coding transcript)
            if cds_start is not None and cds_end is not None and cds_start == cds_end:
                cds_start = None
                cds_end = None

            rows.append((
                transcript_id,
                row["chrom"],
                row["strand"],
                int(row["txStart"]),
                int(row["txEnd"]),
                cds_start,
                cds_end,
                int(row["exonCount"]),
                row["exonStarts"],
                row["exonEnds"],
            ))

        # Bulk insert using execute_values (much faster than individual inserts)
        batch_size = 1000
        insert_sql = f"""INSERT INTO {temp_table_name}
               (transcript_id, chrom, strand, tx_start, tx_end,
                cds_start, cds_end, exon_count, exon_starts, exon_ends)
               VALUES %s
               ON CONFLICT (transcript_id) DO UPDATE SET
                   chrom = EXCLUDED.chrom,
                   strand = EXCLUDED.strand,
                   tx_start = EXCLUDED.tx_start,
                   tx_end = EXCLUDED.tx_end,
                   cds_start = EXCLUDED.cds_start,
                   cds_end = EXCLUDED.cds_end,
                   exon_count = EXCLUDED.exon_count,
                   exon_starts = EXCLUDED.exon_starts,
                   exon_ends = EXCLUDED.exon_ends"""

        for i in tqdm(range(0, len(rows), batch_size), desc=f"Inserting GRCh{genome_version}"):
            batch = rows[i:i + batch_size]
            execute_values(cursor, insert_sql, batch)

        conn.commit()

        # Create index on transcript_id for fast lookups
        logging.info(f"Creating index on {temp_table_name}...")
        cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{temp_table_name}_tid ON {temp_table_name} (transcript_id)")
        conn.commit()

        # Replace the old table with the new one
        logging.info(f"Replacing {table_name} with new data...")
        cursor.execute(f"DROP TABLE IF EXISTS {table_name}")
        cursor.execute(f"ALTER TABLE {temp_table_name} RENAME TO {table_name}")
        cursor.execute(f"ALTER INDEX IF EXISTS idx_{temp_table_name}_tid RENAME TO idx_{table_name}_tid")
        conn.commit()

        logging.info(f"Inserted {len(rows):,d} records into {table_name}")
        total_records += len(rows)

    cursor.close()
    conn.close()

    logging.info(f"Done! Inserted {total_records:,d} total transcript records into the database.")


def get_service_name(tool, genome_version, gene_set="basic"):
    """Cloud Run service name for one (tool, genome, gene set) combination.

    The basic services keep their original un-suffixed names so the API URLs published in the
    docs, and used by external clients, keep working. Only the comprehensive services, which
    are new, carry a suffix.
    """
    return f"{tool}-{genome_version}" + ("" if gene_set == "basic" else f"-{gene_set}")

def get_tag(tool, genome_version, repo_name="gcr.io"):
    if repo_name == "gcr.io":
        return f"us-central1-docker.pkg.dev/spliceai-lookup-412920/docker/{get_service_name(tool, genome_version)}"
    elif repo_name == "dockerhub":
        return f"{DOCKERHUB_REPO}/{get_service_name(tool, genome_version)}"
    else:
        raise ValueError(f"Invalid repo_name arg: {repo_name}")

def current_git_commit():
    """The commit this checkout is on, marked when the working tree is dirty.

    Recorded alongside a built image's digest purely so that the source of a promoted image can be
    identified afterwards. It is deliberately NOT what promotion checks: "-dirty" is the same
    string for every dirty tree, and an unrelated commit between building and promoting would move
    it, so it identifies a build too loosely in one direction and too strictly in the other.
    """
    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain"], capture_output=True, text=True, check=True).stdout.strip()
    return f"{commit}-dirty" if dirty else commit


def run(c):
    """Run a shell command, raising if it fails.

    os.system's exit status used to be discarded, so a failed `gcloud run deploy` looked
    exactly like a successful one: a deploy that never happened was reported green by CI, and
    the missing service was only noticed by hand afterwards. Raising makes the CI job fail on
    the step that actually broke.
    """
    logging.info(c)
    # os.system returns a wait status, not an exit code: the low byte carries the signal that
    # killed the process and the high byte the exit code, so a plain `!= 0` would report the
    # right thing for the wrong reason and mangle the code it prints.
    status = os.system(c)
    exit_code = os.waitstatus_to_exitcode(status)
    if exit_code != 0:
        raise RuntimeError(f"Command failed with exit code {exit_code}: {c}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genome-version", choices=["37", "38"], help="If not specified, command will run for both GRCh37 and GRCh38")
    parser.add_argument("-t", "--tool", choices=["spliceai", "pangolin"], help="If not specified, command will run for both spliceai and pangolin")
    parser.add_argument("-s", "--gene-set", choices=["basic", "comprehensive"],
                        help="Which Gencode gene set's service to deploy. If not specified, deploys both. "
                             "One image serves either, so this only affects the deploy step, not the build.")
    parser.add_argument("-d", "--docker-command", choices=["docker", "podman"], default="docker", help="Whether to use docker or podman to build the image")
    g = parser.add_mutually_exclusive_group()
    g.add_argument("--gencode-version",
                   help="The gencode version to use for the 'update_annotations' command (example: 'v49'). Either this "
                        "or --gencode-gtf must be specified for the 'update_annotations' command")
    g.add_argument("--gencode-gtf",
                   help="Path of the newest 'basic' Gencode GTF file that was downloaded from "
                        "https://www.gencodegenes.org/human/. Either this or --gencode-version must be specified for "
                        "the 'update_annotations' command")

    parser.add_argument("--dev", action="store_true",
                        help="Deploy as a 'dev'-tagged revision with --no-traffic so production keeps "
                             "serving the existing revision. Test against the tagged URL printed by "
                             "gcloud, then promote with --promote. Do not promote by switching traffic "
                             "to the dev revision: DEPLOYMENT=dev is baked into that revision, so "
                             "production would keep using the dev cache namespace.")
    parser.add_argument("--promote", action="store_true",
                        help="Promote the image a --dev deploy already tested: skip the build, read "
                             "docker/<tool>/sha256_grch<hg>_dev.txt, and deploy that exact digest as a "
                             "DEPLOYMENT=prod revision that takes traffic. Use this rather than a plain "
                             "re-run without --dev, which rebuilds: the base image tag and the unpinned "
                             "pip requirements resolve to whatever is current at build time, so the same "
                             "checkout can produce a different image from the one that was tested.")
    parser.add_argument("command", nargs="?", choices=VALID_COMMANDS,
                        help="Command to run. If not specified, it will run 'build' and then 'deploy'")

    args = parser.parse_args()

    if args.promote and args.dev:
        parser.error("--promote and --dev are opposites: --promote deploys the digest that a previous "
                     "--dev run recorded, as production.")

    if args.genome_version:
        genome_versions = [args.genome_version]
    else:
        genome_versions = ["38", "37"]

    if args.tool:
        tools = [args.tool]
    else:
        tools = ["spliceai", "pangolin"]

    if args.gene_set:
        gene_sets = [args.gene_set]
    else:
        gene_sets = ["basic", "comprehensive"]

    if args.gencode_version:
        if not re.fullmatch(r"v\d+", args.gencode_version):
            parser.error("--gencode-version must be of the form 'v46'")
        gencode_version_number = int(args.gencode_version.lstrip("v"))
    else:
        gencode_version_number = None

    if args.command == "update_annotations":
        if not args.gencode_version and not args.gencode_gtf:
            parser.error("Either --gencode-version or --gencode-gtf must be specified for the update_annotations command")

        gencode_gtf_paths = {}
        if args.gencode_version:
            for genome_version in genome_versions:
                for basic_or_comprehensive in "", ".basic":
                    if genome_version == "37":
                        gencode_gtf_url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_version_number}/GRCh37_mapping/gencode.{args.gencode_version}lift37{basic_or_comprehensive}.annotation.gtf.gz"
                    elif genome_version == "38":
                        gencode_gtf_url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_version_number}/gencode.{args.gencode_version}{basic_or_comprehensive}.annotation.gtf.gz"
                    else:
                        parser.error(f"Invalid genome version: {genome_version}")

                    run(f"wget -nc {gencode_gtf_url}")
                    # Pin to v369 on Linux — the "latest" binary requires glibc >=2.32 which is too new for Debian 11 / Ubuntu 20.04.
                    gtf_to_gene_pred_os_dir = "macOSX.x86_64" if platform.system() == "Darwin" else "linux.x86_64.v369"
                    run(f"wget -nc https://hgdownload.soe.ucsc.edu/admin/exe/{gtf_to_gene_pred_os_dir}/gtfToGenePred")
                    run(f"chmod 777 gtfToGenePred")
                    gencode_gtf_paths[(genome_version, basic_or_comprehensive)] = os.path.basename(gencode_gtf_url)
        else:
            if not args.genome_version:
                parser.error("If --gencode-gtf is specified, --genome-version is required")
            if not os.path.isfile(args.gencode_gtf):
                parser.error(f"File not found: {args.gencode_gtf}")
            gencode_gtf_paths[(args.genome_version, "basic")] = args.gencode_gtf

        # Deduplicated: gencode_gtf_paths holds one key per (genome_version, gene set), so
        # iterating it directly ran each rm twice per genome. `rm -f` because a glob that matches
        # nothing is the normal state on a clean checkout and must not abort the run now that
        # run() raises on a non-zero exit.
        for genome_version in sorted({gv for gv, _ in gencode_gtf_paths}):
            # Deliberately still only ".basic.": --gencode-gtf mode regenerates just the basic
            # transcript_annotations file, so widening this to match the comprehensive one too
            # would delete a file that mode never rebuilds. (The comprehensive file from a
            # previous Gencode release is therefore left in place here -- pre-existing, and not
            # something this change set set out to alter.)
            run(f"rm -f ./docker/ref/GRCh{genome_version}/gencode.*.basic.annotation.transcript_annotations.json.gz")
            run(f"rm -f ./docker/spliceai/annotations/GRCh{genome_version}/gencode.*.annotation*.txt.gz")
            run(f"rm -f ./docker/pangolin/annotations/GRCh{genome_version}/gencode.*.annotation*.db")

        for (genome_version, basic_or_comprehensive), gencode_gtf_path in gencode_gtf_paths.items():
            # generate genePred files to use as gene tracks in IGV.js
            if args.gencode_version:
                gene_pred_path = f"gencode.{args.gencode_version}.GRCh{genome_version}.txt"
                run(f"./gtfToGenePred -genePredExt -geneNameAsName2 {gencode_gtf_path} {gene_pred_path}")

                print(f"Reading {gene_pred_path}")
                column_names = [
                    "name",
                    "chrom",
                    "strand",
                    "txStart",
                    "txEnd",
                    "cdsStart",
                    "cdsEnd",
                    "exonCount",
                    "exonStarts",
                    "exonEnds",
                    "score",
                    "name2",
                    "cdsStartStat",
                    "cdsEndStat",
                    "exonFrames",
                ]
                df = pd.read_table(gene_pred_path, names=column_names)
                df["txStart"] = df["txStart"].astype(int)
                df["txEnd"] = df["txEnd"].astype(int)
                filter_exp = (df["txStart"] > 0) & (df["txEnd"] > 0)
                df2 = df[filter_exp]
                if len(df) - len(df2) > 0:
                    print(f"Filtered out {len(df) - len(df2):,d} records from {gene_pred_path}:")
                    print(df[~filter_exp])

                df2 = df2.sort_values(["chrom", "txStart", "txEnd"])
                df2["i"] = df2["name2"].map({name: i for i, name in enumerate(df2.name2.unique())})
                df2 = df2[["i"] + column_names]
                sorted_gene_pred_path = gene_pred_path.replace(".txt", ".sorted.txt")
                df2.to_csv(sorted_gene_pred_path, header=False, index=False, sep="\t")
                run(f"bgzip -f {sorted_gene_pred_path}")
                run(f"tabix -s 3 -b 5 -e 6 -f {sorted_gene_pred_path}.gz")

                run(f"gsutil -m cp {sorted_gene_pred_path}.gz* gs://tgg-viewer/ref/GRCh{genome_version}/gencode_{args.gencode_version}/")

            # generate SpliceAI annotation files
            run(f"python3 ../annotations/generate_transcript_annotation_json.py {gencode_gtf_path}")
            output_json_path = gencode_gtf_path.replace(".gtf.gz", ".transcript_annotations.json.gz")
            run(f"python3 ../annotations/convert_gtf_to_SpliceAI_annotation_input_format.py -a {output_json_path} {gencode_gtf_path}")
            if not os.path.isfile(output_json_path):
                raise ValueError(f"Unable to find {output_json_path}")

            run(f"mv {output_json_path} ./docker/ref/GRCh{genome_version}/")
            run(f"mv {gencode_gtf_path.replace('.gtf.gz', '.txt.gz')} ./docker/spliceai/annotations/GRCh{genome_version}/")

            if genome_version == "37":
                gencode_gtf_path_without_chr_prefix = gencode_gtf_path.replace(".gtf.gz", ".without_chr_prefix.gtf.gz")
                run(f"gunzip -c {gencode_gtf_path} | sed 's/chr//g' | bgzip > {gencode_gtf_path_without_chr_prefix}")
                gencode_gtf_path = gencode_gtf_path_without_chr_prefix

            # generate Pangolin annotation files
            run(f"python3 create_pangolin_db.py {gencode_gtf_path}")
            run(f"mv {gencode_gtf_path.replace('.gtf.gz', '.db')} ./docker/pangolin/annotations/GRCh{genome_version}/")

        if args.gencode_version:
            with open("server.py", "rt") as f:
                server_py = f.readlines()

            updated_line = False
            with open("server.py", "wt") as f:
                for i, line in enumerate(server_py):
                    if line.startswith("GENCODE_VERSION ="):
                        new_gencode_line = f"GENCODE_VERSION = \"{args.gencode_version}\""
                        f.write(f"{new_gencode_line}\n")
                        updated_line = True
                        print(f"Updated server.py line #{i} to {new_gencode_line}")
                    else:
                        f.write(line)

            with open("../index.html", "rt") as f:
                index_html = f.readlines()

            updated_line = False
            with open("../index.html", "wt") as f:
                for i, line in enumerate(index_html):
                    if "const GENCODE_VERSION = " in line:
                        new_gencode_line = f"\tconst GENCODE_VERSION = \"{args.gencode_version}\""
                        f.write(f"{new_gencode_line}\n")
                        updated_line = True
                        print(f"Updated index.html line #{i} to {new_gencode_line}")
                    else:
                        f.write(line)

            if not updated_line:
                print("WARNING: Unable to find GENCODE_VERSION line in index.html")

        return

    if args.command == "update_transcript_tables":
        if not args.gencode_version:
            parser.error("--gencode-version is required for the update_transcript_tables command")

        update_transcript_tables(genome_versions, args.gencode_version)
        return

    if args.command == "test2":
        run(f"gcloud beta code dev")
        return

    if args.command in {"test", "run"}:
        if not args.genome_version:
            parser.error(f"--genome-version is required for the {args.command} command")
        if not args.tool:
            parser.error(f"--tool is required for the {args.command} command")

        tag = get_tag(args.tool, args.genome_version)

        if args.command == "run":
            print("Run this command: ")
            print(f"{args.docker_command} run -it {tag}:latest /bin/bash")
        elif args.command == "test":
            run(f"{args.docker_command} run -p 8080:8080 {tag}:latest")

        return

    if not args.command or args.command in {"build", "deploy"}:
        if args.docker_command == "podman":
            print("WARNING: Google Cloud Run doesn't appear to work with images built using podman. "
                  "Containers may fail to deploy to Google Cloud Run unless they are built using docker.")
            time.sleep(10)

        for genome_version in genome_versions:
            for tool in tools:
                tag = get_tag(tool, genome_version)
                dockerhub_tag = get_tag(tool, genome_version, repo_name="dockerhub")
                # Names the image, and so the per-image buildx cache scope below. The services
                # deployed from it are named per gene set inside the deploy block.
                image_name = get_service_name(tool, genome_version)
                # gunicorn worker processes per instance (the `--workers` count baked into
                # the image). Each worker is single-threaded (TF/torch thread pools are
                # capped to 1 in the Dockerfile) so `workers` alone sets how many model
                # inferences run at once. Lowered 6 -> 3 to ease CPU oversubscription on the
                # 2-CPU instances: at 6 workers, >2 simultaneous cache-miss inferences shared
                # 2 cores so each ~6s prediction stretched ~3x and tripped the 120s gunicorn
                # timeout (SIGABRT -> 500/503). 3 workers is a modest 1.5x over 2 cores, so a
                # single slow inference can't starve the others into the timeout; extra burst
                # concurrency comes from scaling out across instances, not more workers.
                # This is the count baked into the image; the comprehensive services override it
                # per deploy below. Cloud Run's --concurrency is always set to match, since a
                # sync gunicorn worker serves one request at a time and any excess concurrency
                # would queue inside the instance rather than scaling out.
                workers = 3
                min_instances = 0  # if tool == 'pangolin' else 2
                # Raised 3 -> 6 so bursts of cache-miss traffic scale out across instances
                # instead of saturating a few and timing out (504). This is only a ceiling:
                # min_instances=0 means idle services still scale to zero, so the baseline
                # cost is unchanged.
                max_instances = 6
                # Keep dev image digests separate from prod so a stray non-dev
                # deploy from the same checkout can't accidentally promote the
                # dev image.
                sha256_path = f"docker/{tool}/sha256_grch{genome_version}{'_dev' if args.dev or args.promote else ''}.txt"
                if (not args.command or args.command == "build") and not args.promote:
                    if args.docker_command == "podman":
                        run(f"gcloud --project {GCLOUD_PROJECT} auth print-access-token | podman login -u oauth2accesstoken --password-stdin us-central1-docker.pkg.dev")

                    # In CI (USE_BUILDX_CACHE set; see deploy-on-tag.yml) build with buildx and a
                    # per-service GitHub Actions layer cache, so unchanged layers (base image,
                    # tensorflow-cpu, pip deps) restore instead of rebuilding from scratch on the
                    # fresh runner. `scope={image_name}` keeps the 4 matrix builds from clobbering each
                    # other's cache. --load puts the built image into the local docker store so the
                    # push / pull / digest-capture steps below work unchanged. Locally the buildx
                    # container driver and the gha backend aren't set up, so fall back to a plain build.
                    build_cmd = (
                        f"docker buildx build --load "
                        f"--cache-from type=gha,scope={image_name} --cache-to type=gha,mode=max,scope={image_name} "
                        if os.environ.get("USE_BUILDX_CACHE") else f"{args.docker_command} build "
                    )
                    run(f"{build_cmd}-f docker/{tool}/Dockerfile --build-arg=\"WORKERS={workers}\" --build-arg=\"GENOME_VERSION={genome_version}\" -t {tag}:latest -t {dockerhub_tag}:latest .")
                    run(f"{args.docker_command} push {tag}:latest")
                    run(f"{args.docker_command} push {dockerhub_tag}:latest")

                    run(f"{args.docker_command} pull {tag}:latest")
                    run(f"{args.docker_command} inspect --format='{{{{range .RepoDigests}}}}{{{{println .}}}}{{{{end}}}}' {tag}:latest | grep 'us-central1-docker.pkg.dev' | cut -f 2 -d @ > {sha256_path}")  # record the image's sha256

                if not args.command or args.command == "deploy":
                    with open(sha256_path) as f:
                        sha256 = f.read().strip()

                    if not re.match("^sha256:[a-f0-9]{64}$", sha256):
                        raise ValueError(f"Invalid sha256 value found in {sha256_path}: {sha256}")

                    # One image, deployed once per gene set. The image carries both gene sets'
                    # annotation files (the Dockerfile COPYs per genome, not per gene set), so
                    # the two services differ only by the GENE_SET env var that tells server.py
                    # which single annotator to load. Deploying rather than rebuilding keeps the
                    # slow part -- 4 multi-GB image builds -- unchanged by the split.
                    for gene_set in gene_sets:
                        service = get_service_name(tool, genome_version, gene_set)

                        # One stamp per SERVICE, not per image. The two gene sets share an image
                        # but are separate Cloud Run services deployed separately, and a --dev run
                        # can succeed for one and fail for the other. A single per-image stamp
                        # would then have been written by the service that worked and would
                        # authorize promoting the one that never deployed.
                        built_here_path = f"{sha256_path}.built_here.{gene_set}"

                        if args.promote:
                            # Promote only a digest a --dev run on this machine actually deployed
                            # to this service. The stamp is gitignored, so its presence cannot come
                            # from a checkout: no stamp means no local dev deploy of this service,
                            # and a stamp naming a different digest means the file moved on since.
                            built_here = None
                            if os.path.exists(built_here_path):
                                with open(built_here_path) as f:
                                    built_here = f.readline().strip()
                            if built_here != sha256:
                                raise RuntimeError(
                                    f"Refusing to promote {sha256} to {service}: no --dev deploy on this "
                                    f"machine recorded it ({built_here_path} "
                                    f"{'names ' + built_here if built_here else 'is missing'}). Promotion "
                                    f"deploys that digest to production and sends it all traffic, and these "
                                    f"digest files are tracked in git, so a fresh checkout holds whatever "
                                    f"was committed last rather than anything built here. Run the --dev "
                                    f"deploy first, test it, then promote.")

                        # The comprehensive gene set carries several times more transcripts, so
                        # both its annotator and the working memory of a request over a
                        # many-isoform gene are larger, and 3 workers of it on a 4Gi instance
                        # was enough for the kernel to OOM-kill one mid-request. Its share of
                        # traffic is ~2%, so trading per-instance concurrency for headroom costs
                        # little and bursts still scale out to max_instances. WORKERS is read
                        # from the environment by the image's gunicorn command, so this overrides
                        # the value baked in at build time without needing a separate image.
                        service_workers = workers if gene_set == "basic" else 2

                        if args.dev:
                            # Whether the service already exists decides whether --no-traffic can
                            # be passed at all: Cloud Run rejects it when it has to create the
                            # service, since a first revision has nothing to hold traffic back
                            # from. Only the --dev path passes --no-traffic, so only it needs to
                            # ask -- probing on the production path would let a transient gcloud
                            # failure abort a deploy whose behavior does not depend on the answer.
                            #
                            # Only a confirmed "not found" counts as absent. Treating every
                            # non-zero exit as absence would fail open in the worst possible
                            # direction: an expired credential or a transient API error during a
                            # --dev deploy of an EXISTING service would drop --no-traffic and hand
                            # an untested dev revision 100% of production traffic. Anything that is
                            # not a clean "exists" or a clean "not found" aborts the deploy instead.
                            probe = subprocess.run(
                                ["gcloud", f"--project={GCLOUD_PROJECT}", "run", "services", "describe", service,
                                 "--region=us-central1", "--format=value(name)"],
                                capture_output=True, text=True)
                            if probe.returncode == 0:
                                service_exists = True
                            # "Cannot find service [x]" is what gcloud actually prints for an absent
                            # service; the other spellings are kept in case the wording changes.
                            elif re.search(r"cannot find service|NOT_FOUND|could not be found|does not exist",
                                           probe.stderr, re.IGNORECASE):
                                service_exists = False
                            else:
                                raise RuntimeError(
                                    f"Could not determine whether the Cloud Run service {service} exists "
                                    f"(gcloud exit {probe.returncode}). Refusing to deploy, since guessing "
                                    f"wrong would route production traffic to this revision. stderr:\n{probe.stderr}")

                            if not service_exists:
                                # Refuse rather than create it here. Without --no-traffic the created
                                # revision would serve 100% of traffic while stamped DEPLOYMENT=dev,
                                # so every request to it would read and write the "__dev" cache
                                # namespace (get_splicing_scores_cache_key in server.py) and nothing
                                # in the --dev path would ever move traffic off it. That is not
                                # hypothetical for a service the frontend already points at:
                                # index.html hardcodes all eight hostnames, so a service becomes
                                # reachable the moment it exists, and a tag containing "dev" runs
                                # only the --dev half of the deploy workflow, leaving it that way
                                # until someone pushes a tag without "dev" in it.
                                raise RuntimeError(
                                    f"The Cloud Run service {service} does not exist yet, and creating it "
                                    f"from a --dev run would leave its first revision serving all traffic "
                                    f"with DEPLOYMENT=dev. Create it with a production deploy first "
                                    f"(python3 build_and_deploy.py -t {tool} -g {genome_version} "
                                    f"-s {gene_set}), then re-run with --dev.")

                            dev_flags = "--tag dev --no-traffic "
                            traffic_note = " (dev revision, no traffic)"
                        else:
                            dev_flags = ""
                            traffic_note = ""
                        print(f"Deploying {service} (GENE_SET={gene_set}, WORKERS={service_workers}) "
                              f"with image sha256 {sha256}{traffic_note}")
                        # Comma-separated list of IPs to hard-block at the API door
                        # (server.py block_ips). Sourced from the BLOCKED_IPS
                        # env var (set as a GitHub Actions repo Variable, or exported
                        # locally); unset -> empty, which clears any previous value and
                        # disables blocking. The "^@^" prefix tells gcloud to split
                        # env-var assignments on "@" instead of ",", so the commas
                        # between IPs stay inside the single BLOCKED_IPS value.
                        # GENE_SET pins the gene set, and DEPLOYMENT keeps a dev revision's
                        # cache entries off the keys production reads (server.py
                        # get_splicing_scores_cache_key).
                        env_vars = (f'BLOCKED_IPS={os.environ.get("BLOCKED_IPS", "").strip()}'
                                    f'@GENE_SET={gene_set}'
                                    f'@DEPLOYMENT={"dev" if args.dev else "prod"}'
                                    f'@WORKERS={service_workers}')
                        # Attach the Cloud SQL instance explicitly. The original four services
                        # carry this annotation from however they were first set up, and
                        # `gcloud run deploy` preserves it for them, so its absence here went
                        # unnoticed until the split created services that did not have it: they
                        # came up with no database at all, which silently disables response
                        # caching, per-IP rate limiting and transcript-structure enrichment.
                        run(f"""gcloud \
--project {GCLOUD_PROJECT} beta run deploy {service} \
--image {tag}@{sha256} \
--add-cloudsql-instances {GCLOUD_PROJECT}:us-central1:spliceai-lookup-db \
--min-instances {min_instances} \
--service-min-instances {min_instances} \
--max-instances {max_instances} \
--concurrency {service_workers} \
--service-account 1042618492363-compute@developer.gserviceaccount.com \
--execution-environment gen2 \
--region us-central1 \
--update-secrets=DB_PASSWORD=spliceai-lookup-db-password:2 \
--allow-unauthenticated \
--memory 4Gi \
--cpu 2 \
--cpu-boost \
--timeout 900s \
--update-env-vars "^@^{env_vars}" {dev_flags}""")

                        if args.dev:
                            # Promote by re-running this script with --promote, not by switching
                            # traffic to the dev revision. DEPLOYMENT is baked into a revision's
                            # environment, so a traffic-only promotion would leave production
                            # running a revision stamped DEPLOYMENT=dev: it would then read and
                            # write the "__dev"-suffixed cache keys forever, never seeing the
                            # entries production already has (see get_splicing_scores_cache_key).
                            # A DEPLOYMENT=prod revision of the same image digest is otherwise
                            # identical, and takes traffic. --promote is what deploys that exact
                            # digest: it reads sha256_grch<hg>_dev.txt, the image this run just
                            # pushed. Re-running without --dev instead would rebuild, and a
                            # rebuild from the same checkout is not the same image -- the
                            # Dockerfile's base tag and the unpinned Flask/gunicorn requirements
                            # both resolve to whatever is current at that moment -- so production
                            # could end up running something the dev testing never covered.
                            # Written here, after the dev deploy has actually succeeded (run()
                            # raises otherwise), rather than at build time: --promote means "ship
                            # the image dev has been running", so the thing it checks should record
                            # a completed dev deploy, not merely a local build.
                            #
                            # The stamp records the digest rather than the commit that produced it.
                            # Commit identity looks like the more natural choice and is the wrong
                            # one: the build rewrites the tracked sha256_*_dev.txt, and committing
                            # that file between testing a dev revision and promoting it is this
                            # repo's normal workflow, which would move HEAD and make the guard
                            # reject an image dev genuinely ran. The digest does not move when
                            # unrelated commits do. The commit goes on a second line for the
                            # record, and is not what the check compares.
                            with open(built_here_path, "w") as stamp:
                                stamp.write(f"{sha256}\n{current_git_commit()}\n")
                        else:
                            # Required when a previous --dev deploy left the service in manual-traffic mode; otherwise `gcloud run deploy` keeps traffic on the old revision.
                            run(f"gcloud --project {GCLOUD_PROJECT} run services update-traffic {service} "
                                f"--region us-central1 --to-latest")

                                # --add-volume=name=ref,type=cloud-storage,bucket=spliceai-lookup-reference-data,readonly=true \
                # --add-volume-mount=volume=ref,mount-path=/ref \

                    if args.dev:
                        # One hint per image, printed after every service's stamp is written. It
                        # repeats -s only when this run was narrowed to one gene set. A promote
                        # narrowed with -s skips the production digest write below (see the
                        # comment there), which is right after a narrowed dev run but would leave
                        # docker/<tool>/sha256_grch<hg>.txt stale if a run that deployed both
                        # services were followed by two narrowed promotes, one per service.
                        services = " and ".join(get_service_name(tool, genome_version, gene_set) for gene_set in gene_sets)
                        print(f"To promote {services} to production, deploy the digest this run recorded:")
                        print(f"  python3 build_and_deploy.py -t {tool} -g {genome_version}"
                              f"{' -s ' + args.gene_set if args.gene_set else ''} --promote")

                    # Record the promoted digest as the production one. Written once after the
                    # gene-set loop rather than inside it, because the path has no gene-set
                    # component while the two gene sets are separate services: a run narrowed to
                    # one of them with -s leaves the other on its previous digest, so writing from
                    # inside the loop would record a claim about a service this run never touched.
                    # That claim is not merely informational -- a later plain `deploy` reads this
                    # same file (sha256_path above) and would push the untouched service to an
                    # image it was never promoted to. Written after the deploys and traffic
                    # switches have all succeeded, so a failure partway through cannot leave the
                    # repository claiming production runs an image that never got there.
                    if args.promote and not args.gene_set:
                        with open(f"docker/{tool}/sha256_grch{genome_version}.txt", "w") as f:
                            f.write(f"{sha256}\n")

if __name__ == "__main__":
    main()
