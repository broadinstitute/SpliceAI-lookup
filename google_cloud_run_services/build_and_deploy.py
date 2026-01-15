import argparse
import logging
import os
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
                exon_ends TEXT NOT NULL,
                exon_frames TEXT
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
                row["exonFrames"] if pd.notna(row["exonFrames"]) else None,
            ))

        # Bulk insert using execute_values (much faster than individual inserts)
        batch_size = 1000
        insert_sql = f"""INSERT INTO {temp_table_name}
               (transcript_id, chrom, strand, tx_start, tx_end,
                cds_start, cds_end, exon_count, exon_starts, exon_ends, exon_frames)
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
                   exon_ends = EXCLUDED.exon_ends,
                   exon_frames = EXCLUDED.exon_frames"""

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


def get_service_name(tool, genome_version):
    return f"{tool}-{genome_version}"

def get_tag(tool, genome_version, repo_name="gcr.io"):
    if repo_name == "gcr.io":
        return f"us-central1-docker.pkg.dev/spliceai-lookup-412920/docker/{get_service_name(tool, genome_version)}"
    elif repo_name == "dockerhub":
        return f"{DOCKERHUB_REPO}/{get_service_name(tool, genome_version)}"
    else:
        raise ValueError(f"Invalid repo_name arg: {repo_name}")

def run(c):
    logging.info(c)
    os.system(c)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genome-version", choices=["37", "38"], help="If not specified, command will run for both GRCh37 and GRCh38")
    parser.add_argument("-t", "--tool", choices=["spliceai", "pangolin"], help="If not specified, command will run for both spliceai and pangolin")
    parser.add_argument("-d", "--docker-command", choices=["docker", "podman"], default="docker", help="Whether to use docker or podman to build the image")
    g = parser.add_mutually_exclusive_group()
    g.add_argument("--gencode-version",
                   help="The gencode version to use for the 'update_annotations' command (example: 'v49'). Either this "
                        "or --gencode-gtf must be specified for the 'update_annotations' command")
    g.add_argument("--gencode-gtf",
                   help="Path of the newest 'basic' Gencode GTF file that was downloaded from "
                        "https://www.gencodegenes.org/human/. Either this or --gencode-version must be specified for "
                        "the 'update_annotations' command")

    parser.add_argument("command", nargs="?", choices=VALID_COMMANDS,
                        help="Command to run. If not specified, it will run 'build' and then 'deploy'")

    args = parser.parse_args()

    if args.genome_version:
        genome_versions = [args.genome_version]
    else:
        genome_versions = ["38", "37"]

    if args.tool:
        tools = [args.tool]
    else:
        tools = ["spliceai", "pangolin"]

    if args.gencode_version:
        if not re.match("v[0-9][0-9]", args.gencode_version):
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
                    run(f"wget -nc https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/gtfToGenePred")
                    run(f"chmod 777 gtfToGenePred")
                    gencode_gtf_paths[(genome_version, basic_or_comprehensive)] = os.path.basename(gencode_gtf_url)
        else:
            if not args.genome_version:
                parser.error("If --gencode-gtf is specified, --genome-version is required")
            if not os.path.isfile(args.gencode_gtf):
                parser.error(f"File not found: {args.gencode_gtf}")
            gencode_gtf_paths[(args.genome_version, "basic")] = args.gencode_gtf

        for genome_version, _ in gencode_gtf_paths.keys():
            run(f"rm ./docker/ref/GRCh{genome_version}/gencode.*.basic.annotation.transcript_annotations.json.gz")
            run(f"rm ./docker/spliceai/annotations/GRCh{genome_version}/gencode.*.annotation*.txt.gz")
            run(f"rm ./docker/pangolin/annotations/GRCh{genome_version}/gencode.*.annotation*.db")

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
                run(f"gzcat {gencode_gtf_path} | sed 's/chr//g' | bgzip > {gencode_gtf_path_without_chr_prefix}")
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
                service = get_service_name(tool, genome_version)
                concurrency = 6    # if genome_version == '37' else 2
                min_instances = 0  # if tool == 'pangolin' else 2
                max_instances = 3
                if not args.command or args.command == "build":
                    if args.docker_command == "podman":
                        run(f"gcloud --project {GCLOUD_PROJECT} auth print-access-token | podman login -u oauth2accesstoken --password-stdin us-central1-docker.pkg.dev")

                    run(f"{args.docker_command} build -f docker/{tool}/Dockerfile --build-arg=\"CONCURRENCY={concurrency}\" --build-arg=\"GENOME_VERSION={genome_version}\" -t {tag}:latest -t {dockerhub_tag}:latest .")
                    run(f"{args.docker_command} push {tag}:latest")
                    run(f"{args.docker_command} push {dockerhub_tag}:latest")

                    run(f"{args.docker_command} pull {tag}:latest")
                    run(f"{args.docker_command} inspect --format='{{{{index .RepoDigests 0}}}}' {tag}:latest | cut -f 2 -d @ > docker/{tool}/sha256_grch{genome_version}.txt")  # record the image's sha256

                if not args.command or args.command == "deploy":
                    with open(f"docker/{tool}/sha256_grch{genome_version}.txt") as f:
                        sha256 = f.read().strip()

                    if not re.match("^sha256:[a-f0-9]{64}$", sha256):
                        raise ValueError(f"Invalid sha256 value found in docker/{tool}/sha256_grch{genome_version}.txt: {sha256}")

                    print(f"Deploying {service} with image sha256 {sha256}")

                    run(f"""gcloud \
--project {GCLOUD_PROJECT} beta run deploy {service} \
--image {tag}@{sha256} \
--min-instances {min_instances} \
--service-min-instances {min_instances} \
--max-instances {max_instances} \
--concurrency {concurrency} \
--service-account 1042618492363-compute@developer.gserviceaccount.com \
--execution-environment gen2 \
--region us-central1 \
--update-secrets=DB_PASSWORD=spliceai-lookup-db-password:2 \
--allow-unauthenticated \
--memory 4Gi \
--cpu 4
""")

                                # --add-volume=name=ref,type=cloud-storage,bucket=spliceai-lookup-reference-data,readonly=true \
                # --add-volume-mount=volume=ref,mount-path=/ref \

if __name__ == "__main__":
    main()
