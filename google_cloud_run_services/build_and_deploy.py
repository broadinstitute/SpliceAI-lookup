import argparse
import logging
import os
import pandas as pd
import re

logging.basicConfig(level=logging.INFO, format="%(asctime)s: %(message)s")

VALID_COMMANDS = {
    "update_annotations", "build", "deploy", "test", "test2", "run",
}

GCLOUD_PROJECT = "spliceai-lookup-412920"

def get_service_name(tool, genome_version):
    return f"{tool}-{genome_version}"

def get_tag(tool, genome_version):
    return f"us-central1-docker.pkg.dev/spliceai-lookup-412920/docker/{get_service_name(tool, genome_version)}"

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
                   help="The gencode version to use for the 'update_annotations' command (example: 'v46'). Either this "
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
        for genome_version in genome_versions:
            for tool in tools:
                tag = get_tag(tool, genome_version)
                service = get_service_name(tool, genome_version)
                concurrency = 6    # if genome_version == '37' else 2
                min_instances = 0  # if tool == 'pangolin' else 2
                max_instances = 3
                if not args.command or args.command == "build":
                    run(f"{args.docker_command} build -f docker/{tool}/Dockerfile --build-arg=\"CONCURRENCY={concurrency}\" --build-arg=\"GENOME_VERSION={genome_version}\" -t {tag}:latest .")
                    run(f"{args.docker_command} push {tag}:latest 2>&1 | grep digest: | cut -d ' ' -f 3 > docker/{tool}/sha256.txt")

                if not args.command or args.command == "deploy":
                    with open(f"docker/{tool}/sha256.txt") as f:
                        sha256 = f.read().strip()

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
