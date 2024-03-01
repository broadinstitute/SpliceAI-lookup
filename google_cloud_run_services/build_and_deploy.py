import argparse
import json
import logging
import os

logging.basicConfig(level=logging.INFO, format="%(asctime)s: %(message)s")

GCLOUD_PROJECT = "spliceai-lookup-412920"
CONCURRENCY = 2  # max requests per instance

def get_service_name(tool, genome_version):
	return f"{tool}-{genome_version}"

def get_tag(tool, genome_version):
	return f"us-central1-docker.pkg.dev/spliceai-lookup-412920/docker/{get_service_name(tool, genome_version)}:latest"

def run(c):
	logging.info(c)
	os.system(c)

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-g", "--genome-version", choices=["37", "38"], help="If not specified, command will run for both GRCh37 and GRCh38")
	parser.add_argument("-t", "--tool", choices=["spliceai", "pangolin"], help="If not specified, command will run for both spliceai and pangolin")
	parser.add_argument("command", nargs="?", choices=["build", "deploy", "test", "test2", "run"], help="Command to run. If not specified, it will run 'build' and then 'deploy'")

	args = parser.parse_args()

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
			print(f"docker run -it {tag} /bin/bash")
		elif args.command == "test":
			run(f"docker run -p 8080:8080 {tag}")

		return

	if args.genome_version:
		genome_versions = [args.genome_version]
	else:
		genome_versions = ["38", "37"]

	if args.tool:
		tools = [args.tool]
	else:
		tools = ["spliceai", "pangolin"]

	with open("params.json") as f:
		params = json.load(f)

	for genome_version in genome_versions:
		for tool in tools:
			tag = get_tag(tool, genome_version)
			service = get_service_name(tool, genome_version)

			if not args.command or args.command == "build":
				run(f"docker build -f docker/{tool}/Dockerfile --build-arg=\"GENOME_VERSION=GRCh{genome_version}\" -t {tag} .")
				run(f"docker push {tag}")
				run(f"docker pull {tag} 2>&1 | grep Digest | cut -c 9- > docker/{tool}/sha256.txt")

			if not args.command or args.command == "deploy":
				#f"gcloud --project $(GCLOUD_PROJECT) beta run deploy {service} --image {tag} --concurrency 5 --execution-environment gen2 --region us-central1 --set-env-vars \"DB_PASSWORD=${SPLICEAI_LOOKUP_DB_PASSWORD}\" --add-volume=name=ref,type=cloud-storage,bucket=spliceai-lookup-reference-data,readonly=true --add-volume-mount=volume=ref,mount-path=/ref --allow-unauthenticated --memory 8Gi --cpu 4"
				run(f"""gcloud \
--project {GCLOUD_PROJECT} beta run deploy {service} \
--image {tag} \
--max-instances 8 \
--concurrency {CONCURRENCY} \
--execution-environment gen2 \
--region us-central1 \
--set-env-vars "DB_PASSWORD={params['SPLICEAI_LOOKUP_DB_PASSWORD']}" \
--add-volume=name=ref,type=cloud-storage,bucket=spliceai-lookup-reference-data,readonly=true \
--allow-unauthenticated \
--memory 8Gi \
--cpu 4
""")

				# 					--add-volume-mount=volume=ref,mount-path=/ref \

if __name__ == "__main__":
	main()