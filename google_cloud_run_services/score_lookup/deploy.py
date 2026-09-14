"""Deploys the score-lookup function to Cloud Run.

The same dev/production split as build_and_deploy.py uses for the SpliceAI and Pangolin services:

    python3 deploy.py            # production: builds and deploys a revision that takes all traffic
    python3 deploy.py --dev      # builds and deploys a revision tagged 'dev' that takes no traffic,
                                 #   served at https://dev---score-lookup-<hash>-uc.a.run.app
    python3 deploy.py --promote  # moves all traffic to the revision tagged 'dev', without rebuilding

The dev site (and any page opened with api=dev) sends its lookups to the 'dev' tag, so a --dev deploy
can be tried there before --promote puts that exact revision into production.
"""

import argparse
import os
import subprocess
import sys

GCLOUD_PROJECT = "spliceai-lookup-412920"
REGION = "us-central1"
SERVICE = "score-lookup"
SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))


def run(command):
    """Prints and runs a command, exiting if it fails."""
    print(" ".join(command))
    if subprocess.run(command).returncode != 0:
        sys.exit(f"Command failed: {' '.join(command)}")


def service_exists():
    """Returns whether the Cloud Run service has been created."""
    return subprocess.run(
        ["gcloud", "run", "services", "describe", SERVICE, f"--project={GCLOUD_PROJECT}", f"--region={REGION}",
         "--format=value(metadata.name)"],
        capture_output=True).returncode == 0


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--dev", action="store_true", help="Deploy a no-traffic revision tagged 'dev'")
    mode.add_argument("--promote", action="store_true", help="Move all traffic to the revision tagged 'dev'")
    args = parser.parse_args()

    if args.promote:
        run(["gcloud", "run", "services", "update-traffic", SERVICE, f"--project={GCLOUD_PROJECT}",
             f"--region={REGION}", "--to-tags=dev=100"])
        return

    if args.dev and not service_exists():
        # --no-traffic can't be used when creating a service, so its first revision would take all
        # traffic while being the untested dev build.
        sys.exit(f"{SERVICE} doesn't exist yet. Create it with a production deploy first, then re-run with --dev.")

    run([
        "gcloud", "run", "deploy", SERVICE,
        "--quiet",  # answer yes to setup prompts, such as creating the Artifact Registry repo for source builds
        f"--project={GCLOUD_PROJECT}",
        f"--region={REGION}",
        f"--source={SOURCE_DIR}",
        "--function=scoreLookup",
        "--base-image=nodejs22",
        "--allow-unauthenticated",
        # Scale to zero when idle: a cold start costs the next request about a second, which is cheaper
        # than paying for an idle instance around the clock.
        "--min-instances=0",
        # Caps what a flood of requests can cost.
        "--max-instances=20",
        # Each request mostly waits on reads from Cloud Storage, so one instance handles many at once.
        "--concurrency=80",
        "--cpu=1",
        # Measured at about 150 MiB with the five tables' indexes loaded, and under 200 MiB after 150
        # distinct lookups, since lookup.js caps each table's chunk cache at 16 MiB.
        "--memory=512Mi",
        "--cpu-boost",
        "--timeout=60s",
    ] + (["--tag=dev", "--no-traffic"] if args.dev else []))

    if not args.dev:
        # A --dev deploy or --promote leaves traffic pinned to a named revision, and while it is,
        # `gcloud run deploy` gives the new revision no traffic, so send it all to the latest.
        run(["gcloud", "run", "services", "update-traffic", SERVICE, f"--project={GCLOUD_PROJECT}",
             f"--region={REGION}", "--to-latest"])


if __name__ == "__main__":
    main()
