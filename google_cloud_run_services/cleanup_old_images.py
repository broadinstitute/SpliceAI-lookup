#!/usr/bin/env python3

"""Delete old images from one Artifact Registry package, keeping the most recent few.

Every deploy pushes a new image, and nothing used to remove the old ones, so the repo had grown
to 217 images / 278 GB by 2026-09-01. This is run as the last step of each deploy job (see
.github/workflows/deploy-on-tag.yml) so the repo stays small on its own.

What it keeps:

  - the newest --keep images in the package, so a recent deploy can still be rolled back, and
  - every image a Cloud Run revision in some service's traffic block points at, whatever its age.

The second rule is the one that matters. Production and the no-traffic "dev" revision can sit on
very different images: on 2026-09-01 liftover served production from a 2026-05-12 image while its
dev revision ran one built that afternoon. A plain "newest N" rule would have deleted the image
production was serving after a few more dev-only tags, leaving a service that cannot start a new
container. Pinning the traffic block covers both the production and the dev revision of every
service, since `gcloud run deploy --tag dev --no-traffic` leaves the dev revision listed there at
0 percent.

Deleting an image only breaks revisions that are not in any traffic block, i.e. superseded ones
that would have to be rolled back to explicitly.

Examples:
  ./cleanup_old_images.py --package spliceai-38
  ./cleanup_old_images.py --package pangolin-37 --keep 5
  ./cleanup_old_images.py --package liftover --dry-run
"""

import argparse
import json
import subprocess
import sys

GCLOUD_PROJECT = "spliceai-lookup-412920"
REGION = "us-central1"
REPO = f"{REGION}-docker.pkg.dev/{GCLOUD_PROJECT}/docker"


def gcloud_json(args):
    """Run a gcloud command with --format=json and return the parsed output.

    Args:
        args (list): gcloud arguments, without the leading "gcloud" or the project/region flags.

    Returns:
        The parsed JSON, or None if the command failed.
    """
    proc = subprocess.run(
        ["gcloud"] + args + [f"--project={GCLOUD_PROJECT}", "--format=json"],
        capture_output=True, text=True)
    if proc.returncode != 0:
        print(f"  WARNING: `gcloud {' '.join(args)}` failed: {proc.stderr.strip().splitlines()[-1:]}")
        return None
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError:
        print(f"  WARNING: `gcloud {' '.join(args)}` returned output that isn't JSON")
        return None


def pinned_digests():
    """Return the set of image digests that some Cloud Run revision is currently serving.

    Covers every service in the region rather than just the ones built from the package being
    cleaned, so that a service sharing an image (spliceai-38 and spliceai-38-comprehensive run
    the same one) pins it too. An image in here is never deleted, however old it is.

    Returns:
        set: "sha256:..." digest strings. Empty if the services can't be listed, which the
            caller treats as a reason to do nothing rather than to delete unpinned images.
    """
    services = gcloud_json(["run", "services", "list", f"--region={REGION}"])
    if not services:
        return set()

    revision_names = set()
    for service in services:
        for entry in service.get("status", {}).get("traffic", []):
            if entry.get("revisionName"):
                revision_names.add(entry["revisionName"])

    digests = set()
    for name in sorted(revision_names):
        revision = gcloud_json(["run", "revisions", "describe", name, f"--region={REGION}"])
        if not revision:
            continue
        for container in revision.get("spec", {}).get("containers", []):
            if "@sha256:" in container.get("image", ""):
                digests.add(container["image"].split("@", 1)[1])
    return digests


def images_to_delete(images, pinned, keep):
    """Pick the images to delete: everything except the newest `keep` and the pinned ones.

    Tagged images are kept too. Nothing here tags anything but "latest", which is always on a
    recent image anyway, but an image someone tagged by hand is one someone wanted to keep.

    Args:
        images (list): image dicts as `gcloud artifacts docker images list` returns them.
        pinned (set): digests that a Cloud Run revision is serving, from pinned_digests().
        keep (int): how many of the newest images to keep regardless.

    Returns:
        list: the image dicts to delete, oldest first.
    """
    newest_first = sorted(images, key=lambda i: i["createTime"], reverse=True)
    keepers = {i["version"] for i in newest_first[:keep]}
    keepers |= {i["version"] for i in newest_first if i["version"] in pinned or i.get("tags")}
    return [i for i in reversed(newest_first) if i["version"] not in keepers]


def delete_images(package, images):
    """Delete the given images by digest. Returns the number that could not be deleted.

    Failures are retried once, because a multi-arch image's child manifests can only be deleted
    after the parent manifest that references them, and nothing in the order they are deleted in
    puts the parent first: the two are pushed together and share a timestamp.
    """
    remaining = list(images)
    for attempt in (1, 2):
        if attempt == 2:
            if not remaining:
                break
            print(f"  retrying {len(remaining)} image(s) that failed the first pass")
        still_failing = []
        for image in remaining:
            reference = f"{REPO}/{package}@{image['version']}"
            proc = subprocess.run(
                ["gcloud", "artifacts", "docker", "images", "delete", reference,
                 f"--project={GCLOUD_PROJECT}", "--quiet"],
                capture_output=True, text=True)
            if proc.returncode != 0:
                still_failing.append(image)
                if attempt == 2:
                    error = proc.stderr.strip().splitlines()
                    print(f"  FAILED {image['version'][:19]}: {error[-1] if error else 'unknown error'}")
        remaining = still_failing
    return len(remaining)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-p", "--package", required=True,
                        help=f"Image package to clean up within {REPO}, eg. spliceai-38")
    parser.add_argument("-k", "--keep", type=int, default=3,
                        help="How many of the newest images to keep, on top of every image a "
                             "Cloud Run revision is serving")
    parser.add_argument("-n", "--dry-run", action="store_true",
                        help="Print what would be deleted without deleting it")
    args = parser.parse_args()

    if args.keep < 1:
        parser.error("--keep must be at least 1")

    images = gcloud_json(["artifacts", "docker", "images", "list", f"{REPO}/{args.package}",
                          "--include-tags"])
    if images is None:
        print("ERROR: couldn't list the images, so nothing was deleted")
        return 1
    print(f"{REPO}/{args.package}: {len(images)} images")

    # An empty set means the Cloud Run API didn't answer, not that nothing is deployed. Deleting
    # on the strength of that would take the newest-N rule as the only protection and could
    # delete an image production is serving, so stop instead.
    pinned = pinned_digests()
    if not pinned:
        print("ERROR: couldn't read what Cloud Run is serving, so nothing was deleted")
        return 1
    print(f"{len(pinned)} image(s) are being served by a Cloud Run revision and will be kept")

    doomed = images_to_delete(images, pinned, args.keep)
    if not doomed:
        print("nothing to delete")
        return 0

    for image in doomed:
        print(f"  {'would delete' if args.dry_run else 'deleting'} "
              f"{image['createTime'][:10]}  {image['version'][:19]}")
    if args.dry_run:
        print(f"dry run: {len(doomed)} image(s) would be deleted, "
              f"{len(images) - len(doomed)} kept")
        return 0

    failures = delete_images(args.package, doomed)
    print(f"deleted {len(doomed) - failures} of {len(doomed)} image(s), "
          f"{len(images) - len(doomed)} kept")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
