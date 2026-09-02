"""Tests for cleanup_old_images.py's keep/delete decision.

images_to_delete is the whole risk in that script: get it wrong and a deploy silently removes
the image a Cloud Run service is serving, leaving a service that can't start a container. These
tests cover the case that made the pinning necessary, drawn from what the repo actually looked
like on 2026-09-01: liftover served production from a 2026-05-12 image while three newer ones
sat above it, so the "newest N" rule alone would have deleted the live one.

Run with:  python3 -m unittest cleanup_old_images_tests -v
"""

import unittest
from unittest import mock

import cleanup_old_images
from cleanup_old_images import images_to_delete, pinned_digests


def image(digest, create_time, tags=None):
    """One entry shaped like `gcloud artifacts docker images list --include-tags` returns."""
    return {"version": f"sha256:{digest}", "createTime": create_time, "tags": tags or []}


# Newest first. The oldest is the one production is serving in the tests below.
IMAGES = [
    image("ee", "2026-09-01T18:42:15Z"),
    image("dd", "2026-09-01T18:42:13Z"),
    image("cc", "2026-09-01T18:42:03Z"),
    image("bb", "2026-08-31T15:04:09Z"),
    image("aa", "2026-05-12T14:09:53Z"),
]


def digests(images):
    return {i["version"] for i in images}


class ImagesToDeleteTest(unittest.TestCase):

    def test_keeps_the_newest_images(self):
        deleted = images_to_delete(IMAGES, pinned={"sha256:ee"}, keep=3)
        self.assertEqual(digests(deleted), {"sha256:bb", "sha256:aa"})

    def test_keeps_a_pinned_image_however_old_it_is(self):
        """The liftover case: production sits on the oldest image of the five."""
        deleted = images_to_delete(IMAGES, pinned={"sha256:aa"}, keep=3)
        self.assertEqual(digests(deleted), {"sha256:bb"})
        self.assertNotIn("sha256:aa", digests(deleted))

    def test_keeps_every_pinned_image_when_prod_and_dev_differ(self):
        """A dev-only tag leaves production and the dev revision on different images."""
        deleted = images_to_delete(IMAGES, pinned={"sha256:aa", "sha256:ee"}, keep=1)
        self.assertEqual(digests(deleted), {"sha256:dd", "sha256:cc", "sha256:bb"})

    def test_keeps_tagged_images(self):
        tagged = [image("aa", "2026-05-12T14:09:53Z", tags=["latest"]) if i["version"] == "sha256:aa"
                  else i for i in IMAGES]
        deleted = images_to_delete(tagged, pinned=set(), keep=2)
        self.assertEqual(digests(deleted), {"sha256:cc", "sha256:bb"})

    def test_deletes_nothing_when_keep_covers_everything(self):
        self.assertEqual(images_to_delete(IMAGES, pinned=set(), keep=len(IMAGES)), [])

    def test_deletes_oldest_first(self):
        deleted = images_to_delete(IMAGES, pinned=set(), keep=1)
        self.assertEqual([i["createTime"] for i in deleted],
                         sorted(i["createTime"] for i in deleted))

    def test_ignores_the_order_the_images_arrive_in(self):
        shuffled = [IMAGES[2], IMAGES[0], IMAGES[4], IMAGES[1], IMAGES[3]]
        self.assertEqual(digests(images_to_delete(shuffled, pinned={"sha256:ee"}, keep=3)),
                         {"sha256:bb", "sha256:aa"})


class PinnedDigestsTest(unittest.TestCase):
    """A revision that can't be described may be the only one serving some image, so a partial
    answer has to come back as None rather than as the digests that could be read: main() would
    otherwise take the partial set as complete and delete the unread revision's image."""

    SERVICES = [{"status": {"traffic": [{"revisionName": "spliceai-38-00001"},
                                        {"revisionName": "spliceai-38-00002"}]}}]

    @staticmethod
    def fake_gcloud_json(failing_revision):
        """gcloud_json stand-in: lists SERVICES and describes each revision as serving an image
        named after it, except failing_revision, whose describe fails."""
        def fake(args):
            if args[:3] == ["run", "services", "list"]:
                return PinnedDigestsTest.SERVICES
            name = args[3]
            if name == failing_revision:
                return None
            return {"spec": {"containers": [{"image": f"{cleanup_old_images.REPO}/spliceai-38@sha256:{name[-1]}"}]}}
        return fake

    def test_returns_every_served_digest(self):
        with mock.patch("cleanup_old_images.gcloud_json", self.fake_gcloud_json(failing_revision=None)):
            self.assertEqual(pinned_digests(), {"sha256:1", "sha256:2"})

    def test_returns_none_when_one_revision_cannot_be_described(self):
        with mock.patch("cleanup_old_images.gcloud_json", self.fake_gcloud_json(failing_revision="spliceai-38-00002")):
            self.assertIsNone(pinned_digests())

    def test_returns_none_when_services_cannot_be_listed(self):
        with mock.patch("cleanup_old_images.gcloud_json", lambda args: None):
            self.assertIsNone(pinned_digests())


if __name__ == "__main__":
    unittest.main()
