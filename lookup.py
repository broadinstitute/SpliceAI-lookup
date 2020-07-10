import json
import os
import re
from flask import Flask, request, Response
from flask_cors import CORS
from spliceai.utils import Annotator, get_delta_scores

ANNOTATOR = {
    "grch37": Annotator(os.path.expanduser("~/p1/ref/GRCh37/hg19.fa"), "grch37"),
    "grch38": Annotator(os.path.expanduser("~/p1/ref/GRCh38/hg38.fa"), "grch38"),
}

DISTANCE = 50  # maximum distance between the variant and gained/lost splice site, defaults to 50
MASK = 0  # mask scores representing annotated acceptor/donor gain and unannotated acceptor/donor loss, defaults to 0

app = Flask(__name__, template_folder='.')
CORS(app)

VARIANT_RE = re.compile(
    "(chr)?(?P<chrom>[0-9XYMTt]{1,2})"
    "[: -]+"
    "(?P<pos>[0-9]{1,9})"
    "[: -]+"
    "(?P<ref>[ACGT]+)"
    "[: ->]+"
    "(?P<alt>[ACGT]+)"
)


def parse_variant(variant_str):
    match = VARIANT_RE.match(variant_str)
    if not match:
        raise ValueError(f"Unable to parse variant: {variant_str}")

    return match['chrom'], int(match['pos']), match['ref'], match['alt']


class VariantRecord:
    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alts = [alt]

    def __repr__(self):
        return f"{self.chrom}-{self.pos}-{self.ref}-{self.alts[0]}"


print(get_delta_scores(VariantRecord(*parse_variant("8:140300615 C>G")), ANNOTATOR["grch38"], DISTANCE, MASK))
print(get_delta_scores(VariantRecord(*parse_variant("2-179531962-C-A")), ANNOTATOR["grch37"], DISTANCE, MASK))
print(get_delta_scores(VariantRecord(*parse_variant("2-179532167-A-G")), ANNOTATOR["grch37"], DISTANCE, MASK))
print(get_delta_scores(VariantRecord(*parse_variant("2-179529170-GACAGTTAAGAATGTACCTTTGACAGGTACA-G")), ANNOTATOR["grch37"], DISTANCE, MASK))

