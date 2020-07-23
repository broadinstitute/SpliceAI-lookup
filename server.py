import json
import markdown2
import os
import re
import sys
from datetime import datetime
from flask import Flask, request, Response
from flask_cors import CORS
from spliceai.utils import Annotator, get_delta_scores

app = Flask(__name__, template_folder='.')
CORS(app)

SPLICEAI_ANNOTATOR = {
    "37": Annotator(os.path.expanduser("~/hg19.fa"), "grch37"),
    "38": Annotator(os.path.expanduser("~/hg38.fa"), "grch38"),
}

SPLICEAI_MAX_DISTANCE_LIMIT = 20000
SPLICEAI_DEFAULT_DISTANCE = 50  # maximum distance between the variant and gained/lost splice site, defaults to 50
SPLICEAI_DEFAULT_MASK = 0  # mask scores representing annotated acceptor/donor gain and unannotated acceptor/donor loss, defaults to 0

SPLICEAI_SCORE_FIELDS = "ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL".split("|")

SPLICEAI_EXAMPLE = f"/spliceai/?hg=38&distance=50&mask=0&variant=chr8-140300615-C-G"

VARIANT_RE = re.compile(
    "(chr)?(?P<chrom>[0-9XYMTt]{1,2})"
    "[-\s:]+"
    "(?P<pos>[0-9]{1,9})"
    "[-\s:]+"
    "(?P<ref>[ACGT]+)"
    "[-\s:>]+"
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


def process_variant(variant, genome_version, spliceai_distance, spliceai_mask):
    try:
        chrom, pos, ref, alt = parse_variant(variant)
    except ValueError as e:
        return {
            "variant": variant,
            "error": f"ERROR: {e}",
        }

    record = VariantRecord(chrom, pos, ref, alt)
    try:
        scores = get_delta_scores(
            record,
            SPLICEAI_ANNOTATOR[genome_version],
            spliceai_distance,
            spliceai_mask)
    except Exception as e:
        return {
            "variant": variant,
            "error": f"ERROR: {type(e)}: {e}",
        }

    if len(scores) == 0:
        return {
            "variant": variant,
            "error": f"ERROR: Unable to compute scores for {variant}. Please check that the genome version and reference allele are correct, and the variant is either exonic or intronic in Gencode v24.",
        }

    scores = [s[s.index("|")+1:] for s in scores]  # drop allele field

    return {
        "variant": variant,
        "genome_version": genome_version,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "scores": scores,
    }


@app.route("/spliceai/", methods=['POST', 'GET'])
def run_spliceai():

    # check params
    params = {}
    if request.values:
        params.update(request.values)

    if 'variant' not in params:
        params.update(request.get_json(force=True, silent=True) or {})

    variant = params.get('variant', '')
    variant = variant.strip().strip("'").strip('"').strip(",")
    if not variant:
        return f'"variant" not specified. For example: {SPLICEAI_EXAMPLE}\n', 400

    if not isinstance(variant, str):
        return f'"variant" value must be a string rather than a {type(variant)}.\n', 400

    genome_version = params.get("hg")
    if not genome_version:
        return f'"hg" not specified. The URL must include an "hg" arg: hg=37 or hg=38. For example: {SPLICEAI_EXAMPLE}\n', 400

    if genome_version not in ("37", "38"):
        return f'Invalid "hg" value: "{genome_version}". The value must be either "37" or "38". For example: {SPLICEAI_EXAMPLE}\n', 400

    spliceai_distance = params.get("distance", SPLICEAI_DEFAULT_DISTANCE)
    try:
        spliceai_distance = int(spliceai_distance)
    except Exception as e:
        return f'Invalid "distance": "{spliceai_distance}". The value must be an integer.\n', 400

    if spliceai_distance > SPLICEAI_MAX_DISTANCE_LIMIT:
        return f'Invalid "distance": "{spliceai_distance}". The value must be < {SPLICEAI_MAX_DISTANCE_LIMIT}.\n', 400

    spliceai_mask = params.get("mask", str(SPLICEAI_DEFAULT_MASK))
    if spliceai_mask not in ("0", "1"):
        return f'Invalid "mask" value: "{spliceai_mask}". The value must be either "0" or "1". For example: {SPLICEAI_EXAMPLE}\n', 400

    start_time = datetime.now()
    print(f"======================", flush=True)
    print(start_time.strftime("%m/%d/%Y %H:%M:%S"), flush=True)
    print(f"{request.remote_addr}", flush=True)
    print(f"Processing {variant}  with hg={genome_version}, distance={spliceai_distance}, mask={spliceai_mask}", flush=True)

    results = process_variant(variant, genome_version, spliceai_distance, spliceai_mask)

    print(f"Done processing variant: {variant}", flush=True)
    print(f"Results: {results}", flush=True)
    print(f"This took " + str(datetime.now() - start_time), flush=True)

    status = 400 if results.get("error") else 200
    return Response(json.dumps(results), status=status, mimetype='application/json')


@app.route('/', defaults={'path': ''})
@app.route('/<path:path>/')
def catch_all(path):
    with open("README.md") as f:
        return markdown2.markdown(f.read())


print("Initialization completed.", flush=True)

if __name__ == "__main__":
    app.run(debug=False, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))
