import json
import os
import re
from flask import Flask, request, Response
from flask_cors import CORS
from spliceai.utils import Annotator, get_delta_scores

SPLICEAI_ANNOTATOR = {
    "37": Annotator(os.path.expanduser("hg19.fa"), "grch37"),
    "38": Annotator(os.path.expanduser("hg38.fa"), "grch38"),
}

SPLICEAI_DEFAULT_DISTANCE = 50  # maximum distance between the variant and gained/lost splice site, defaults to 50
SPLICEAI_DEFAULT_MASK = 0  # mask scores representing annotated acceptor/donor gain and unannotated acceptor/donor loss, defaults to 0

SPLICEAI_SCORE_FIELDS = "ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL".split("|")

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


SPLICEAI_EXAMPLE = f"/spliceai/?hg=38&variants=chr8-140300615-C-G"

@app.route("/spliceai/?", methods=['POST', 'GET'])
def get_spliceai_scores():

    # check params
    params = {}
    if request.values:
        params.update(request.values)

    if 'variants' not in params:
        params.update(request.get_json(force=True, silent=True) or {})

    genome_version = params.get("hg")
    if not genome_version:
        return f'"hg" not specified. The URL must include an "hg" arg: hg=37 or hg=38. For example: {SPLICEAI_EXAMPLE}\n', 400

    if genome_version not in ("37", "38"):
        return f'Invalid "hg" value: "{genome_version}". The value must be either "37" or "38". For example: {SPLICEAI_EXAMPLE}\n', 400

    variants = params.get('variants')
    if not variants:
        return f'"variants" not specified. The URL must include a "variants" arg. For example: {SPLICEAI_EXAMPLE}\n', 400

    if isinstance(variants, str):
        variants = variants.split(",")
    else:
        return f'"variants" value must be a string rather than a {type(variants)}.\n', 400

    spliceai_distance = params.get("distance", SPLICEAI_DEFAULT_DISTANCE)
    try:
        spliceai_distance = int(spliceai_distance)
    except Exception as e:
        return f'Invalid "distance" value: "{spliceai_distance}". The value must be an integer.\n', 400

    # parse and perform liftover
    results = []
    for variant in variants:
        variant = variant.strip().strip("'").strip('"').strip(",")
        if not variant:
            continue

        try:
            chrom, pos, ref, alt = parse_variant(variant)
        except ValueError as e:
            results.append({"variant": variant, "error": str(e)})
            continue

        record = VariantRecord(chrom, pos, ref, alt)
        try:
            scores = get_delta_scores(
                record,
                SPLICEAI_ANNOTATOR[genome_version],
                spliceai_distance,
                SPLICEAI_DEFAULT_MASK)
        except Exception as e:
            results.append({"variant": variant, "error": f"{type(e)}: {e}"})
            continue

        if len(scores) == 0:
            results.append({"variant": variant, "error": f"unable to compute scores for {variant}"})
            continue

        parsed_scores = []
        for score in scores:
            parsed_scores.append(dict(zip(SPLICEAI_SCORE_FIELDS, score.split("|"))))
        results.append({"variant": variant, "scores": parsed_scores})

    return Response(json.dumps(results),  mimetype='application/json')


@app.route('/?', defaults={'path': ''})
@app.route('/<path:path>/?')
def catch_all(path):

    return f"""<html>
<head>
<title>TGG APIs</title>
</head>
<body style="font-family: monospace">
This server provides the following APIs:<br/>
<br />

GET {SPLICEAI_EXAMPLE} <br />
<b>hg</b> can be: {' or '.join(SPLICEAI_ANNOTATOR.keys())} <br />
<b>variants</b> can have the format "chrom:pos ref&gt;alt" or "chrom-pos-ref-alt" or "chrom pos ref alt" <br />
<br />
<b>spliceai API response:</b> a json list that's the same length as the input "variants" list and has the splice AI scores or error message for each variant.<br/>
<br />
</body>
</html>"""


if __name__ == "__main__":
    app.run(debug=False, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))
