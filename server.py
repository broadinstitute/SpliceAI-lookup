import json
import markdown2
import os
import re
import sys
from datetime import datetime
from flask import Flask, request, Response
from flask_cors import CORS
from spliceai.utils import Annotator, get_delta_scores

SPLICEAI_ANNOTATOR = {
    "37": Annotator(os.path.expanduser("~/hg19.fa"), "grch37"),
    "38": Annotator(os.path.expanduser("~/hg38.fa"), "grch38"),
}

SPLICEAI_MAX_INPUT_VARIANTS = 100
SPLICEAI_MAX_DISTANCE_LIMIT = 20000
SPLICEAI_DEFAULT_DISTANCE = 50  # maximum distance between the variant and gained/lost splice site, defaults to 50
SPLICEAI_DEFAULT_MASK = 0  # mask scores representing annotated acceptor/donor gain and unannotated acceptor/donor loss, defaults to 0

SPLICEAI_SCORE_FIELDS = "ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL".split("|")

app = Flask(__name__, template_folder='.')
CORS(app)

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


def get_ucsc_link(genome_version, chrom, pos):
    genome_version = genome_version.replace('37', '19')
    chrom = chrom.replace('chr', '')
    return f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg{genome_version}&position=chr{chrom}:{pos}"


class VariantRecord:
    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alts = [alt]

    def __repr__(self):
        return f"{self.chrom}-{self.pos}-{self.ref}-{self.alts[0]}"


SPLICEAI_EXAMPLE = f"/spliceai/?hg=38&variants=chr8-140300615-C-G"

@app.route("/spliceai/", methods=['POST', 'GET'])
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

    if len(variants) > SPLICEAI_MAX_INPUT_VARIANTS:
        return f'{SPLICEAI_MAX_INPUT_VARIANTS} variant limit exceeded. The provided list contains {len(variants)} variants.\n', 400


    spliceai_distance = params.get("distance", SPLICEAI_DEFAULT_DISTANCE)
    try:
        spliceai_distance = int(spliceai_distance)
    except Exception as e:
        return f'Invalid "distance": "{spliceai_distance}". The value must be an integer.\n', 400

    if spliceai_distance > SPLICEAI_MAX_DISTANCE_LIMIT:
        return f'Invalid "distance": "{spliceai_distance}". The value must be < {SPLICEAI_MAX_DISTANCE_LIMIT}.\n', 400

    # parse and perform liftover
    results = []
    for variant in variants:
        variant = variant.strip().strip("'").strip('"').strip(",")
        if not variant:
            continue

        start_time = datetime.now()
        print(start_time.strftime("%d/%m/%Y %H:%M:%S ") + f"Processing {variant}  with hg={genome_version}, distance={spliceai_distance}", file=sys.stderr)
        try:
            chrom, pos, ref, alt = parse_variant(variant)
        except ValueError as e:
            results.append({
                "variant": variant,
                "url": "#",
                "error": str(e),
            })
            continue

        record = VariantRecord(chrom, pos, ref, alt)
        try:
            scores = get_delta_scores(
                record,
                SPLICEAI_ANNOTATOR[genome_version],
                spliceai_distance,
                SPLICEAI_DEFAULT_MASK)
        except Exception as e:
            results.append({
                "variant": variant,
                "url": get_ucsc_link(genome_version, chrom, pos),
                "error": f"{type(e)}: {e}",
            })
            continue

        if len(scores) == 0:
            results.append({
                "variant": variant,
                "url": get_ucsc_link(genome_version, chrom, pos),
                "error": f"Unable to compute scores for {variant}. Please check that the genome version and reference allele are correct, and the variant is either exonic or intronic.",
            })
            continue

        parsed_scores = []
        for score_fields in scores:
            score_dict = dict(zip(SPLICEAI_SCORE_FIELDS, score_fields.split("|")))
            for score_type in "DG", "DL", "AG", "AL":
                try:
                    score_dict[f"{score_type}_url"] = get_ucsc_link(genome_version, chrom, pos + int(score_dict[f"DP_{score_type}"]))
                except Exception as e:
                    print("{type(e)}: {e}", file=sys.stderr)
                    score_dict[f"{score_type}_url"] = "#"

            parsed_scores.append(score_dict)

        results.append({
            "variant": variant,
            "url": get_ucsc_link(genome_version, chrom, pos),
            "scores": parsed_scores,
        })

        print(f"Done processing variant: {variant}. This took " + str(datetime.now() - start_time), file=sys.stderr)

    return Response(json.dumps(results), mimetype='application/json')


@app.route('/', defaults={'path': ''})
@app.route('/<path:path>/')
def catch_all(path):
    with open("README.md") as f:
        return markdown2.markdown(f.read())


print("Initialization completed.", file=sys.stderr)

if __name__ == "__main__":
    app.run(debug=False, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))
