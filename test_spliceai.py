import json
import os
import re
from flask import Flask, request, Response
from flask_cors import CORS
from spliceai.utils import Annotator, get_delta_scores
from server import SPLICEAI_ANNOTATOR, SPLICEAI_DEFAULT_DISTANCE, SPLICEAI_DEFAULT_MASK, VariantRecord, parse_variant


print(get_delta_scores(VariantRecord(*parse_variant("8:140300615 C>G")), SPLICEAI_ANNOTATOR["38"], SPLICEAI_DEFAULT_DISTANCE, SPLICEAI_DEFAULT_MASK))
print(get_delta_scores(VariantRecord(*parse_variant("2-179531962-C-A")), SPLICEAI_ANNOTATOR["37"], SPLICEAI_DEFAULT_DISTANCE, SPLICEAI_DEFAULT_MASK))
print(get_delta_scores(VariantRecord(*parse_variant("2-179532167-A-G")), SPLICEAI_ANNOTATOR["37"], SPLICEAI_DEFAULT_DISTANCE, SPLICEAI_DEFAULT_MASK))
print(get_delta_scores(VariantRecord(*parse_variant("2-179529170-GACAGTTAAGAATGTACCTTTGACAGGTACA-G")), SPLICEAI_ANNOTATOR["37"], SPLICEAI_DEFAULT_DISTANCE, SPLICEAI_DEFAULT_MASK))

