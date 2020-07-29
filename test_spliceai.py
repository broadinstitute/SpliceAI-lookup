import unittest
from server import SPLICEAI_ANNOTATOR, SPLICEAI_DEFAULT_DISTANCE, SPLICEAI_DEFAULT_MASK, VariantRecord, parse_variant, process_variant


class Test(unittest.TestCase):

    def test_parse_variant(self):
        self.assertEqual(parse_variant("chr3:12345 A>G"), ("3", 12345, "A", "G"))
        self.assertEqual(parse_variant("3:12345:A:G"), ("3", 12345, "A", "G"))
        self.assertEqual(parse_variant("chrX:12345:A:G"), ("X", 12345, "A", "G"))
        self.assertEqual(parse_variant("chrY:12345:A:G"), ("Y", 12345, "A", "G"))
        with self.assertRaises(ValueError):
            parse_variant("Z:12345:A:G")

    def test_spliceai_results(self):
        # from test_data/spliceai_scores.raw.indel.hg38_subset.vcf.gz
        # 1       69091   .       A       AA      .       .       SpliceAI=AA|OR4F5|0.00|0.00|0.03|0.00|-15|42|2|24
        # 1       69124   .       GATT    G       .       .       SpliceAI=G|OR4F5|0.00|0.02|0.00|0.06|18|9|27|-31

        variant = "1 69091 A AA"
        for distance in SPLICEAI_DEFAULT_DISTANCE, SPLICEAI_DEFAULT_DISTANCE - 1:
            result = process_variant(variant, "38", distance, SPLICEAI_DEFAULT_MASK)
            self.assertEqual(result['variant'], variant)
            self.assertEqual(result['chrom'], "1")
            self.assertEqual(result['pos'], 69091)
            self.assertEqual(result['ref'], "A")
            self.assertEqual(result['alt'], "AA")
            self.assertEqual(result['genome_version'], "38")
            self.assertEqual(result['source'], "lookup" if distance == SPLICEAI_DEFAULT_DISTANCE else "computed")

        #print(get_delta_scores(VariantRecord(*parse_variant("2-179531962-C-A")), SPLICEAI_ANNOTATOR["37"], SPLICEAI_DEFAULT_DISTANCE, SPLICEAI_DEFAULT_MASK))
        #print(get_delta_scores(VariantRecord(*parse_variant("2-179532167-A-G")), SPLICEAI_ANNOTATOR["37"], SPLICEAI_DEFAULT_DISTANCE, SPLICEAI_DEFAULT_MASK))
        #print(get_delta_scores(VariantRecord(*parse_variant("2-179529170-GACAGTTAAGAATGTACCTTTGACAGGTACA-G")), SPLICEAI_ANNOTATOR["37"], SPLICEAI_DEFAULT_DISTANCE, SPLICEAI_DEFAULT_MASK))

