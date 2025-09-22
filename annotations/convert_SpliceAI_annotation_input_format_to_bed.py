"""
The original SpliceAI annotation format is:

#NAME   CHROM   STRAND  TX_START        TX_END  EXON_START      EXON_END
OR4F5   1       +       69090   70008   69090,  70008,
OR4F16  1       -       685715  686654  685715, 686654,
...

Convert it to BED format so that it can be viewed in IGV.
"""


import argparse
import gzip
import os
import re

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", type=int, help="Number of lines to process (for testing)")
    parser.add_argument("-o", "--output-prefix", help="Output prefix for genePred file")
    parser.add_argument("spliceai_annotation_table")
    args = parser.parse_args()

    if not args.output_prefix:
        args.output_prefix = re.sub("(.tsv|.txt)(.gz)?$", "", os.path.basename(args.spliceai_annotation_table))

    output_path = f"{args.output_prefix}.bed"

    line_count = 0
    fopen = gzip.open if args.spliceai_annotation_table.endswith("gz") else open
    with fopen(args.spliceai_annotation_table, "rt") as f:
        with open(output_path, "wt") as out:
            header = f.readline().strip().split("\t")
            if header != ["#NAME", "CHROM", "STRAND", "TX_START", "TX_END", "EXON_START", "EXON_END"]:
                raise ValueError(f"Unexpected header: {header}")

            for i, line in enumerate(f):
                line_count += 1
                fields = line.strip().split("\t")
                if len(fields) != 7:
                    raise ValueError(f"Expected 7 fields, got {len(fields)}: {fields}")

                name, chrom, strand, tx_start, tx_end, exon_starts, exon_ends = fields


                exon_starts = exon_starts.strip(",").split(",")
                exon_ends = exon_ends.strip(",").split(",")
                if exon_starts.count(",") != exon_ends.count(","):
                    raise ValueError(f"Mismatch in the number of exon starts and ends: {fields}")

                exon_sizes = [str(int(end) - int(start)) for start, end in zip(exon_starts, exon_ends)]
                exon_starts = [str(int(start) - int(tx_start)) for start in exon_starts]

                score = item_rgb = "."
                out.write("\t".join([
                    chrom, tx_start, tx_end, name, score, strand,
                    tx_start, tx_end, item_rgb,
                    str(len(exon_sizes)), ",".join(exon_sizes), ",".join(exon_starts),
                ]) + "\n")

                if args.n is not None and i > args.n:
                    break

    os.system(f"bgzip -f {output_path}")
    os.system(f"tabix -f {output_path}.gz")

    print(f"Wrote {line_count:,d} lines to {output_path}.gz")

if __name__ == "__main__":
    main()
