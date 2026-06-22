#!/usr/bin/env python3

import re
import argparse

def add_prefix_to_gtf(input_gtf, output_gtf, prefix):
    gene_pattern = re.compile(r'gene_id "([^"]+)"')
    transcript_pattern = re.compile(r'transcript_id "([^"]+)"')

    with open(input_gtf, "r") as fin, open(output_gtf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            line = gene_pattern.sub(
                lambda m: f'gene_id "{prefix}{m.group(1)}"',
                line
            )

            line = transcript_pattern.sub(
                lambda m: f'transcript_id "{prefix}{m.group(1)}"',
                line
            )

            fout.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add prefix to gene_id and transcript_id in a GTF file"
    )
    parser.add_argument("input_gtf", help="Input GTF file")
    parser.add_argument("output_gtf", help="Output GTF file")
    parser.add_argument("--prefix", required=True,
                        help="Prefix to add")

    args = parser.parse_args()

    add_prefix_to_gtf(args.input_gtf, args.output_gtf, args.prefix)
