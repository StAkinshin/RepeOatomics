#!/usr/bin/env python3

import argparse

def add_prefix_to_fasta_headers(input_fasta, output_fasta, prefix):
    with open(input_fasta) as fin, open(output_fasta, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                header = line[1:].strip()

                # TE_00000003#unknown -> PREFIX_TE_00000003#unknown
                if "#" in header:
                    seq_id, rest = header.split("#", 1)
                    new_header = f">{prefix}{seq_id}#{rest}\n"
                else:
                    new_header = f">{prefix}{header}\n"

                fout.write(new_header)
            else:
                fout.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta")
    parser.add_argument("output_fasta")
    parser.add_argument("--prefix", required=True)

    args = parser.parse_args()

    add_prefix_to_fasta_headers(
        args.input_fasta,
        args.output_fasta,
        args.prefix
    )
