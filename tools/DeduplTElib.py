#!/usr/bin/env python3

import argparse
import hashlib
import re
from collections import OrderedDict

from Bio import SeqIO


def dedup_fasta(fasta_in, fasta_out, mapping_out):
    """
    Deduplicate FASTA by exact sequence identity.
    Keep the first encountered ID as canonical.
    """

    hash_to_canonical = OrderedDict()
    id_map = {}

    n_total = 0
    n_unique = 0

    with open(fasta_out, "w") as fout:

        for rec in SeqIO.parse(fasta_in, "fasta"):

            n_total += 1

            te_id = rec.id.split("#")[0]
            seq = str(rec.seq).upper()

            seq_hash = hashlib.md5(seq.encode()).hexdigest()

            if seq_hash not in hash_to_canonical:

                hash_to_canonical[seq_hash] = te_id

                SeqIO.write(rec, fout, "fasta")

                n_unique += 1

            id_map[te_id] = hash_to_canonical[seq_hash]

    with open(mapping_out, "w") as out:

        out.write("old_id\tcanonical_id\n")

        for old_id in sorted(id_map):
            out.write(f"{old_id}\t{id_map[old_id]}\n")

    print(f"[INFO] FASTA records: {n_total}")
    print(f"[INFO] Unique sequences: {n_unique}")
    print(f"[INFO] Duplicates removed: {n_total - n_unique}")

    return id_map


def remap_gtf(gtf_in, gtf_out, id_map):
    """
    Replace gene_id with canonical ID.
    Preserve coordinates and all records.
    Add original_gene_id attribute.
    """

    gene_id_re = re.compile(r'gene_id "([^"]+)"')

    n_records = 0
    n_remapped = 0

    with open(gtf_in) as fin, open(gtf_out, "w") as fout:

        for line in fin:

            if line.startswith("#"):
                fout.write(line)
                continue

            n_records += 1

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                fout.write(line)
                continue

            attrs = fields[8]

            m = gene_id_re.search(attrs)

            if not m:
                fout.write(line)
                continue

            old_id = m.group(1)
            canonical_id = id_map.get(old_id, old_id)

            if canonical_id != old_id:
                n_remapped += 1

            attrs = gene_id_re.sub(
                f'gene_id "{canonical_id}"',
                attrs,
                count=1
            )

            if 'original_gene_id "' not in attrs:

                attrs = attrs.rstrip().rstrip(";")
                attrs += f'; original_gene_id "{old_id}";'

            fields[8] = attrs

            fout.write("\t".join(fields) + "\n")

    print(f"[INFO] GTF records processed: {n_records}")
    print(f"[INFO] Records remapped: {n_remapped}")


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Deduplicate TE FASTA by exact sequence identity "
            "and remap gene_id in GTF."
        )
    )

    parser.add_argument(
        "-f",
        "--fasta",
        required=True,
        help="Input FASTA file"
    )

    parser.add_argument(
        "-g",
        "--gtf",
        required=True,
        help="Input GTF file"
    )

    parser.add_argument(
        "-o",
        "--out-prefix",
        required=True,
        help="Output prefix"
    )

    args = parser.parse_args()

    fasta_out = f"{args.out_prefix}.dedup.fa"
    gtf_out = f"{args.out_prefix}.remapped.gtf"
    mapping_out = f"{args.out_prefix}.id_mapping.tsv"

    print("[INFO] Deduplicating FASTA...")

    id_map = dedup_fasta(
        args.fasta,
        fasta_out,
        mapping_out
    )

    print("[INFO] Remapping GTF...")

    remap_gtf(
        args.gtf,
        gtf_out,
        id_map
    )

    print("[INFO] Done.")
    print(f"[INFO] FASTA:   {fasta_out}")
    print(f"[INFO] GTF:     {gtf_out}")
    print(f"[INFO] Mapping: {mapping_out}")


if __name__ == "__main__":
    main()
