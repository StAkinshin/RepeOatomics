#!/usr/bin/env python3

import argparse
import os
import re
from collections import defaultdict


def parse_gtf_attributes(attr_string):
    """
    Разбор атрибутов GTF:
    gene_id "xxx"; class_id "DNA";
    """
    attrs = {}

    for item in attr_string.strip().split(";"):
        item = item.strip()

        if not item:
            continue

        m = re.match(r'(\S+)\s+"([^"]+)"', item)

        if m:
            attrs[m.group(1)] = m.group(2)

    return attrs


def parse_gtf_and_bin(gtf_file, bin_size, out_dir, metric="density"):

    # gene_id -> информация о TE
    te_records = {}

    chrom_sizes = defaultdict(int)

    print(f"Чтение файла {gtf_file}...")

    with open(gtf_file, "r", encoding="utf-8") as f:

        for line in f:

            if line.startswith("#") or not line.strip():
                continue

            parts = line.rstrip().split("\t")

            if len(parts) < 9:
                continue

            seqid = parts[0]

            start = int(parts[3]) - 1
            end = int(parts[4])

            attrs = parse_gtf_attributes(parts[8])

            gene_id = attrs.get("gene_id")

            if gene_id is None:
                continue

            te_type = attrs.get("class_id", "Unknown")
            te_type_safe = (
                te_type.replace("/", "_")
                .replace("?", "unknown")
                .replace(" ", "_")
            )

            if gene_id not in te_records:

                te_records[gene_id] = {
                    "seqid": seqid,
                    "start": start,
                    "end": end,
                    "type": te_type_safe
                }

            else:
                te_records[gene_id]["start"] = min(
                    te_records[gene_id]["start"],
                    start
                )

                te_records[gene_id]["end"] = max(
                    te_records[gene_id]["end"],
                    end
                )

            if end > chrom_sizes[seqid]:
                chrom_sizes[seqid] = end

    # Тип -> хромосома -> интервалы
    te_data = defaultdict(lambda: defaultdict(list))

    for gene_id, info in te_records.items():

        te_data[info["type"]][info["seqid"]].append(
            (info["start"], info["end"])
        )

    os.makedirs(out_dir, exist_ok=True)

    print(
        f"Найдено {len(te_data)} классов TE. "
        f"Размер бина: {bin_size} bp"
    )

    for te_type, seqs in te_data.items():

        outfile = os.path.join(out_dir, f"{te_type}.txt")

        with open(outfile, "w", encoding="utf-8") as out:

            for seqid, chrom_len in chrom_sizes.items():

                num_bins = (
                    chrom_len // bin_size
                    + (1 if chrom_len % bin_size else 0)
                )

                bin_coverage = [0] * num_bins
                bin_counts = [0] * num_bins

                for start, end in seqs.get(seqid, []):

                    start_bin = start // bin_size
                    end_bin = (end - 1) // bin_size

                    for b in range(start_bin, end_bin + 1):

                        bin_start = b * bin_size
                        bin_end = min(
                            (b + 1) * bin_size,
                            chrom_len
                        )

                        overlap_start = max(start, bin_start)
                        overlap_end = min(end, bin_end)

                        overlap = max(
                            0,
                            overlap_end - overlap_start
                        )

                        if overlap > 0:

                            bin_coverage[b] += overlap
                            bin_counts[b] += 1

                for b in range(num_bins):

                    bin_start = b * bin_size
                    bin_end = min(
                        (b + 1) * bin_size,
                        chrom_len
                    )

                    actual_bin_size = (
                        bin_end - bin_start
                    )

                    if actual_bin_size <= 0:
                        continue

                    if metric == "density":

                        value = round(
                            bin_coverage[b] /
                            actual_bin_size,
                            5
                        )

                    elif metric == "count":

                        value = bin_counts[b]

                    elif metric == "bp":

                        value = bin_coverage[b]

                    out.write(
                        f"{seqid} "
                        f"{bin_start} "
                        f"{bin_end} "
                        f"{value}\n"
                    )

        print(f"  -> Создан файл: {outfile}")

    print("Готово!")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Подготовка TE из GTF для Circos"
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Входной GTF файл"
    )

    parser.add_argument(
        "-b",
        "--bin_size",
        type=int,
        default=100000,
        help="Размер бина (bp)"
    )

    parser.add_argument(
        "-o",
        "--out_dir",
        default="circos_tracks",
        help="Каталог результатов"
    )

    parser.add_argument(
        "-m",
        "--metric",
        choices=["density", "count", "bp"],
        default="count"
    )

    args = parser.parse_args()

    parse_gtf_and_bin(
        args.input,
        args.bin_size,
        args.out_dir,
        args.metric
    )
