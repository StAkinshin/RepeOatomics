import argparse
import os
import re
from collections import defaultdict

def parse_gff3_and_bin(gff_file, bin_size, out_dir, metric="density"):
    # Словарь: Тип_ТЭ -> Хромосома -> список координат [(start, end), ...]
    te_data = defaultdict(lambda: defaultdict(list))
    # Словарь для хранения примерной длины хромосом (по максимальной координате)
    chrom_sizes = defaultdict(int)

    print(f"Чтение файла {gff_file}...")
    with open(gff_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
                
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            seqid = parts[0]
            col3_type = parts[2]
            
            # GFF3 использует 1-based координаты включительно [start, end].
            # Переводим в 0-based полуоткрытые интервалы [start_0, end_0) для удобства расчетов и Circos.
            start_0 = int(parts[3]) - 1
            end_0 = int(parts[4])
            attributes = parts[8]

            # Извлекаем класс транспозона из атрибута classification (если есть)
            match = re.search(r'classification=([^;]+)', attributes)
            if match:
                te_type = match.group(1)
            else:
                te_type = col3_type # Фоллбэк на 3 колонку, если нет классификации

            # Делаем имя безопасным для сохранения в качестве файла (заменяем слеши)
            te_type_safe = te_type.replace('/', '_').replace('?', 'unknown').replace(' ', '_')

            te_data[te_type_safe][seqid].append((start_0, end_0))
            
            # Обновляем длину хромосомы
            if end_0 > chrom_sizes[seqid]:
                chrom_sizes[seqid] = end_0

    # Создаем выходную директорию
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    print(f"Найдено {len(te_data)} уникальных типов транспозонов. Расчет бинов (размер: {bin_size} bp)...")

    # Для каждого типа транспозона создаем отдельный файл
    for te_type, seqs in te_data.items():
        out_filename = os.path.join(out_dir, f"{te_type}.txt")
        
        with open(out_filename, 'w', encoding='utf-8') as out:
            for seqid, length in chrom_sizes.items():
                # Количество бинов для данной хромосомы
                num_bins = (length // bin_size) + (1 if length % bin_size != 0 else 0)

                # Массивы для хранения покрытия и количества
                bin_coverage = [0] * num_bins
                bin_counts = [0] * num_bins

                # Считаем пересечение транспозонов с бинами
                for (start, end) in seqs[seqid]:
                    start_bin = start // bin_size
                    end_bin = (end - 1) // bin_size 

                    for b in range(start_bin, end_bin + 1):
                        bin_start = b * bin_size
                        bin_end = min((b + 1) * bin_size, length)

                        # Находим длину перекрытия текущего TE с текущим бином
                        overlap_start = max(start, bin_start)
                        overlap_end = min(end, bin_end)
                        overlap = max(0, overlap_end - overlap_start)

                        if overlap > 0:
                            bin_coverage[b] += overlap
                            bin_counts[b] += 1

                # Записываем результат в формате для Circos: seqid start end value
                for b in range(num_bins):
                    bin_start = b * bin_size
                    bin_end = min((b + 1) * bin_size, length)
                    actual_bin_size = bin_end - bin_start

                    if actual_bin_size <= 0: continue

                    if metric == "density":
                        # Доля бина, покрытая транспозонами (0.0 - 1.0)
                        val = bin_coverage[b] / actual_bin_size
                        val = round(val, 5)
                    elif metric == "count":
                        # Количество транспозонов, попадающих в бин
                        val = bin_counts[b]
                    else: # bp
                        # Точное число пар нуклеотидов в бине
                        val = bin_coverage[b]

                    out.write(f"{seqid} {bin_start} {bin_end} {val}\n")

        print(f"  -> Создан файл: {out_filename}")
        
    print("Готово!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Подготовка данных о TE из EDTA (GFF3) для Circos")
    parser.add_argument("-i", "--input", required = True, help = "Входной GFF3 файл от EDTA")
    parser.add_argument("-b", "--bin_size", type = int, default = 100000, help = "Размер бина (в парах нуклеотидов). По умолчанию: 100000 (100 kb)")
    parser.add_argument("-o", "--out_dir", default = "circos_tracks", help = "Папка для сохранения результатов")
    parser.add_argument("-m", "--metric", choices = ["density", "count", "bp"], default = "count", 
                        help="Какое значение вычислять: density (доля покрытия 0-1), count (кол-во TE), bp (кол-во перекрытых пар нуклеотидов)")

    args = parser.parse_args()
    parse_gff3_and_bin(args.input, args.bin_size, args.out_dir, args.metric)