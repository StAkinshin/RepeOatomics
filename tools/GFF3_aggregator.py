import argparse
import os
import glob
import tempfile
import shutil

def aggregate_gff3(input_files, output_file):
    # Множество для хранения уникальных мета-заголовков, чтобы избежать их дублирования
    seen_headers = set()
    
    # Флаг наличия FASTA-блока хотя бы в одном из файлов
    has_fasta = False

    print(f"Начало объединения {len(input_files)} файлов...")

    with open(output_file, 'w', encoding='utf-8') as out_f:
        # 1. Записываем обязательный заголовок версии только один раз
        out_f.write("##gff-version 3\n")
        
        # 2. Создаем временный файл для FASTA-последовательностей.
        # Это нужно, чтобы не загружать огромные геномы в оперативную память.
        with tempfile.NamedTemporaryFile('w+', encoding='utf-8', delete=False) as temp_fasta:
            temp_fasta_name = temp_fasta.name
            
            for infile in input_files:
                if not os.path.isfile(infile):
                    print(f"Предупреждение: Файл {infile} не найден. Пропуск.")
                    continue
                
                print(f"Обработка: {infile}")
                in_fasta_mode = False
                last_char = '\n'
                last_char_fasta = '\n'

                with open(infile, 'r', encoding='utf-8') as in_f:
                    for line in in_f:
                        # Если мы уже дошли до блока FASTA в текущем файле
                        if in_fasta_mode:
                            temp_fasta.write(line)
                            last_char_fasta = line[-1]
                            continue

                        # Игнорируем версию (мы ее уже записали)
                        if line.startswith("##gff-version"):
                            continue
                        
                        # Если встретили начало FASTA-блока
                        if line.startswith("##FASTA"):
                            in_fasta_mode = True
                            has_fasta = True
                            continue
                            
                        # Обработка других прагм и метаданных (начинаются с ##)
                        if line.startswith("##"):
                            # Директива ### используется для закрытия блоков, её можно дублировать
                            if line.strip() != "###":
                                if line not in seen_headers:
                                    seen_headers.add(line)
                                    out_f.write(line)
                                continue
                        
                        # Записываем сами данные признаков (колонки GFF3) или обычные комментарии (#)
                        out_f.write(line)
                        last_char = line[-1]
                
                # Добавляем перенос строки, если файл неожиданно закончился без него
                if not in_fasta_mode and last_char != '\n':
                    out_f.write('\n')
                elif in_fasta_mode and last_char_fasta != '\n':
                    temp_fasta.write('\n')

        # 3. Если были найдены секции FASTA, добавляем их в самый конец объединенного файла
        if has_fasta:
            print("Добавление объединенных FASTA-последовательностей в конец файла...")
            out_f.write("##FASTA\n")
            # Перемещаем указатель в начало временного файла перед чтением
            temp_fasta.seek(0)
            with open(temp_fasta_name, 'r', encoding='utf-8') as temp_fasta_read:
                # Быстрое копирование содержимого
                shutil.copyfileobj(temp_fasta_read, out_f)
        
        # Удаляем временный файл
        os.remove(temp_fasta_name)

    print(f"\nУспешно! Все данные сохранены в файл: {output_file}")

def get_files_from_args(inputs):
    """Раскрывает маски (например *.gff3) и собирает список файлов"""
    files = []
    for item in inputs:
        # glob.glob раскрывает шаблоны вроде *.gff3
        matched = glob.glob(item)
        if matched:
            files.extend(matched)
        else:
            files.append(item)
    return sorted(list(set(files)))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Аггрегация нескольких GFF3 файлов в один с сохранением структуры.")
    parser.add_argument("-i", "--input", nargs='+', required=True, 
                        help="Пути к входным GFF3 файлам. Можно использовать маски (например: data/*.gff3)")
    parser.add_argument("-o", "--output", required=True, 
                        help="Путь к итоговому (объединенному) GFF3 файлу")

    args = parser.parse_args()
    
    input_files_list = get_files_from_args(args.input)
    
    if not input_files_list:
        print("Ошибка: Входные файлы не найдены!")
    else:
        aggregate_gff3(input_files_list, args.output)