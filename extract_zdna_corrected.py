#!/usr/bin/env python3
"""
Извлечение Z-DNA структур из результатов Z-Hunt с исправлениями
Обработка файлов .probability и фильтрация по Z-score
"""

import os
import json
import pandas as pd
from collections import defaultdict

def extract_zdna_structures(zhunt_dir="z_hunt_results", min_zscore=300, max_zscore=400):
    """
    Извлекает Z-DNA структуры из файлов Z-Hunt .probability
    """
    print("🧬 ИЗВЛЕЧЕНИЕ Z-DNA СТРУКТУР (Z-score {}-{})".format(min_zscore, max_zscore))
    print("=" * 60)
    
    all_structures = []
    chromosome_stats = {}
    
    # Обрабатываем каждый файл .probability
    for filename in os.listdir(zhunt_dir):
        if filename.endswith(".probability"):
            filepath = os.path.join(zhunt_dir, filename)
            chrom = filename.replace(".fa.probability", "")
            
            print(f"📊 Обрабатываем {filepath}...")
            
            structures = []
            zscores = []
            
            try:
                with open(filepath, 'r') as f:
                    for line_num, line in enumerate(f, 1):
                        line = line.strip()
                        if not line or line.startswith('#'):
                            continue
                        
                        parts = line.split()
                        if len(parts) < 4:
                            continue
                        
                        try:
                            position = int(parts[0])
                            zscore = float(parts[3])  # Z-score в 4-й колонке
                            
                            # Фильтруем по Z-score
                            if min_zscore <= zscore <= max_zscore:
                                structures.append({
                                    'chromosome': chrom,
                                    'start': position,
                                    'end': position + 11,  # Окно Z-Hunt ~12bp
                                    'zscore': zscore
                                })
                                zscores.append(zscore)
                                
                        except (ValueError, IndexError) as e:
                            # Пропускаем проблемные строки
                            continue
                            
            except Exception as e:
                print(f"❌ Ошибка при обработке {filepath}: {e}")
                continue
            
            # Сохраняем статистику по хромосоме
            if zscores:
                chromosome_stats[chrom] = {
                    'count': len(structures),
                    'min_zscore': min(zscores),
                    'max_zscore': max(zscores),
                    'avg_zscore': sum(zscores) / len(zscores)
                }
                
                print(f"✅ Найдено {len(structures)} Z-DNA структур")
            else:
                print(f"⚠️ Структуры не найдены")
                chromosome_stats[chrom] = {
                    'count': 0,
                    'min_zscore': 0,
                    'max_zscore': 0,
                    'avg_zscore': 0
                }
            
            all_structures.extend(structures)
    
    return all_structures, chromosome_stats

def save_results(structures, stats, output_dir="results"):
    """
    Сохраняет результаты в файлы
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Сохраняем структуры в текстовый файл
    structures_file = os.path.join(output_dir, "zdna_structures_corrected.txt")
    with open(structures_file, 'w') as f:
        f.write("chromosome\tstart\tend\tzscore\n")
        for struct in structures:
            f.write(f"{struct['chromosome']}\t{struct['start']}\t{struct['end']}\t{struct['zscore']:.1f}\n")
    
    # Сохраняем сводку в JSON
    summary_file = os.path.join(output_dir, "zdna_summary_corrected.json")
    summary = {
        'total_structures': len(structures),
        'chromosome_stats': stats,
        'filtering': {
            'min_zscore': 300,
            'max_zscore': 400
        }
    }
    
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    return structures_file, summary_file

def main():
    # Извлекаем структуры
    structures, stats = extract_zdna_structures()
    
    print("=" * 60)
    print("🎉 ОБЩИЕ РЕЗУЛЬТАТЫ:")
    print(f"📊 Всего найдено Z-DNA структур: {len(structures)}")
    
    print("📈 Статистика по хромосомам:")
    for chrom, stat in sorted(stats.items()):
        print(f"   {chrom}: {stat['count']} структур, Z-score: {stat['min_zscore']:.1f}-{stat['max_zscore']:.1f} (avg: {stat['avg_zscore']:.1f})")
    
    # Сохраняем результаты
    structures_file, summary_file = save_results(structures, stats)
    
    print("📄 Результаты сохранены:")
    print(f"   📊 Z-DNA структуры: {structures_file}")
    print(f"   📋 Сводка: {summary_file}")
    print("✅ Извлечение завершено!")

if __name__ == "__main__":
    main()