#!/usr/bin/env python3
"""
🧬 Исправленное извлечение Z-DNA структур из .probability файлов
Формат файла: позиция, score1, score2, z_score, последовательность
"""

import os
import re
import json

def extract_zdna_from_probability(file_path, min_zscore=300, max_zscore=400):
    """Извлекаем Z-DNA структуры из .probability файла"""
    zdna_structures = []
    
    if not os.path.exists(file_path):
        print(f"⚠️  Файл не найден: {file_path}")
        return zdna_structures
    
    print(f"📊 Обрабатываем {file_path}...")
    
    with open(file_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            try:
                # Парсим строку: позиция, score1, score2, z_score, последовательность
                parts = line.split()
                if len(parts) < 5:
                    continue
                
                position = int(parts[0])
                score1 = float(parts[1])
                score2 = float(parts[2])
                z_score_str = parts[3]
                sequence = parts[4]
                
                # Обрабатываем z_score (может быть в научной нотации)
                z_score = float(z_score_str)
                
                # Проверяем, попадает ли в диапазон Z-score
                if min_zscore <= z_score <= max_zscore:
                    zdna_structures.append({
                        'position': position,
                        'z_score': z_score,
                        'score1': score1,
                        'score2': score2,
                        'sequence': sequence,
                        'length': len(sequence)
                    })
                    
            except (ValueError, IndexError) as e:
                if line_num <= 10:  # Показываем только первые несколько ошибок
                    print(f"⚠️  Ошибка в строке {line_num}: {line[:50]}... - {e}")
                continue
    
    print(f"✅ Найдено {len(zdna_structures)} Z-DNA структур")
    return zdna_structures

def main():
    # Настройки
    min_zscore = 300
    max_zscore = 400
    z_hunt_dir = "z_hunt_results"
    output_dir = "results"
    
    print(f"🧬 ИЗВЛЕЧЕНИЕ Z-DNA СТРУКТУР (Z-score {min_zscore}-{max_zscore})")
    print("=" * 60)
    
    # Создаем папку результатов
    os.makedirs(output_dir, exist_ok=True)
    
    # Список хромосом
    chromosomes = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX', 'chrY']
    
    all_zdna = []
    chromosome_stats = {}
    
    for chrom in chromosomes:
        prob_file = os.path.join(z_hunt_dir, f"{chrom}.fa.probability")
        
        if os.path.exists(prob_file):
            zdna_structures = extract_zdna_from_probability(prob_file, min_zscore, max_zscore)
            all_zdna.extend(zdna_structures)
            
            # Добавляем информацию о хромосоме к каждой структуре
            for struct in zdna_structures:
                struct['chromosome'] = chrom
            
            chromosome_stats[chrom] = {
                'count': len(zdna_structures),
                'avg_zscore': sum(s['z_score'] for s in zdna_structures) / len(zdna_structures) if zdna_structures else 0,
                'max_zscore': max(s['z_score'] for s in zdna_structures) if zdna_structures else 0,
                'min_zscore': min(s['z_score'] for s in zdna_structures) if zdna_structures else 0
            }
        else:
            print(f"⚠️  Файл не найден: {prob_file}")
            chromosome_stats[chrom] = {'count': 0, 'avg_zscore': 0, 'max_zscore': 0, 'min_zscore': 0}
    
    print("\n" + "=" * 60)
    print(f"🎉 ОБЩИЕ РЕЗУЛЬТАТЫ:")
    print(f"📊 Всего найдено Z-DNA структур: {len(all_zdna)}")
    
    # Статистика по хромосомам
    print(f"\n📈 Статистика по хромосомам:")
    for chrom, stats in chromosome_stats.items():
        if stats['count'] > 0:
            print(f"   {chrom}: {stats['count']} структур, Z-score: {stats['min_zscore']:.1f}-{stats['max_zscore']:.1f} (avg: {stats['avg_zscore']:.1f})")
        else:
            print(f"   {chrom}: 0 структур")
    
    # Сохраняем результаты
    output_file = os.path.join(output_dir, "zdna_structures_corrected.txt")
    with open(output_file, 'w') as f:
        f.write("# Z-DNA структуры (Z-score 300-400)\n")
        f.write("# Хромосома\tПозиция\tZ-score\tScore1\tScore2\tДлина\tПоследовательность\n")
        
        for struct in sorted(all_zdna, key=lambda x: (x['chromosome'], x['position'])):
            f.write(f"{struct['chromosome']}\t{struct['position']}\t{struct['z_score']:.3f}\t")
            f.write(f"{struct['score1']:.3f}\t{struct['score2']:.3f}\t{struct['length']}\t{struct['sequence']}\n")
    
    # Сохраняем JSON сводку
    summary = {
        'total_structures': len(all_zdna),
        'parameters': {
            'min_zscore': min_zscore,
            'max_zscore': max_zscore
        },
        'chromosome_stats': chromosome_stats,
        'processing_info': {
            'processed_files': len([c for c in chromosomes if chromosome_stats[c]['count'] >= 0]),
            'successful_files': len([c for c in chromosomes if chromosome_stats[c]['count'] > 0])
        }
    }
    
    summary_file = os.path.join(output_dir, "zdna_summary_corrected.json")
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    
    print(f"\n📄 Результаты сохранены:")
    print(f"   📊 Z-DNA структуры: {output_file}")
    print(f"   📋 Сводка: {summary_file}")
    print("\n✅ Извлечение завершено!")

if __name__ == "__main__":
    main() 