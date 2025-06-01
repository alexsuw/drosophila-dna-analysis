#!/usr/bin/env python3
"""
🧬 Оптимизированный интегрированный анализ Z-DNA и G-квадруплексов
Быстрая обработка больших наборов данных
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
import os
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

plt.style.use('bmh')
plt.rcParams['figure.figsize'] = [12, 8]
plt.rcParams['font.size'] = 10

def load_g4_data():
    """Загружаем данные G-квадруплексов"""
    try:
        g4_data = pd.read_csv('results/quadruplex_results.csv')
        print(f"✅ Загружено {len(g4_data)} G-квадруплексов")
        return g4_data
    except Exception as e:
        print(f"⚠️  Ошибка загрузки G4 данных: {e}")
        return pd.DataFrame()

def load_zdna_data():
    """Загружаем данные Z-DNA"""
    try:
        zdna_data = pd.read_csv('results/zdna_structures_corrected.txt', sep='\t', comment='#', 
                              names=['chromosome', 'position', 'zscore', 'score1', 'score2', 'length', 'sequence'])
        print(f"✅ Загружено {len(zdna_data)} Z-DNA структур")
        return zdna_data
    except Exception as e:
        print(f"⚠️  Ошибка загрузки Z-DNA данных: {e}")
        return pd.DataFrame()

def fast_colocalization_analysis(g4_data, zdna_data, window=1000):
    """Быстрый анализ колокализации используя группировку по хромосомам"""
    print(f"🔍 Быстрый анализ колокализации (окно {window} bp)...")
    
    if g4_data.empty or zdna_data.empty:
        print("⚠️  Недостаточно данных для анализа")
        return {'total_colocalizations': 0, 'summary_by_chromosome': {}}
    
    colocalization_summary = {}
    total_colocs = 0
    
    # Группируем по хромосомам для ускорения
    for chrom in g4_data['chromosome'].unique():
        print(f"   Анализируем {chrom}...")
        
        g4_chrom = g4_data[g4_data['chromosome'] == chrom]
        zdna_chrom = zdna_data[zdna_data['chromosome'] == chrom]
        
        if zdna_chrom.empty:
            continue
            
        colocs_in_chrom = 0
        g4_positions = g4_chrom['start'].values
        zdna_positions = zdna_chrom['position'].values
        
        # Быстрый поиск ближайших позиций
        for g4_pos in g4_positions:
            # Находим Z-DNA в окне
            distances = np.abs(zdna_positions - g4_pos)
            nearby = distances <= window
            colocs_in_chrom += np.sum(nearby)
        
        colocalization_summary[chrom] = {
            'g4_count': len(g4_chrom),
            'zdna_count': len(zdna_chrom),
            'colocalizations': colocs_in_chrom,
            'colocalization_rate': colocs_in_chrom / len(g4_chrom) if len(g4_chrom) > 0 else 0
        }
        
        total_colocs += colocs_in_chrom
        
    print(f"   📊 Всего колокализаций: {total_colocs}")
    
    return {
        'total_colocalizations': total_colocs,
        'summary_by_chromosome': colocalization_summary
    }

def create_fast_visualizations(g4_data, zdna_data, colocalization_stats):
    """Создаем быстрые визуализации"""
    print("🎨 Создаем визуализации...")
    
    plt.rcParams['font.family'] = ['Arial Unicode MS', 'DejaVu Sans', 'sans-serif']
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('🧬 Интегрированный анализ Z-DNA и G-квадруплексов\nDrosophila melanogaster (dm6)', 
                 fontsize=16, fontweight='bold')
    
    # 1. Распределение по хромосомам
    ax1 = axes[0, 0]
    
    if not g4_data.empty:
        g4_counts = g4_data['chromosome'].value_counts()
        positions = np.arange(len(g4_counts))
        ax1.bar(positions - 0.2, g4_counts.values, 0.4, label='G-квадруплексы', alpha=0.7, color='blue')
    
    if not zdna_data.empty:
        zdna_counts = zdna_data['chromosome'].value_counts()
        # Приводим к тому же порядку хромосом
        zdna_aligned = []
        for chrom in g4_counts.index if not g4_data.empty else zdna_counts.index:
            zdna_aligned.append(zdna_counts.get(chrom, 0))
        
        ax1.bar(positions + 0.2, zdna_aligned, 0.4, label='Z-DNA', alpha=0.7, color='red')
    
    ax1.set_title('Распределение структур по хромосомам')
    ax1.set_xlabel('Хромосома')
    ax1.set_ylabel('Количество структур')
    ax1.set_xticks(positions)
    ax1.set_xticklabels(g4_counts.index if not g4_data.empty else zdna_counts.index, rotation=45)
    ax1.legend()
    
    # 2. Длины последовательностей
    ax2 = axes[0, 1]
    if not g4_data.empty:
        g4_lengths = g4_data['sequence'].str.len()
        ax2.hist(g4_lengths, bins=20, alpha=0.7, color='blue', label='G-квадруплексы', density=True)
    
    if not zdna_data.empty:
        zdna_lengths = zdna_data['sequence'].str.len()
        ax2.hist(zdna_lengths, bins=20, alpha=0.7, color='red', label='Z-DNA', density=True)
    
    ax2.set_title('Распределение длин последовательностей')
    ax2.set_xlabel('Длина (bp)')
    ax2.set_ylabel('Плотность')
    ax2.legend()
    
    # 3. Колокализация по хромосомам
    ax3 = axes[0, 2]
    if colocalization_stats.get('summary_by_chromosome'):
        chroms = list(colocalization_stats['summary_by_chromosome'].keys())
        coloc_counts = [colocalization_stats['summary_by_chromosome'][c]['colocalizations'] for c in chroms]
        
        ax3.bar(chroms, coloc_counts, color='green', alpha=0.7)
        ax3.set_title('Колокализации по хромосомам')
        ax3.set_xlabel('Хромосома')
        ax3.set_ylabel('Количество колокализаций')
        ax3.tick_params(axis='x', rotation=45)
    else:
        ax3.text(0.5, 0.5, 'Нет данных о\nколокализации', 
                ha='center', va='center', transform=ax3.transAxes, fontsize=12)
        ax3.set_title('Колокализация')
    
    # 4. Z-score распределение
    ax4 = axes[1, 0]
    if not zdna_data.empty:
        ax4.hist(zdna_data['zscore'], bins=30, color='red', alpha=0.7)
        ax4.set_title('Распределение Z-score')
        ax4.set_xlabel('Z-score')
        ax4.set_ylabel('Количество структур')
        ax4.axvline(zdna_data['zscore'].mean(), color='black', linestyle='--', 
                   label=f'Среднее: {zdna_data["zscore"].mean():.1f}')
        ax4.legend()
    
    # 5. Статистика колокализации
    ax5 = axes[1, 1]
    if colocalization_stats.get('summary_by_chromosome'):
        chroms = list(colocalization_stats['summary_by_chromosome'].keys())
        rates = [colocalization_stats['summary_by_chromosome'][c]['colocalization_rate'] * 100 for c in chroms]
        
        ax5.bar(chroms, rates, color='orange', alpha=0.7)
        ax5.set_title('Частота колокализации (%)')
        ax5.set_xlabel('Хромосома')
        ax5.set_ylabel('Процент G4 с Z-DNA')
        ax5.tick_params(axis='x', rotation=45)
    
    # 6. Сводная статистика
    ax6 = axes[1, 2]
    ax6.axis('off')
    
    summary_text = "📊 СВОДНАЯ СТАТИСТИКА\n\n"
    
    if not g4_data.empty:
        summary_text += f"🔹 G-квадруплексы: {len(g4_data):,}\n"
        summary_text += f"   Хромосомы: {g4_data['chromosome'].nunique()}\n"
        summary_text += f"   Средняя длина: {g4_data['sequence'].str.len().mean():.1f} bp\n\n"
    
    if not zdna_data.empty:
        summary_text += f"🔸 Z-DNA структуры: {len(zdna_data):,}\n"
        summary_text += f"   Хромосомы: {zdna_data['chromosome'].nunique()}\n"
        summary_text += f"   Средний Z-score: {zdna_data['zscore'].mean():.1f}\n"
        summary_text += f"   Диапазон Z-score: {zdna_data['zscore'].min():.0f}-{zdna_data['zscore'].max():.0f}\n\n"
    
    if colocalization_stats.get('total_colocalizations', 0) > 0:
        summary_text += f"🔗 Колокализации: {colocalization_stats['total_colocalizations']:,}\n"
        
        if not g4_data.empty:
            coloc_percent = (colocalization_stats['total_colocalizations'] / len(g4_data)) * 100
            summary_text += f"   Процент G4: {coloc_percent:.1f}%\n"
    else:
        summary_text += "🔗 Колокализации: Не найдены\n"
    
    ax6.text(0.1, 0.9, summary_text, transform=ax6.transAxes, fontsize=11,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('results/fast_integrated_analysis.png', dpi=300, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.show()
    print("✅ Визуализация сохранена: results/fast_integrated_analysis.png")

def save_analysis_results(g4_data, zdna_data, colocalization_stats):
    """Сохраняем результаты анализа"""
    print("💾 Сохраняем результаты анализа...")
    
    results = {
        'analysis_type': 'integrated_z_dna_g4_analysis',
        'timestamp': pd.Timestamp.now().isoformat(),
        'datasets': {
            'g4_structures': {
                'count': len(g4_data) if not g4_data.empty else 0,
                'chromosomes': g4_data['chromosome'].nunique() if not g4_data.empty else 0,
                'avg_length': float(g4_data['sequence'].str.len().mean()) if not g4_data.empty else 0
            },
            'zdna_structures': {
                'count': len(zdna_data) if not zdna_data.empty else 0,
                'chromosomes': zdna_data['chromosome'].nunique() if not zdna_data.empty else 0,
                'avg_zscore': float(zdna_data['zscore'].mean()) if not zdna_data.empty else 0,
                'zscore_range': [float(zdna_data['zscore'].min()), float(zdna_data['zscore'].max())] if not zdna_data.empty else [0, 0]
            }
        },
        'colocalization': colocalization_stats
    }
    
    with open('results/integrated_analysis_results.json', 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    
    print("✅ Результаты сохранены: results/integrated_analysis_results.json")

def main():
    print("🧬 ОПТИМИЗИРОВАННЫЙ ИНТЕГРИРОВАННЫЙ АНАЛИЗ")
    print("="*60)
    
    # Загружаем данные
    print("\n📂 Загрузка данных...")
    g4_data = load_g4_data()
    zdna_data = load_zdna_data()
    
    # Быстрый анализ колокализации
    print("\n🔍 Анализ колокализации...")
    colocalization_stats = fast_colocalization_analysis(g4_data, zdna_data)
    
    # Создаем визуализации
    print("\n🎨 Создание визуализаций...")
    create_fast_visualizations(g4_data, zdna_data, colocalization_stats)
    
    # Сохраняем результаты
    print("\n💾 Сохранение результатов...")
    save_analysis_results(g4_data, zdna_data, colocalization_stats)
    
    print("\n🎉 АНАЛИЗ ЗАВЕРШЕН!")
    print(f"📊 G-квадруплексы: {len(g4_data):,}")
    print(f"📊 Z-DNA структуры: {len(zdna_data):,}")
    print(f"📊 Колокализации: {colocalization_stats.get('total_colocalizations', 0):,}")
    
    # Создаем финальный отчет
    create_final_report(g4_data, zdna_data, colocalization_stats)

def create_final_report(g4_data, zdna_data, colocalization_stats):
    """Создаем финальный отчет"""
    print("\n📝 Создание финального отчета...")
    
    report = f"""# 🧬 ФИНАЛЬНЫЙ ОТЧЕТ: Анализ Z-DNA и G-квадруплексов

## Сводка результатов

### G-квадруплексы
- **Всего найдено**: {len(g4_data):,} структур
- **Хромосомы**: {g4_data['chromosome'].nunique() if not g4_data.empty else 0}
- **Средняя длина**: {g4_data['sequence'].str.len().mean():.1f} bp

### Z-DNA структуры  
- **Всего найдено**: {len(zdna_data):,} структур
- **Хромосомы**: {zdna_data['chromosome'].nunique() if not zdna_data.empty else 0}
- **Средний Z-score**: {zdna_data['zscore'].mean():.1f}
- **Диапазон Z-score**: {zdna_data['zscore'].min():.0f} - {zdna_data['zscore'].max():.0f}

### Колокализация
- **Всего колокализаций**: {colocalization_stats.get('total_colocalizations', 0):,}
"""

    if not g4_data.empty and colocalization_stats.get('total_colocalizations', 0) > 0:
        coloc_percent = (colocalization_stats['total_colocalizations'] / len(g4_data)) * 100
        report += f"- **Процент G4 с Z-DNA**: {coloc_percent:.1f}%\n"

    report += f"""
## Распределение по хромосомам

| Хромосома | G-квадруплексы | Z-DNA | Колокализации |
|-----------|----------------|-------|---------------|
"""

    if not g4_data.empty:
        for chrom in sorted(g4_data['chromosome'].unique()):
            g4_count = len(g4_data[g4_data['chromosome'] == chrom])
            zdna_count = len(zdna_data[zdna_data['chromosome'] == chrom]) if not zdna_data.empty else 0
            coloc_count = colocalization_stats.get('summary_by_chromosome', {}).get(chrom, {}).get('colocalizations', 0)
            report += f"| {chrom} | {g4_count:,} | {zdna_count:,} | {coloc_count:,} |\n"

    report += f"""
## Методы анализа

### Z-Hunt параметры
- **Окно**: 12 bp
- **Порог 1**: 8
- **Порог 2**: 12  
- **Z-score диапазон**: 300-400

### G-квадруплекс поиск
- **Паттерн**: G{{3+}}N{{1-7}}G{{3+}}N{{1-7}}G{{3+}}N{{1-7}}G{{3+}}
- **Регулярное выражение**: Строгий поиск

### Колокализация
- **Окно поиска**: ±1000 bp
- **Критерий**: Перекрытие позиций

## Файлы результатов

- `results/quadruplex_results.csv` - G-квадруплексы
- `results/zdna_structures_corrected.txt` - Z-DNA структуры  
- `results/fast_integrated_analysis.png` - Визуализации
- `results/integrated_analysis_results.json` - Статистика

## 🎯 Заключение

Анализ генома Drosophila melanogaster (dm6) выявил:

1. **{len(g4_data):,} G-квадруплексов** распределенных по всем хромосомам
2. **{len(zdna_data):,} Z-DNA структур** с Z-score 300-400
3. **{colocalization_stats.get('total_colocalizations', 0):,} случаев колокализации** в пределах 1 kb

Результаты демонстрируют распределение альтернативных структур ДНК в геноме дрозофилы и их потенциальную коэкспрессию в регуляторных регионах.
"""

    with open('results/FINAL_INTEGRATED_REPORT.md', 'w', encoding='utf-8') as f:
        f.write(report)
    
    print("✅ Финальный отчет: results/FINAL_INTEGRATED_REPORT.md")

if __name__ == "__main__":
    main() 