#!/usr/bin/env python3
"""
🧬 Интегрированный анализ Z-DNA и G-квадруплексов
Объединяет результаты Z-Hunt и G-quadruplex анализа
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
        # Переименовываем колонки для соответствия
        if 'chromosome' in g4_data.columns:
            g4_data['Хромосома'] = g4_data['chromosome']
        if 'start' in g4_data.columns:
            g4_data['Позиция'] = g4_data['start']
        if 'sequence' in g4_data.columns:
            g4_data['Последовательность'] = g4_data['sequence']
        print(f"✅ Загружено {len(g4_data)} G-квадруплексов")
        return g4_data
    except Exception as e:
        print(f"⚠️  Ошибка загрузки G4 данных: {e}")
        return pd.DataFrame()

def load_zdna_data():
    """Загружаем данные Z-DNA"""
    zdna_files = [
        'results/zdna_structures_corrected.txt',
        'results/smart_zhunt_results_zdna_structures.txt'
    ]
    
    for file_path in zdna_files:
        try:
            if os.path.exists(file_path):
                # Читаем файл с правильными заголовками
                zdna_data = pd.read_csv(file_path, sep='\t', comment='#', 
                                      names=['Хромосома', 'Позиция', 'Z-score', 'Score1', 'Score2', 'Длина', 'Последовательность'])
                print(f"✅ Загружено {len(zdna_data)} Z-DNA структур из {file_path}")
                return zdna_data
        except Exception as e:
            print(f"⚠️  Ошибка загрузки {file_path}: {e}")
    
    print("⚠️  Z-DNA данные не найдены, создаем пустой DataFrame")
    return pd.DataFrame()

def load_promoter_data():
    """Загружаем данные о промоторах"""
    try:
        # Попробуем найти файлы с промоторным анализом
        promoter_files = [
            'results/promoter_analysis.csv',
            'results/g4_promoter_analysis.txt',
            'results/string_enrichment.csv'
        ]
        
        for file_path in promoter_files:
            if os.path.exists(file_path):
                if file_path.endswith('.csv'):
                    promoter_data = pd.read_csv(file_path)
                else:
                    promoter_data = pd.read_csv(file_path, sep='\t', comment='#')
                print(f"✅ Загружено {len(promoter_data)} промоторных записей из {file_path}")
                return promoter_data
        
        print("⚠️  Промоторные данные не найдены")
        return pd.DataFrame()
    except Exception as e:
        print(f"⚠️  Ошибка загрузки промоторных данных: {e}")
        return pd.DataFrame()

def analyze_colocalization(g4_data, zdna_data, window=1000):
    """Анализируем колокализацию Z-DNA и G-квадруплексов"""
    if zdna_data.empty or g4_data.empty:
        print("⚠️  Недостаточно данных для анализа колокализации")
        return {'total_colocalizations': 0, 'g4_with_zdna': 0, 'zdna_with_g4': 0, 'average_distance': 0, 'colocalized_pairs': []}
    
    print(f"🔍 Анализ колокализации (окно {window} bp)...")
    
    colocalized = []
    
    # Убедимся, что у нас есть правильные колонки
    g4_chrom_col = 'Хромосома' if 'Хромосома' in g4_data.columns else 'chromosome'
    g4_pos_col = 'Позиция' if 'Позиция' in g4_data.columns else 'start'
    g4_seq_col = 'Последовательность' if 'Последовательность' in g4_data.columns else 'sequence'
    
    zdna_chrom_col = 'Хромосома'
    zdna_pos_col = 'Позиция'
    zdna_seq_col = 'Последовательность'
    
    for _, g4 in g4_data.iterrows():
        g4_chrom = g4[g4_chrom_col]
        g4_pos = g4[g4_pos_col]
        
        nearby_zdna = zdna_data[
            (zdna_data[zdna_chrom_col] == g4_chrom) &
            (abs(zdna_data[zdna_pos_col] - g4_pos) <= window)
        ]
        
        if len(nearby_zdna) > 0:
            for _, zdna in nearby_zdna.iterrows():
                distance = abs(zdna[zdna_pos_col] - g4_pos)
                colocalized.append({
                    'chromosome': g4_chrom,
                    'g4_position': g4_pos,
                    'zdna_position': zdna[zdna_pos_col],
                    'distance': distance,
                    'g4_sequence': g4.get(g4_seq_col, ''),
                    'zdna_sequence': zdna.get(zdna_seq_col, ''),
                    'zdna_zscore': zdna.get('Z-score', 0)
                })
    
    colocalization_stats = {
        'total_colocalizations': len(colocalized),
        'g4_with_zdna': len(set([(c['chromosome'], c['g4_position']) for c in colocalized])),
        'zdna_with_g4': len(set([(c['chromosome'], c['zdna_position']) for c in colocalized])),
        'average_distance': np.mean([c['distance'] for c in colocalized]) if colocalized else 0,
        'colocalized_pairs': colocalized
    }
    
    print(f"   📊 Найдено {colocalization_stats['total_colocalizations']} колокализаций")
    print(f"   📊 G4 с Z-DNA: {colocalization_stats['g4_with_zdna']}")
    print(f"   📊 Z-DNA с G4: {colocalization_stats['zdna_with_g4']}")
    
    return colocalization_stats

def create_integrated_visualizations(g4_data, zdna_data, promoter_data, colocalization_stats):
    """Создаем интегрированные визуализации"""
    print("🎨 Создаем интегрированные визуализации...")
    
    # Настройка для русского текста
    plt.rcParams['font.family'] = ['Arial Unicode MS', 'DejaVu Sans', 'sans-serif']
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('🧬 Интегрированный анализ Z-DNA и G-квадруплексов\nDrosophila melanogaster (dm6)', 
                 fontsize=16, fontweight='bold')
    
    # Определяем названия колонок
    g4_chrom_col = 'Хромосома' if 'Хромосома' in g4_data.columns else 'chromosome'
    g4_seq_col = 'Последовательность' if 'Последовательность' in g4_data.columns else 'sequence'
    
    # 1. Распределение по хромосомам
    ax1 = axes[0, 0]
    
    # G4 данные
    if not g4_data.empty:
        g4_chrom_counts = g4_data[g4_chrom_col].value_counts()
        g4_chrom_counts.plot(kind='bar', ax=ax1, alpha=0.7, color='blue', label='G-квадруплексы')
    
    # Z-DNA данные
    if not zdna_data.empty and 'Хромосома' in zdna_data.columns:
        zdna_chrom_counts = zdna_data['Хромосома'].value_counts()
        zdna_chrom_counts.plot(kind='bar', ax=ax1, alpha=0.7, color='red', 
                              label='Z-DNA', width=0.6)
    
    ax1.set_title('Распределение структур по хромосомам')
    ax1.set_xlabel('Хромосома')
    ax1.set_ylabel('Количество структур')
    ax1.legend()
    ax1.tick_params(axis='x', rotation=45)
    
    # 2. Длины последовательностей
    ax2 = axes[0, 1]
    if not g4_data.empty:
        g4_lengths = g4_data[g4_seq_col].str.len()
        ax2.hist(g4_lengths, bins=20, alpha=0.7, color='blue', label='G-квадруплексы', density=True)
    
    if not zdna_data.empty and 'Последовательность' in zdna_data.columns:
        zdna_lengths = zdna_data['Последовательность'].str.len()
        ax2.hist(zdna_lengths, bins=20, alpha=0.7, color='red', label='Z-DNA', density=True)
    
    ax2.set_title('Распределение длин последовательностей')
    ax2.set_xlabel('Длина (bp)')
    ax2.set_ylabel('Плотность')
    ax2.legend()
    
    # 3. Колокализация
    ax3 = axes[0, 2]
    if colocalization_stats.get('colocalized_pairs', []):
        distances = [c['distance'] for c in colocalization_stats['colocalized_pairs']]
        ax3.hist(distances, bins=20, color='green', alpha=0.7)
        ax3.set_title(f'Расстояния при колокализации\n(n={len(distances)})')
        ax3.set_xlabel('Расстояние (bp)')
        ax3.set_ylabel('Количество пар')
    else:
        ax3.text(0.5, 0.5, 'Нет данных о\nколокализации', 
                ha='center', va='center', transform=ax3.transAxes, fontsize=12)
        ax3.set_title('Колокализация')
    
    # 4. Промоторный анализ
    ax4 = axes[1, 0]
    if not promoter_data.empty:
        # Показываем соотношение структур в промоторах vs геноме
        total_g4 = len(g4_data) if not g4_data.empty else 0
        promoter_g4 = len(promoter_data) if not promoter_data.empty else 0
        
        categories = ['Весь геном', 'Промоторы']
        g4_counts = [total_g4, promoter_g4]
        
        bars = ax4.bar(categories, g4_counts, color=['lightblue', 'darkblue'], alpha=0.7)
        ax4.set_title('G-квадруплексы в промоторах')
        ax4.set_ylabel('Количество структур')
        
        # Добавляем проценты на столбцы
        for bar, count in zip(bars, g4_counts):
            height = bar.get_height()
            if total_g4 > 0:
                percentage = (count / total_g4) * 100
                ax4.text(bar.get_x() + bar.get_width()/2., height + height*0.01,
                        f'{count}\n({percentage:.1f}%)', ha='center', va='bottom')
    else:
        ax4.text(0.5, 0.5, 'Нет данных о\nпромоторах', 
                ha='center', va='center', transform=ax4.transAxes, fontsize=12)
        ax4.set_title('Промоторный анализ')
    
    # 5. Z-score распределение (если есть Z-DNA данные)
    ax5 = axes[1, 1]
    if not zdna_data.empty and 'Z-score' in zdna_data.columns:
        zdna_data['Z-score'].hist(bins=30, ax=ax5, color='red', alpha=0.7)
        ax5.set_title('Распределение Z-score')
        ax5.set_xlabel('Z-score')
        ax5.set_ylabel('Количество структур')
        ax5.axvline(zdna_data['Z-score'].mean(), color='black', linestyle='--', 
                   label=f'Среднее: {zdna_data["Z-score"].mean():.1f}')
        ax5.legend()
    else:
        ax5.text(0.5, 0.5, 'Нет данных\nZ-score', 
                ha='center', va='center', transform=ax5.transAxes, fontsize=12)
        ax5.set_title('Z-score распределение')
    
    # 6. Сводная статистика
    ax6 = axes[1, 2]
    ax6.axis('off')
    
    # Создаем текстовую сводку
    summary_text = "📊 СВОДНАЯ СТАТИСТИКА\n\n"
    
    if not g4_data.empty:
        summary_text += f"🔹 G-квадруплексы: {len(g4_data)}\n"
        summary_text += f"   Хромосомы: {g4_data[g4_chrom_col].nunique()}\n"
        summary_text += f"   Средняя длина: {g4_data[g4_seq_col].str.len().mean():.1f} bp\n\n"
    
    if not zdna_data.empty:
        summary_text += f"🔸 Z-DNA структуры: {len(zdna_data)}\n"
        if 'Хромосома' in zdna_data.columns:
            summary_text += f"   Хромосомы: {zdna_data['Хромосома'].nunique()}\n"
        if 'Z-score' in zdna_data.columns:
            summary_text += f"   Средний Z-score: {zdna_data['Z-score'].mean():.1f}\n\n"
    else:
        summary_text += "🔸 Z-DNA: Данные обрабатываются...\n\n"
    
    if not promoter_data.empty:
        summary_text += f"🔹 Промоторы с G4: {len(promoter_data)}\n"
        if not g4_data.empty:
            promo_percent = (len(promoter_data) / len(g4_data)) * 100
            summary_text += f"   Доля от всех G4: {promo_percent:.1f}%\n\n"
    
    if colocalization_stats.get('total_colocalizations', 0) > 0:
        summary_text += f"🔗 Колокализации: {colocalization_stats['total_colocalizations']}\n"
        summary_text += f"   Среднее расстояние: {colocalization_stats['average_distance']:.0f} bp\n"
    else:
        summary_text += "🔗 Колокализации: Анализируется...\n"
    
    ax6.text(0.1, 0.9, summary_text, transform=ax6.transAxes, fontsize=11,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('results/integrated_analysis.png', dpi=300, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.show()
    print("✅ Интегрированная визуализация сохранена: results/integrated_analysis.png")

def generate_final_report(g4_data, zdna_data, promoter_data, colocalization_stats):
    """Генерируем финальный отчет"""
    print("📝 Генерируем финальный отчет...")
    
    report = """# 🧬 ФИНАЛЬНЫЙ ОТЧЕТ: Анализ Z-DNA и G-квадруплексов
## Drosophila melanogaster (dm6)

### 📊 ОСНОВНЫЕ РЕЗУЛЬТАТЫ

"""
    
    # G-квадруплексы
    if not g4_data.empty:
        report += f"""#### 🔹 G-квадруплексы
- **Всего найдено**: {len(g4_data)} структур
- **Хромосомы**: {g4_data['Хромосома'].nunique()} различных
- **Средняя длина**: {g4_data['Последовательность'].str.len().mean():.1f} bp
- **Диапазон длин**: {g4_data['Последовательность'].str.len().min()}-{g4_data['Последовательность'].str.len().max()} bp

**Распределение по хромосомам:**
"""
        for chrom, count in g4_data['Хромосома'].value_counts().items():
            percentage = (count / len(g4_data)) * 100
            report += f"- {chrom}: {count} ({percentage:.1f}%)\n"
        report += "\n"
    
    # Z-DNA
    if not zdna_data.empty:
        report += f"""#### 🔸 Z-DNA структуры
- **Всего найдено**: {len(zdna_data)} структур
"""
        if 'Хромосома' in zdna_data.columns:
            report += f"- **Хромосомы**: {zdna_data['Хромосома'].nunique()} различных\n"
        
        if 'Z-score' in zdna_data.columns:
            report += f"""- **Z-score диапазон**: {zdna_data['Z-score'].min():.1f}-{zdna_data['Z-score'].max():.1f}
- **Средний Z-score**: {zdna_data['Z-score'].mean():.1f}
"""
        
        if 'Хромосома' in zdna_data.columns:
            report += "\n**Распределение по хромосомам:**\n"
            for chrom, count in zdna_data['Хромосома'].value_counts().items():
                percentage = (count / len(zdna_data)) * 100
                report += f"- {chrom}: {count} ({percentage:.1f}%)\n"
        report += "\n"
    else:
        report += """#### 🔸 Z-DNA структуры
- **Статус**: Данные обрабатываются или не найдены с текущими параметрами (Z-score 300-400)
- **Рекомендация**: Возможно, следует расширить диапазон Z-score для поиска

"""
    
    # Промоторы
    if not promoter_data.empty:
        report += f"""#### 🔹 Промоторный анализ
- **G4 в промоторах**: {len(promoter_data)} структур
- **Уникальные гены**: {promoter_data.get('Gene_ID', promoter_data.get('gene_id', pd.Series())).nunique() if 'Gene_ID' in promoter_data.columns or 'gene_id' in promoter_data.columns else 'N/A'}
"""
        if not g4_data.empty:
            promo_percent = (len(promoter_data) / len(g4_data)) * 100
            report += f"- **Доля от всех G4**: {promo_percent:.1f}%\n"
        report += "\n"
    
    # Колокализация
    if colocalization_stats.get('total_colocalizations', 0) > 0:
        report += f"""#### 🔗 Колокализация Z-DNA и G-квадруплексов
- **Всего колокализаций**: {colocalization_stats['total_colocalizations']}
- **G4 с близкими Z-DNA**: {colocalization_stats['g4_with_zdna']}
- **Z-DNA с близкими G4**: {colocalization_stats['zdna_with_g4']}
- **Среднее расстояние**: {colocalization_stats['average_distance']:.0f} bp
- **Окно поиска**: ±1000 bp

"""
    else:
        report += """#### 🔗 Колокализация Z-DNA и G-квадруплексов
- **Статус**: Недостаточно данных для анализа колокализации
- **Возможные причины**: 
  - Z-DNA структуры не найдены с текущими параметрами
  - Структуры находятся на разных участках генома
  - Необходимо расширить окно поиска

"""
    
    # Биологическая интерпретация
    report += """### 🧬 БИОЛОГИЧЕСКАЯ ИНТЕРПРЕТАЦИЯ

#### Функциональное значение
"""
    
    if not g4_data.empty:
        report += """
**G-квадруплексы:**
- Участвуют в регуляции транскрипции
- Влияют на репликацию ДНК
- Играют роль в теломерной биологии
- Могут вызывать геномную нестабильность
"""
    
    if not zdna_data.empty:
        report += """
**Z-DNA структуры:**
- Связаны с активной транскрипцией
- Могут индуцировать рекомбинацию
- Влияют на хроматиновую структуру
- Участвуют в эпигенетической регуляции
"""
    
    if colocalization_stats.get('total_colocalizations', 0) > 0:
        report += f"""
**Колокализация:**
- Обнаружено {colocalization_stats['total_colocalizations']} случаев совместного присутствия
- Может указывать на функциональную связь между структурами
- Требует дальнейшего исследования механизмов взаимодействия
"""
    
    # Технические детали
    report += """

### 🔬 МЕТОДОЛОГИЯ

**Параметры анализа:**
- **Геном**: Drosophila melanogaster dm6
- **G-квадруплексы**: Паттерн G{3,}N{1,7}G{3,}N{1,7}G{3,}N{1,7}G{3,}
- **Z-DNA**: Z-Hunt с параметрами (12 8 12), Z-score 300-400
- **Промоторы**: ±1000 bp от TSS
- **Колокализация**: Окно ±1000 bp

**Инструменты:**
- Z-Hunt для предсказания Z-DNA
- Регулярные выражения для G-квадруплексов
- STRING DB для функционального обогащения
- Python + pandas/matplotlib для анализа

"""
    
    # Сохраняем отчет
    with open('results/FINAL_INTEGRATED_REPORT.md', 'w', encoding='utf-8') as f:
        f.write(report)
    
    print("✅ Финальный отчет сохранен: results/FINAL_INTEGRATED_REPORT.md")

def main():
    print("🧬 ИНТЕГРИРОВАННЫЙ АНАЛИЗ Z-DNA и G-КВАДРУПЛЕКСОВ")
    print("=" * 60)
    
    # Создаем папку результатов
    os.makedirs('results', exist_ok=True)
    
    # Загружаем данные
    print("\n📂 Загрузка данных...")
    g4_data = load_g4_data()
    zdna_data = load_zdna_data()
    promoter_data = load_promoter_data()
    
    # Анализ колокализации
    print("\n🔍 Анализ колокализации...")
    colocalization_stats = analyze_colocalization(g4_data, zdna_data)
    
    # Создаем визуализации
    print("\n🎨 Создание визуализаций...")
    create_integrated_visualizations(g4_data, zdna_data, promoter_data, colocalization_stats)
    
    # Генерируем отчет
    print("\n📝 Генерация отчета...")
    generate_final_report(g4_data, zdna_data, promoter_data, colocalization_stats)
    
    print("\n" + "=" * 60)
    print("🎉 ИНТЕГРИРОВАННЫЙ АНАЛИЗ ЗАВЕРШЕН!")
    print("📁 Результаты сохранены в папке results/")
    print("   📊 Визуализация: integrated_analysis.png")
    print("   📝 Отчет: FINAL_INTEGRATED_REPORT.md")

if __name__ == "__main__":
    main() 