#!/usr/bin/env python3
"""
🚀 SMART PARALLEL Z-HUNT WITH PROGRESS MONITORING
- Асинхронный анализ всех хромосом
- Мониторинг прогресса в реальном времени  
- БЕЗ удаления промежуточных файлов
- Показ статуса каждого процесса
"""

import subprocess
import time
import os
import sys
import threading
from pathlib import Path
import concurrent.futures
import json
from datetime import datetime
import psutil

class ZHuntProgressMonitor:
    def __init__(self, work_dir):
        self.work_dir = Path(work_dir)
        self.processes = {}
        self.progress_data = {}
        self.start_time = time.time()
        self.monitoring = True
        
    def start_monitoring(self):
        """Запуск мониторинга в отдельном потоке"""
        monitor_thread = threading.Thread(target=self._monitor_loop, daemon=True)
        monitor_thread.start()
        
    def _monitor_loop(self):
        """Основной цикл мониторинга"""
        while self.monitoring:
            self._update_progress()
            self._display_status()
            time.sleep(2)  # Обновление каждые 2 секунды
            
    def _update_progress(self):
        """Обновление данных прогресса"""
        for chr_name in self.processes:
            if chr_name in self.processes:
                # Проверяем размер Z-SCORE файла
                zscore_file = self.work_dir / f"{chr_name}.fa.Z-SCORE"
                prob_file = self.work_dir / f"{chr_name}.fa.probability"
                
                if zscore_file.exists():
                    size_mb = zscore_file.stat().st_size / (1024 * 1024)
                    self.progress_data[chr_name] = {
                        'status': 'calculating_zscore',
                        'size_mb': size_mb,
                        'file_exists': True
                    }
                elif prob_file.exists():
                    size_mb = prob_file.stat().st_size / (1024 * 1024)
                    self.progress_data[chr_name] = {
                        'status': 'completed',
                        'size_mb': size_mb,
                        'file_exists': True
                    }
                else:
                    self.progress_data[chr_name] = {
                        'status': 'starting',
                        'size_mb': 0,
                        'file_exists': False
                    }
    
    def _display_status(self):
        """Отображение статуса в терминале"""
        os.system('clear')
        elapsed = time.time() - self.start_time
        
        print("🚀" + "=" * 70)
        print(f"   SMART Z-HUNT PARALLEL ANALYSIS - {elapsed:.0f}s")
        print("🚀" + "=" * 70)
        print()
        
        for chr_name, data in self.progress_data.items():
            status_icon = "🔬" if data['status'] == 'calculating_zscore' else "✅" if data['status'] == 'completed' else "⏳"
            print(f"{status_icon} {chr_name:8s} | {data['status']:20s} | {data['size_mb']:6.1f} MB")
        
        print()
        print("💡 Файлы НЕ удаляются - все результаты сохраняются!")
        print("📊 Обновление каждые 2 секунды...")
        
    def register_process(self, chr_name, process):
        """Регистрация нового процесса"""
        self.processes[chr_name] = process
        self.progress_data[chr_name] = {
            'status': 'starting',
            'size_mb': 0,
            'file_exists': False
        }
        
    def stop_monitoring(self):
        """Остановка мониторинга"""
        self.monitoring = False

def split_genome_by_chromosome(fasta_file, output_dir):
    """Разделение генома по хромосомам"""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    print(f"📄 Разделяем {fasta_file} по хромосомам...")
    
    chromosomes = {}
    current_chr = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Сохраняем предыдущую хромосому
                if current_chr:
                    chr_file = output_dir / f"{current_chr}.fa"
                    with open(chr_file, 'w') as cf:
                        cf.write(f">{current_chr}\n")
                        cf.write(''.join(current_seq) + '\n')
                    chromosomes[current_chr] = str(chr_file)
                    print(f"   ✅ {current_chr}: {len(''.join(current_seq))} bp")
                
                # Начинаем новую хромосому
                current_chr = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        # Сохраняем последнюю хромосому
        if current_chr:
            chr_file = output_dir / f"{current_chr}.fa"
            with open(chr_file, 'w') as cf:
                cf.write(f">{current_chr}\n")
                cf.write(''.join(current_seq) + '\n')
            chromosomes[current_chr] = str(chr_file)
            print(f"   ✅ {current_chr}: {len(''.join(current_seq))} bp")
    
    print(f"✅ Разделено на {len(chromosomes)} хромосом")
    return chromosomes

def run_zhunt_on_chromosome(chr_name, chr_file, work_dir, monitor, use_rust=False):
    """Запуск Z-Hunt на одной хромосоме"""
    start_time = time.time()
    
    if use_rust:
        # Попытка использовать Rust версию
        cmd = ["./tools/zhunt-rust/target/release/zhunt", chr_file]
    else:
        # Стандартная C версия
        cmd = ["./tools/zhunt/zhunt2", "12", "8", "12", chr_file]
    
    print(f"🔬 Запуск {chr_name} ({'Rust' if use_rust else 'C'})...")
    
    try:
        # Запускаем процесс
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        monitor.register_process(chr_name, process)
        
        # Ждем завершения
        stdout, stderr = process.communicate()
        
        elapsed = time.time() - start_time
        
        if process.returncode == 0:
            # Проверяем результаты
            prob_file = Path(work_dir) / f"{chr_name}.fa.probability"
            zscore_file = Path(work_dir) / f"{chr_name}.fa.Z-SCORE"
            
            prob_size = prob_file.stat().st_size / (1024 * 1024) if prob_file.exists() else 0
            zscore_size = zscore_file.stat().st_size / (1024 * 1024) if zscore_file.exists() else 0
            
            return {
                'chromosome': chr_name,
                'success': True,
                'time': elapsed,
                'prob_size_mb': prob_size,
                'zscore_size_mb': zscore_size,
                'prob_file': str(prob_file),
                'zscore_file': str(zscore_file)
            }
        else:
            return {
                'chromosome': chr_name,
                'success': False,
                'time': elapsed,
                'error': stderr,
                'stdout': stdout
            }
            
    except Exception as e:
        return {
            'chromosome': chr_name,
            'success': False,
            'error': str(e),
            'time': time.time() - start_time
        }

def extract_zdna_results(results, output_file, min_zscore=300, max_zscore=400):
    """Извлечение Z-DNA структур из .probability файлов"""
    print(f"🧬 Извлекаем Z-DNA структуры (Z-score {min_zscore}-{max_zscore})...")
    
    zdna_regions = []
    
    for result in results:
        if not result['success']:
            continue
            
        prob_file = result['prob_file']
        if not os.path.exists(prob_file):
            continue
            
        chr_name = result['chromosome']
        print(f"   📊 Обрабатываем {chr_name}...")
        
        with open(prob_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                    
                try:
                    parts = line.split()
                    if len(parts) >= 3:
                        start_pos = int(parts[0])
                        end_pos = int(parts[1])
                        zscore = float(parts[2])
                        
                        if min_zscore <= zscore <= max_zscore:
                            zdna_regions.append({
                                'chromosome': chr_name,
                                'start': start_pos,
                                'end': end_pos,
                                'zscore': zscore,
                                'length': end_pos - start_pos + 1
                            })
                except (ValueError, IndexError):
                    continue
    
    # Сохраняем результаты
    with open(output_file, 'w') as f:
        f.write("chromosome\tstart\tend\tzscore\tlength\n")
        for region in zdna_regions:
            f.write(f"{region['chromosome']}\t{region['start']}\t{region['end']}\t{region['zscore']}\t{region['length']}\n")
    
    print(f"✅ Найдено {len(zdna_regions)} Z-DNA структур")
    print(f"📄 Результаты сохранены в {output_file}")
    
    return zdna_regions

def main():
    if len(sys.argv) != 3:
        print("Использование: python smart_zhunt_parallel.py <genome.fa> <output.txt>")
        sys.exit(1)
    
    genome_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(genome_file):
        print(f"❌ Файл генома не найден: {genome_file}")
        sys.exit(1)
    
    print("🚀" + "=" * 70)
    print("   SMART Z-HUNT PARALLEL ANALYSIS")
    print("🚀" + "=" * 70)
    print("✅ Файлы НЕ удаляются!")
    print("📊 Мониторинг прогресса в реальном времени!")
    print("⚡ Максимальная загрузка CPU!")
    print()
    
    # Создаем рабочую директорию
    work_dir = Path("z_hunt_results")
    work_dir.mkdir(exist_ok=True)
    
    # Разделяем геном
    start_time = time.time()
    chromosomes = split_genome_by_chromosome(genome_file, work_dir)
    
    # Фильтруем большие хромосомы (> 1MB)
    large_chromosomes = {}
    for chr_name, chr_file in chromosomes.items():
        size = os.path.getsize(chr_file) / (1024 * 1024)
        if size > 1.0:
            large_chromosomes[chr_name] = chr_file
        else:
            print(f"⏭️  Пропускаем маленькую хромосому {chr_name} ({size:.1f}MB)")
    
    print(f"\n🚀 Запуск анализа {len(large_chromosomes)} хромосом...")
    
    # Создаем монитор прогресса
    monitor = ZHuntProgressMonitor(work_dir)
    monitor.start_monitoring()
    
    # Определяем количество параллельных процессов
    max_workers = min(len(large_chromosomes), psutil.cpu_count())
    print(f"💻 Используем {max_workers} параллельных процессов")
    
    # Проверяем доступность Rust версии
    use_rust = os.path.exists("tools/zhunt-rust/target/release/zhunt")
    if use_rust:
        print("🦀 Используем БЫСТРУЮ Rust версию!")
    else:
        print("🔧 Используем стандартную C версию")
    
    # Запускаем анализ в параллель
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Отправляем все задачи
        future_to_chr = {
            executor.submit(run_zhunt_on_chromosome, chr_name, chr_file, work_dir, monitor, use_rust): chr_name
            for chr_name, chr_file in large_chromosomes.items()
        }
        
        # Собираем результаты
        for future in concurrent.futures.as_completed(future_to_chr):
            result = future.result()
            results.append(result)
    
    # Останавливаем мониторинг
    monitor.stop_monitoring()
    
    # Финальные результаты
    successful_results = [r for r in results if r['success']]
    
    if successful_results:
        print("\n🎉 АНАЛИЗ ЗАВЕРШЕН!")
        print(f"✅ Успешно обработано: {len(successful_results)}/{len(large_chromosomes)} хромосом")
        
        # Извлекаем Z-DNA структуры
        zdna_file = output_file.replace('.txt', '_zdna_structures.txt')
        zdna_regions = extract_zdna_results(successful_results, zdna_file)
        
        # Сохраняем сводку результатов
        summary = {
            'timestamp': datetime.now().isoformat(),
            'total_time_seconds': time.time() - start_time,
            'chromosomes_processed': len(successful_results),
            'chromosomes_total': len(large_chromosomes),
            'zdna_regions_found': len(zdna_regions),
            'results': successful_results
        }
        
        summary_file = output_file.replace('.txt', '_summary.json')
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        print(f"📊 Сводка сохранена в {summary_file}")
        print(f"🧬 Z-DNA структуры в {zdna_file}")
        print(f"💾 Все промежуточные файлы сохранены в {work_dir}/")
        print(f"⏱️  Общее время: {summary['total_time_seconds']:.1f} секунд")
        
    else:
        print("\n❌ Анализ не удался!")
        sys.exit(1)

if __name__ == "__main__":
    main() 