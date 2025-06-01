#!/usr/bin/env python3
"""
🛑 ПЕРЕЗАПУСК УМНОГО Z-HUNT АНАЛИЗА
- Останавливает все текущие процессы
- Запускает новый умный параллельный анализ
"""

import subprocess
import os
import time
import sys

def stop_all_zhunt_processes():
    """Остановка всех процессов Z-Hunt"""
    print("🛑 Останавливаем все процессы Z-Hunt...")
    
    try:
        # Ищем все процессы zhunt
        result = subprocess.run(['pgrep', '-f', 'zhunt'], capture_output=True, text=True)
        
        if result.stdout.strip():
            pids = result.stdout.strip().split('\n')
            print(f"   Найдено {len(pids)} процессов Z-Hunt")
            
            # Останавливаем все процессы
            for pid in pids:
                try:
                    subprocess.run(['kill', '-TERM', pid], check=True)
                    print(f"   ✅ Остановлен процесс {pid}")
                except subprocess.CalledProcessError:
                    print(f"   ⚠️  Не удалось остановить процесс {pid}")
            
            # Ждем завершения процессов
            print("   ⏳ Ждем завершения процессов...")
            time.sleep(3)
            
            # Принудительно убиваем оставшиеся процессы
            result = subprocess.run(['pgrep', '-f', 'zhunt'], capture_output=True, text=True)
            if result.stdout.strip():
                remaining_pids = result.stdout.strip().split('\n')
                for pid in remaining_pids:
                    try:
                        subprocess.run(['kill', '-KILL', pid], check=True)
                        print(f"   💥 Принудительно остановлен процесс {pid}")
                    except subprocess.CalledProcessError:
                        pass
        else:
            print("   ✅ Процессы Z-Hunt не найдены")
            
    except Exception as e:
        print(f"   ⚠️  Ошибка при остановке процессов: {e}")

def install_dependencies():
    """Установка недостающих зависимостей"""
    print("📦 Проверяем зависимости...")
    
    try:
        import psutil
        print("   ✅ psutil установлен")
    except ImportError:
        print("   📦 Устанавливаем psutil...")
        subprocess.run([sys.executable, '-m', 'pip', 'install', 'psutil'], check=True)

def main():
    print("🚀" + "=" * 60)
    print("   ПЕРЕЗАПУСК УМНОГО Z-HUNT АНАЛИЗА")
    print("🚀" + "=" * 60)
    print()
    
    # Остановка текущих процессов
    stop_all_zhunt_processes()
    
    # Установка зависимостей
    install_dependencies()
    
    # Проверка файлов
    genome_file = "data/genome/dm6.fa"
    output_file = "results/smart_zhunt_results.txt"
    
    if not os.path.exists(genome_file):
        print(f"❌ Файл генома не найден: {genome_file}")
        sys.exit(1)
    
    print(f"✅ Файл генома найден: {genome_file}")
    print(f"📄 Результаты будут сохранены в: {output_file}")
    print()
    
    # Запуск умного анализа
    print("🚀 ЗАПУСКАЕМ УМНЫЙ ПАРАЛЛЕЛЬНЫЙ Z-HUNT!")
    print("📊 Мониторинг прогресса будет показан в реальном времени")
    print("💾 Все файлы будут сохранены в директории z_hunt_results/")
    print()
    
    # Запускаем новый скрипт
    cmd = [sys.executable, "scripts/smart_zhunt_parallel.py", genome_file, output_file]
    
    try:
        subprocess.run(cmd, check=True)
    except KeyboardInterrupt:
        print("\n⏹️  Анализ прерван пользователем")
        stop_all_zhunt_processes()
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Ошибка при запуске анализа: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 