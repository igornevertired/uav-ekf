#!/usr/bin/env python3
"""
Главный файл проекта.

Один основной режим (по тезисам):
Navigation EKF (БИНС+ГНСС) → навигационное решение → Parameter EKF (идентификация параметров).

Для исследования влияния точности навигации используется флаг --study, который запускает
серию симуляций для разных классов точности (commercial/tactical/navigation).
"""

import sys
import os
import argparse
import numpy as np

# Добавляем путь к feb в PYTHONPATH для корректных импортов
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Проверяем откуда запущен скрипт
if __package__ is None or __package__ == '':
    # Запущен напрямую из feb/
    from simulation.simulator import IntegratedSimulator
    from config.aircraft_params import AircraftParams
    from config.navigation_params import NavigationParams
else:
    # Запущен как модуль feb.main
    from .simulation.simulator import IntegratedSimulator
    from .config.aircraft_params import AircraftParams
    from .config.navigation_params import NavigationParams


def print_banner():
    """Выводит заголовок программы."""
    print("\n" + "=" * 80)
    print("ИДЕНТИФИКАЦИЯ АЭРОДИНАМИЧЕСКИХ ПАРАМЕТРОВ САМОЛЕТА")
    print("С НАВИГАЦИОННЫМ КОМПЛЕКСОМ БИНС+ГНСС")
    print("=" * 80 + "\n")


def run_single(T: float, accuracy_class: str, visualize: bool, seed: int | None):
    """Один прогон интегрированной системы (БИНС+ГНСС → идентификация)."""
    simulator = IntegratedSimulator(
        dt=0.02,
        accuracy_class=accuracy_class,
        use_navigation=True,
        use_visualization=visualize,
        use_nav_for_identification=True,
        rng_seed=seed,
    )
    success = simulator.run(T=T, profile_type='enhanced')
    if success:
        simulator.print_summary()


def run_study(T: float, visualize: bool, mc: int, seed: int | None):
    """Серия прогонов для разных классов точности навигации."""
    print("\n>>> ИССЛЕДОВАНИЕ ВЛИЯНИЯ ТОЧНОСТИ НАВИГАЦИИ")
    print("    Прогоны: commercial / tactical / navigation\n")

    results = {}
    for accuracy_class in ['commercial', 'tactical', 'navigation']:
        print(f"\n--- Класс: {accuracy_class.upper()} ---")
        mean_errors = []
        nav_pos_rms = []
        nav_vel_rms = []
        for k in range(mc):
            run_seed = None if seed is None else int(seed) + k
            simulator = IntegratedSimulator(
                dt=0.02,
                accuracy_class=accuracy_class,
                use_navigation=True,
                use_visualization=visualize,
                use_nav_for_identification=True,
                rng_seed=run_seed,
            )
            success = simulator.run(T=T, profile_type='enhanced')
            if not success:
                continue
            res = simulator.get_results()
            ident = res.get('identification', {})
            nav = res.get('navigation', {})
            if ident:
                mean_errors.append(float(ident.get('mean_error', np.nan)))
            if nav:
                nav_pos_rms.append(float(nav.get('position_rms', np.nan)))
                nav_vel_rms.append(float(nav.get('velocity_rms', np.nan)))

        results[accuracy_class] = {
            'mean_errors': mean_errors,
            'nav_pos_rms': nav_pos_rms,
            'nav_vel_rms': nav_vel_rms
        }

    print("\n" + "=" * 80)
    print("СРАВНЕНИЕ КЛАССОВ ТОЧНОСТИ")
    print("=" * 80)
    for accuracy_class, res in results.items():
        print(f"\n{accuracy_class.upper()}:")
        errs = np.array(res['mean_errors'], dtype=float)
        if errs.size > 0:
            print(f"  Идентификация - mean_error: {np.nanmean(errs):.2f}% ± {np.nanstd(errs):.2f}% (N={errs.size})")
        pos = np.array(res['nav_pos_rms'], dtype=float)
        vel = np.array(res['nav_vel_rms'], dtype=float)
        if pos.size > 0:
            print(f"  Навигация - RMS позиции: {np.nanmean(pos):.2f} м")
        if vel.size > 0:
            print(f"  Навигация - RMS скорости: {np.nanmean(vel):.3f} м/с")
    print("\n" + "=" * 80)


def main():
    """Главная функция."""
    parser = argparse.ArgumentParser(
        description='Двухуровневая фильтрация: БИНС+ГНСС → идентификация параметров'
    )
    parser.add_argument(
        '--time',
        type=float,
        default=60.0,
        help='Время симуляции в секундах (по умолчанию: 60)'
    )
    parser.add_argument(
        '--accuracy',
        type=str,
        choices=['commercial', 'tactical', 'navigation'],
        default='commercial',
        help='Класс точности навигационного оборудования'
    )
    parser.add_argument(
        '--study',
        action='store_true',
        help='Запустить исследование по 3 классам точности'
    )
    parser.add_argument(
        '--mc',
        type=int,
        default=10,
        help='Число прогонов Monte-Carlo на класс (по умолчанию: 10)'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=123,
        help='Seed для воспроизводимости (по умолчанию: 123). Для случайных прогонов укажите -1'
    )
    parser.add_argument(
        '--visualize',
        action='store_true',
        help='Включить real-time графики (медленнее)'
    )
    
    args = parser.parse_args()
    
    print_banner()
    seed = None if args.seed == -1 else args.seed
    if args.study:
        run_study(T=args.time, visualize=args.visualize, mc=args.mc, seed=seed)
    else:
        run_single(T=args.time, accuracy_class=args.accuracy, visualize=args.visualize, seed=seed)


if __name__ == "__main__":
    main()
