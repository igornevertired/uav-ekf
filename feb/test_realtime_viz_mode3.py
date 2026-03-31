#!/usr/bin/env python3
"""
Real-time графики для режима 3 (идентификация по навигационному решению).

Запуск:
  python test_realtime_viz_mode3.py
"""

import sys
import os

# Добавляем путь к university (чтобы импортировать feb как пакет)
university_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, university_dir)

from feb.simulation.simulator import IntegratedSimulator


def main():
    print("=" * 80)
    print("REAL-TIME ВИЗУАЛИЗАЦИЯ — РЕЖИМ 3 (ПО НАВИГАЦИОННОМУ РЕШЕНИЮ)")
    print("=" * 80)
    print("Графики такие же, как в режиме 2, но измерения для идентификации берутся")
    print("из навигационного решения (БИНС+ГНСС+Navigation EKF).")
    print("-" * 80)

    simulator = IntegratedSimulator(
        dt=0.02,
        accuracy_class="commercial",
        use_navigation=True,
        use_visualization=True,
        use_nav_for_identification=True,
    )

    success = simulator.run(T=60.0, profile_type="enhanced")
    if success:
        simulator.print_summary()
        input("\nНажмите Enter чтобы закрыть окно...")


if __name__ == "__main__":
    main()

