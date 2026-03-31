"""
Анализатор результатов и проведение исследований.
"""

import numpy as np
from typing import Dict, List


class ResultsAnalyzer:
    """Класс для анализа результатов симуляций."""
    
    def __init__(self):
        self.results = []
    
    def analyze_identification(self, param_ekf) -> Dict:
        """
        Анализирует результаты идентификации параметров.
        
        Returns:
            Словарь со статистикой
        """
        if len(param_ekf.param_history) == 0:
            return {}
        
        final_params = param_ekf.param_history[-1]
        true_params = param_ekf.params.TRUE_PARAMS
        
        errors = 100 * np.abs(final_params - true_params) / np.abs(true_params)
        
        return {
            'final_params': final_params,
            'errors_percent': errors,
            'mean_error': np.mean(errors),
            'max_error': np.max(errors),
            'convergence_time': param_ekf.time_history[-1] if len(param_ekf.time_history) > 0 else 0
        }
    
    def analyze_navigation(self, nav_ekf) -> Dict:
        """
        Анализирует точность навигационного решения.
        
        Returns:
            Словарь со статистикой навигации
        """
        solution = nav_ekf.get_navigation_solution()
        
        return {
            'position_std': solution['position_std'],
            'velocity_std': solution['velocity_std'],
            'attitude_std': solution['attitude_std'],
            'position_rms': np.linalg.norm(solution['position_std']),
            'velocity_rms': np.linalg.norm(solution['velocity_std'])
        }
    
    def print_summary(self, ident_results: Dict, nav_results: Dict):
        """Выводит сводку результатов."""
        print("=" * 80)
        print("СВОДКА РЕЗУЛЬТАТОВ")
        print("=" * 80)
        
        if ident_results:
            print("\nИДЕНТИФИКАЦИЯ ПАРАМЕТРОВ:")
            print(f"  Средняя ошибка: {ident_results['mean_error']:.2f}%")
            print(f"  Максимальная ошибка: {ident_results['max_error']:.2f}%")
            print(f"  Время симуляции: {ident_results['convergence_time']:.1f} с")
        
        if nav_results:
            print("\nНАВИГАЦИОННОЕ РЕШЕНИЕ:")
            print(f"  RMS позиции: {nav_results['position_rms']:.2f} м")
            print(f"  RMS скорости: {nav_results['velocity_rms']:.3f} м/с")
        
        print("=" * 80)
