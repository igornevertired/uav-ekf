"""
Графики для анализа результатов исследования.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class AnalysisPlots:
    """Класс для построения аналитических графиков."""
    
    @staticmethod
    def plot_final_parameters(param_ekf):
        """Строит финальные графики параметров."""
        if len(param_ekf.param_history) == 0:
            return
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle('ФИНАЛЬНЫЙ АНАЛИЗ ИДЕНТИФИКАЦИИ', fontsize=14, fontweight='bold')
        
        time_data = np.array(param_ekf.time_history)
        params = np.array(param_ekf.param_history)
        
        # Графики параметров
        for i in range(5):
            ax = axes.flat[i]
            ax.plot(time_data, params[:, i], 'b-', linewidth=1.5, label='Оценка')
            ax.axhline(y=param_ekf.params.TRUE_PARAMS[i], color='r', linestyle='--', label='Истинное')
            ax.set_xlabel('Время, с')
            ax.set_ylabel(param_ekf.params.PARAM_NAMES[i])
            ax.set_title(f'{param_ekf.params.PARAM_NAMES[i]} - Идентификация')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
    
    @staticmethod
    def plot_navigation_trajectory(nav_ekf, ins_model):
        """Строит 3D траекторию полета."""
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Здесь можно добавить траекторию из истории
        ax.set_xlabel('X, м')
        ax.set_ylabel('Y, м')
        ax.set_zlabel('Z, м')
        ax.set_title('ТРАЕКТОРИЯ ПОЛЕТА')
        plt.show()
    
    @staticmethod
    def plot_navigation_errors(nav_ekf):
        """Строит графики ошибок навигации."""
        fig, axes = plt.subplots(3, 1, figsize=(12, 10))
        fig.suptitle('ОШИБКИ НАВИГАЦИОННОГО РЕШЕНИЯ', fontsize=14, fontweight='bold')
        
        time_data = np.array(nav_ekf.time_history) if len(nav_ekf.time_history) > 0 else [0]
        
        axes[0].set_ylabel('Ошибка позиции, м')
        axes[0].set_title('ПОЗИЦИЯ')
        axes[0].grid(True)
        
        axes[1].set_ylabel('Ошибка скорости, м/с')
        axes[1].set_title('СКОРОСТЬ')
        axes[1].grid(True)
        
        axes[2].set_xlabel('Время, с')
        axes[2].set_ylabel('Ошибка ориентации, °')
        axes[2].set_title('ОРИЕНТАЦИЯ')
        axes[2].grid(True)
        
        plt.tight_layout()
        plt.show()
