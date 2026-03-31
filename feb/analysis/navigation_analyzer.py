"""
Анализатор навигационных данных в реальном времени.

Предоставляет функции для анализа и визуализации данных БИНС и ГНСС.
"""

import numpy as np
import matplotlib.pyplot as plt


class NavigationDataAnalyzer:
    """
    Класс для анализа данных навигационной системы.
    
    Позволяет:
    - Извлекать данные БИНС и ГНСС из истории симуляции
    - Анализировать дрейф БИНС
    - Оценивать эффективность коррекции
    - Визуализировать траектории
    """
    
    @staticmethod
    def get_ins_data(history: dict) -> dict:
        """
        Извлекает данные БИНС из истории.
        
        Args:
            history: История симуляции
        
        Returns:
            dict: Данные БИНС по времени
        """
        if history['ins_raw'] is None or len(history['ins_raw']) == 0:
            return None
        
        times = np.array(history['time'])
        
        positions = np.array([h['position'] for h in history['ins_raw']])
        velocities = np.array([h['velocity'] for h in history['ins_raw']])
        attitudes = np.array([h['attitude'] for h in history['ins_raw']])
        gyro_biases = np.array([h['gyro_bias'] for h in history['ins_raw']])
        accel_biases = np.array([h['accel_bias'] for h in history['ins_raw']])
        
        return {
            'time': times,
            'position': positions,
            'velocity': velocities,
            'attitude': attitudes,
            'gyro_bias': gyro_biases,
            'accel_bias': accel_biases
        }
    
    @staticmethod
    def get_gnss_data(history: dict) -> dict:
        """
        Извлекает данные ГНСС из истории.
        
        Args:
            history: История симуляции
        
        Returns:
            dict: Данные ГНСС по времени
        """
        if history['gnss_measurement'] is None or len(history['gnss_measurement']) == 0:
            return None
        
        times = []
        positions = []
        velocities = []
        
        for i, meas in enumerate(history['gnss_measurement']):
            if meas['available']:
                times.append(history['time'][i])
                positions.append(meas['position'])
                velocities.append(meas['velocity'])
        
        return {
            'time': np.array(times),
            'position': np.array(positions),
            'velocity': np.array(velocities)
        }
    
    @staticmethod
    def get_corrected_data(history: dict) -> dict:
        """
        Извлекает скорректированные данные из истории.
        
        Args:
            history: История симуляции
        
        Returns:
            dict: Скорректированные данные по времени
        """
        if history['nav_corrected'] is None or len(history['nav_corrected']) == 0:
            return None
        
        times = []
        positions = []
        velocities = []
        attitudes = []
        
        for i, corr in enumerate(history['nav_corrected']):
            if corr is not None:
                times.append(history['time'][i])
                positions.append(corr['position'])
                velocities.append(corr['velocity'])
                attitudes.append(corr['attitude'])
        
        return {
            'time': np.array(times),
            'position': np.array(positions),
            'velocity': np.array(velocities),
            'attitude': np.array(attitudes)
        }
    
    @staticmethod
    def compute_drift(ins_data: dict, true_trajectory: dict = None) -> dict:
        """
        Вычисляет дрейф БИНС.
        
        Args:
            ins_data: Данные БИНС
            true_trajectory: Истинная траектория (если известна)
        
        Returns:
            dict: Статистика дрейфа
        """
        if ins_data is None:
            return None
        
        times = ins_data['time']
        positions = ins_data['position']
        
        # Начальная позиция
        pos_0 = positions[0]
        
        # Дрейф относительно начала
        drift = np.linalg.norm(positions - pos_0, axis=1)
        
        # Скорость дрейфа
        drift_rate = np.gradient(drift, times)
        
        result = {
            'times': times,
            'drift_magnitude': drift,
            'drift_rate': drift_rate,
            'final_drift': drift[-1],
            'max_drift_rate': np.max(np.abs(drift_rate))
        }
        
        if true_trajectory is not None:
            # Ошибка относительно истины
            true_pos = true_trajectory['position']
            errors = np.linalg.norm(positions - true_pos, axis=1)
            result['position_error'] = errors
            result['rms_error'] = np.sqrt(np.mean(errors**2))
        
        return result
    
    @staticmethod
    def plot_trajectories(history: dict, save_path: str = None):
        """
        Визуализирует траектории БИНС, ГНСС и скорректированные.
        
        Args:
            history: История симуляции
            save_path: Путь для сохранения графика
        """
        ins_data = NavigationDataAnalyzer.get_ins_data(history)
        gnss_data = NavigationDataAnalyzer.get_gnss_data(history)
        corr_data = NavigationDataAnalyzer.get_corrected_data(history)
        
        if ins_data is None:
            print("Нет данных навигации для визуализации")
            return
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # 1. Траектория в 3D (проекция X-Z)
        ax = axes[0, 0]
        ax.plot(ins_data['position'][:, 0], ins_data['position'][:, 2], 
                'r-', linewidth=1.5, label='БИНС (raw)', alpha=0.7)
        
        if gnss_data is not None:
            ax.scatter(gnss_data['position'][:, 0], gnss_data['position'][:, 2],
                      c='g', s=20, label='ГНСС', alpha=0.5)
        
        if corr_data is not None:
            ax.plot(corr_data['position'][:, 0], corr_data['position'][:, 2],
                   'b-', linewidth=2, label='Скорректированная', alpha=0.8)
        
        ax.set_xlabel('X (м)')
        ax.set_ylabel('Z - высота (м)')
        ax.set_title('Траектория полета (вид сбоку)', fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 2. Дрейф позиции
        ax = axes[0, 1]
        drift = np.linalg.norm(ins_data['position'] - ins_data['position'][0], axis=1)
        ax.plot(ins_data['time'], drift / 1000, 'r-', linewidth=2, label='БИНС дрейф')
        
        if corr_data is not None:
            corr_drift = np.linalg.norm(corr_data['position'] - ins_data['position'][0], axis=1)
            ax.plot(corr_data['time'], corr_drift / 1000, 'b-', linewidth=2, label='После коррекции')
        
        ax.set_xlabel('Время (с)')
        ax.set_ylabel('Дрейф (км)')
        ax.set_title('Дрейф позиции во времени', fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 3. Bias датчиков
        ax = axes[1, 0]
        gyro_bias_mag = np.linalg.norm(ins_data['gyro_bias'], axis=1)
        accel_bias_mag = np.linalg.norm(ins_data['accel_bias'], axis=1)
        
        ax.plot(ins_data['time'], np.rad2deg(gyro_bias_mag) * 3600, 'r-', 
                linewidth=1.5, label='Gyro bias (deg/hr)')
        ax.plot(ins_data['time'], accel_bias_mag * 1000, 'b-',
                linewidth=1.5, label='Accel bias (mg)')
        
        ax.set_xlabel('Время (с)')
        ax.set_ylabel('Величина bias')
        ax.set_title('Дрейф смещений датчиков', fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 4. Статистика
        ax = axes[1, 1]
        ax.axis('off')
        
        stats_text = "СТАТИСТИКА НАВИГАЦИИ\n"
        stats_text += "=" * 40 + "\n\n"
        
        final_drift = drift[-1]
        stats_text += f"Финальный дрейф БИНС: {final_drift:.1f} м\n"
        
        if corr_data is not None:
            final_corr_drift = corr_drift[-1]
            improvement = (1 - final_corr_drift / final_drift) * 100
            stats_text += f"После коррекции: {final_corr_drift:.1f} м\n"
            stats_text += f"Улучшение: {improvement:.1f}%\n\n"
        
        if gnss_data is not None:
            stats_text += f"ГНСС измерений: {len(gnss_data['time'])}\n"
            stats_text += f"Частота: {len(gnss_data['time']) / ins_data['time'][-1]:.1f} Гц\n\n"
        
        stats_text += f"Финальный gyro bias: {np.rad2deg(gyro_bias_mag[-1]) * 3600:.2f} °/hr\n"
        stats_text += f"Финальный accel bias: {accel_bias_mag[-1] * 1000:.2f} mg"
        
        ax.text(0.1, 0.9, stats_text, transform=ax.transAxes,
               fontsize=10, verticalalignment='top', family='monospace',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.suptitle('АНАЛИЗ НАВИГАЦИОННОЙ СИСТЕМЫ', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=150)
            print(f"График сохранен: {save_path}")
        
        return fig
    
    @staticmethod
    def print_realtime_summary(ins_data: dict, gnss_data: dict, corr_data: dict):
        """
        Выводит сводку навигационных данных.
        
        Args:
            ins_data: Данные БИНС
            gnss_data: Данные ГНСС
            corr_data: Скорректированные данные
        """
        print("\n" + "=" * 80)
        print("СВОДКА НАВИГАЦИОННЫХ ДАННЫХ")
        print("=" * 80)
        
        if ins_data is not None:
            print(f"\nБИНС:")
            print(f"  Время работы: {ins_data['time'][-1]:.1f} с")
            print(f"  Точек данных: {len(ins_data['time'])}")
            
            drift = np.linalg.norm(ins_data['position'] - ins_data['position'][0], axis=1)
            print(f"  Начальный дрейф: {drift[0]:.2f} м")
            print(f"  Финальный дрейф: {drift[-1]:.1f} м")
            print(f"  Скорость дрейфа: {drift[-1] / ins_data['time'][-1]:.2f} м/с")
        
        if gnss_data is not None:
            print(f"\nГНСС:")
            print(f"  Измерений: {len(gnss_data['time'])}")
            print(f"  Частота: {len(gnss_data['time']) / ins_data['time'][-1]:.1f} Гц")
        
        if corr_data is not None:
            print(f"\nСКОРРЕКТИРОВАННЫЕ ДАННЫЕ:")
            print(f"  Точек данных: {len(corr_data['time'])}")
            
            corr_drift = np.linalg.norm(corr_data['position'] - ins_data['position'][0], axis=1)
            print(f"  Финальный дрейф: {corr_drift[-1]:.1f} м")
            
            if ins_data is not None:
                improvement = (1 - corr_drift[-1] / drift[-1]) * 100
                print(f"  Улучшение: {improvement:.1f}%")
        
        print("=" * 80)
