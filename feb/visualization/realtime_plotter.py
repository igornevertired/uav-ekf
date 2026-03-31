"""
Визуализация в реальном времени для идентификации параметров.

Этот модуль создает интерактивные графики, обновляемые во время симуляции.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.gridspec as gridspec


class RealtimePlotter:
    """
    Класс для визуализации результатов идентификации в реальном времени.
    
    График только ДОПОЛНЯЕТСЯ: новые точки добавляются справа, ось X от 0 до текущего времени.
    Никакого скользящего окна — данные накапливаются слева направо.
    
    Отображает:
    - Оценки EKF (синяя линия) и эталон из модели (красная горизонтальная)
    - Ошибки идентификации для каждого параметра
    - Состояния самолета (V, alpha, q, theta)
    - Перегрузку nz
    
    «Эталон (из модели)» — известные значения аэродинамических параметров,
    заложенные в модель самолёта; EKF к ним сходится в процессе идентификации.
    """
    
    def __init__(self, figsize=(16, 10)):
        """
        Инициализация визуализатора.
        
        Args:
            figsize: Размер фигуры
        """
        plt.ion()  # Интерактивный режим
        
        self.fig = plt.figure(figsize=figsize)
        self.title_text = self.fig.suptitle("ИДЕНТИФИКАЦИЯ ПАРАМЕТРОВ В РЕАЛЬНОМ ВРЕМЕНИ - Cessna 172R Skyhawk", 
                         fontsize=14, fontweight='bold')
        
        # Создаем сетку графиков: 3 строки x 3 столбца
        self.gs = gridspec.GridSpec(3, 3, figure=self.fig, hspace=0.35, wspace=0.3)
        
        # Параметры (5 штук) - займем 2 столбца, 2 строки
        self.ax_params = []
        param_names = ['CLα', 'CLδe', 'Cmα', 'Cmq̄', 'Cmδe']
        
        # Размещаем параметры: 3 сверху, 2 снизу
        positions = [(0, 0), (0, 1), (0, 2),  # верхний ряд
                     (1, 0), (1, 1)]           # средний ряд (2 левых)
        
        for idx, (row, col) in enumerate(positions):
            ax = self.fig.add_subplot(self.gs[row, col])
            ax.set_title(param_names[idx], fontweight='bold')
            ax.set_xlabel('Время (с)')
            ax.grid(True, alpha=0.3)
            self.ax_params.append(ax)
        
        # Ошибки параметров (средний ряд, правый столбец)
        self.ax_errors = self.fig.add_subplot(self.gs[1, 2])
        self.ax_errors.set_title('Ошибки параметров (%)', fontweight='bold')
        self.ax_errors.set_xlabel('Время (с)')
        self.ax_errors.set_ylabel('Ошибка (%)')
        self.ax_errors.grid(True, alpha=0.3)
        self.ax_errors.set_ylim(0, 50)
        
        # Состояния самолета (нижний ряд, левый и средний)
        self.ax_states = self.fig.add_subplot(self.gs[2, :2])
        self.ax_states.set_title('Состояния самолета', fontweight='bold')
        self.ax_states.set_xlabel('Время (с)')
        self.ax_states.grid(True, alpha=0.3)
        
        # Перегрузка (нижний ряд, правый)
        self.ax_nz = self.fig.add_subplot(self.gs[2, 2])
        self.ax_nz.set_title('Перегрузка nz', fontweight='bold')
        self.ax_nz.set_xlabel('Время (с)')
        self.ax_nz.set_ylabel('nz (g)')
        self.ax_nz.grid(True, alpha=0.3)
        
        # Линии для обновления
        self.lines_params = []
        self.lines_true = []
        self.lines_errors = []
        
        self.update_interval = 5  # Обновлять каждые N шагов
        self.update_counter = 0
        
        # График просто растёт вправо: ось X от 0 до текущего времени (без скользящего окна)
        plt.show(block=False)
        
        print("Визуализация готова (график дополняется справа, без скользящего окна)")
    
    def update(self, param_ekf, x_true, t_current, true_state_history=None):
        """
        Обновляет графики.
        
        Args:
            param_ekf: EKF для идентификации параметров
            x_true: Истинное состояние (текущее)
            t_current: Текущее время
            true_state_history: Список истинных состояний по шагам [ [V,α,q,θ], ... ]
                               для отображения эталона на графике состояний.
        """
        self.update_counter += 1
        
        # Обновляем только каждый N-й кадр для производительности
        if self.update_counter % self.update_interval != 0:
            return
        
        # Получаем историю
        if len(param_ekf.param_history) == 0:
            return
        
        # Синхронизируем длины массивов
        min_len = min(len(param_ekf.time_history), 
                      len(param_ekf.param_history), 
                      len(param_ekf.state_history))
        
        if min_len == 0:
            return
        
        times = np.array(param_ekf.time_history[:min_len])
        param_hist = np.array(param_ekf.param_history[:min_len])
        state_hist = np.array(param_ekf.state_history[:min_len])
        
        # Истинные значения
        true_params = param_ekf.params.TRUE_PARAMS
        param_names = param_ekf.params.PARAM_NAMES
        
        # Обновляем заголовок
        current_errors = 100 * np.abs(param_hist[-1, :] - true_params) / np.abs(true_params)
        mean_error = np.mean(current_errors)
        self.title_text.set_text(
            f"ИДЕНТИФИКАЦИЯ ПАРАМЕТРОВ В РЕАЛЬНОМ ВРЕМЕНИ - Cessna 172R Skyhawk | "
            f"t = {t_current:.1f}с | Средняя ошибка: {mean_error:.2f}%"
        )
        
        # График только ДОПОЛНЯЕТСЯ: ось X от 0 до текущего времени (никакого сдвига окна)
        x_min = 0
        x_max = max(t_current + 1.0, 2.0)  # справа — текущее время + запас
        
        # 1. Параметры: дописываем данные (рисуем всю накопленную историю)
        for i, ax in enumerate(self.ax_params):
            ax.clear()
            ax.plot(times, param_hist[:, i], 'b-', linewidth=1.5, label='Оценка EKF')
            # Эталон = известные значения из модели самолёта (к ним сходится оценка)
            ax.axhline(true_params[i], color='r', linestyle='--', linewidth=2, 
                       label='Эталон (из модели)')
            ax.axvline(t_current, color='green', linestyle=':', linewidth=1, alpha=0.6)
            ax.set_title(f'{param_names[i]}', fontweight='bold')
            ax.set_xlabel('Время (с)')
            ax.grid(True, alpha=0.3)
            ax.legend(loc='upper right', fontsize=8)
            ax.set_xlim(x_min, x_max)
        
        # 2. Ошибки: дописываем данные
        self.ax_errors.clear()
        errors = 100 * np.abs(param_hist - true_params) / np.abs(true_params)
        for i in range(5):
            self.ax_errors.plot(times, errors[:, i], label=param_names[i], linewidth=1.5)
        self.ax_errors.axvline(t_current, color='green', linestyle=':', linewidth=1, alpha=0.6)
        self.ax_errors.set_title('Ошибки параметров (%)', fontweight='bold')
        self.ax_errors.set_xlabel('Время (с)')
        self.ax_errors.set_ylabel('Ошибка (%)')
        self.ax_errors.legend(loc='upper right', fontsize=8)
        self.ax_errors.grid(True, alpha=0.3)
        self.ax_errors.set_ylim(0, min(100, max(50, np.max(errors) * 1.1 + 5)))
        self.ax_errors.set_xlim(x_min, x_max)
        
        # 3. Состояния: оценка EKF и при наличии — истина (эталон)
        self.ax_states.clear()
        self.ax_states.plot(times, state_hist[:, 0], 'b-', linewidth=1.5, label='Оценка EKF')
        self.ax_states.plot(times, np.rad2deg(state_hist[:, 1]), 'b-', linewidth=1.2)
        self.ax_states.plot(times, np.rad2deg(state_hist[:, 2]), 'b-', linewidth=1.2)
        self.ax_states.plot(times, np.rad2deg(state_hist[:, 3]), 'b-', linewidth=1.2)
        if true_state_history is not None and len(true_state_history) >= min_len:
            true_arr = np.array(true_state_history[:min_len])
            self.ax_states.plot(times, true_arr[:, 0], 'r--', linewidth=1, alpha=0.8, label='Истина (эталон)')
            self.ax_states.plot(times, np.rad2deg(true_arr[:, 1]), 'r--', linewidth=0.8, alpha=0.8)
            self.ax_states.plot(times, np.rad2deg(true_arr[:, 2]), 'r--', linewidth=0.8, alpha=0.8)
            self.ax_states.plot(times, np.rad2deg(true_arr[:, 3]), 'r--', linewidth=0.8, alpha=0.8)
        self.ax_states.axvline(t_current, color='green', linestyle=':', linewidth=1, alpha=0.6)
        self.ax_states.set_title('Состояния самолета', fontweight='bold')
        self.ax_states.set_xlabel('Время (с)')
        self.ax_states.legend(loc='upper right', fontsize=8)
        self.ax_states.grid(True, alpha=0.3)
        self.ax_states.set_xlim(x_min, x_max)
        
        # 4. Перегрузка nz: дописываем данные
        if len(param_ekf.measurement_history) > 0:
            meas_len = min(len(param_ekf.measurement_history), len(times))
            if meas_len > 0:
                meas_hist = np.array(param_ekf.measurement_history[:meas_len])
                times_meas = times[:meas_len]
                self.ax_nz.clear()
                self.ax_nz.plot(times_meas, meas_hist[:, 4], 'g-', linewidth=1.5)
                self.ax_nz.axvline(t_current, color='green', linestyle=':', linewidth=1, alpha=0.6)
                self.ax_nz.set_title('Перегрузка nz', fontweight='bold')
                self.ax_nz.set_xlabel('Время (с)')
                self.ax_nz.set_ylabel('nz (g)')
                self.ax_nz.grid(True, alpha=0.3)
                self.ax_nz.set_xlim(x_min, x_max)
        
        # Финальные ошибки (текстом)
        if t_current > 50:  # Показываем под конец симуляции
            final_errors = errors[-1, :]
            mean_error = np.mean(final_errors)
            
            text_str = f'Средняя ошибка: {mean_error:.2f}%\n'
            text_str += f'Макс ошибка: {np.max(final_errors):.2f}%'
            
            self.ax_errors.text(0.02, 0.98, text_str,
                              transform=self.ax_errors.transAxes,
                              verticalalignment='top',
                              bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                              fontsize=9)
        
        # Обновляем отображение графиков
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()
        plt.pause(0.01)  # Небольшая пауза для плавного обновления
    
    def close(self):
        """Закрывает окно визуализации."""
        plt.ioff()
        # plt.close(self.fig)  # Оставляем окно открытым для анализа
    
    def show_final(self):
        """Показывает финальный результат."""
        plt.ioff()
        plt.show()
