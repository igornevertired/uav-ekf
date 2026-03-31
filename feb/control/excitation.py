"""
Профили управляющих сигналов для идентификации параметров.

Этот модуль содержит различные профили возбуждающих сигналов (input signals)
для эффективной идентификации аэродинамических параметров.

Принцип: богатое возбуждение (persistent excitation) необходимо для того,
чтобы фильтр Калмана мог различить влияние разных параметров на динамику.
"""

import numpy as np
from typing import Callable
from ..config.aircraft_params import AircraftParams


class ExcitationProfile:
    """
    Класс для генерации профилей возбуждающих управляющих сигналов.
    
    Различные типы сигналов:
    - Ступенчатые (doublets)
    - Синусоидальные (multi-sine)
    - Комбинированные
    - 3-2-1-1 маневры
    """
    
    def __init__(self, delta_e_trim: float = 0.0):
        """
        Инициализация генератора профилей.
        
        Args:
            delta_e_trim: Балансировочное отклонение руля высоты (рад)
        """
        self.delta_e_trim = delta_e_trim
        self.current_profile = 'enhanced'  # По умолчанию
    
    def enhanced_profile(self, t: float) -> float:
        """
        Улучшенный профиль с богатым возбуждением.
        
        Сочетает ступенчатые изменения и синусоиды разных частот.
        
        Args:
            t: Время (сек)
        
        Returns:
            Отклонение руля высоты (рад)
        """
        base = self.delta_e_trim
        
        if t < 2:
            return base
        elif t < 4:
            return base + np.deg2rad(-2.0)
        elif t < 6:
            return base + np.deg2rad(2.0)
        elif t < 7:
            return base
        elif t < 8:
            return base + np.deg2rad(-3.0)
        elif t < 10:
            return base + np.deg2rad(3.0)
        elif t < 11:
            return base
        elif t < 12:
            return base + np.deg2rad(-3.0)
        elif t < 14:
            return base + np.deg2rad(3.0)
        else:
            # Многочастотная синусоида
            return base + np.deg2rad(
                2.2 * np.sin(0.5 * t) +
                1.1 * np.sin(1.5 * t) +
                0.4 * np.sin(3.0 * t)
            )
    
    def doublet_profile(self, t: float, amplitude: float = 2.0, duration: float = 1.0) -> float:
        """
        Профиль с дублетами (doublet).
        
        Args:
            t: Время (сек)
            amplitude: Амплитуда (градусы)
            duration: Длительность одного импульса (сек)
        
        Returns:
            Отклонение руля (рад)
        """
        base = self.delta_e_trim
        amp = np.deg2rad(amplitude)
        
        # Начинаем маневры после 2 секунд
        if t < 2:
            return base
        
        t_rel = t - 2
        period = 2 * duration + 2.0  # Два импульса + пауза
        
        phase = t_rel % period
        
        if phase < duration:
            return base + amp
        elif phase < 2 * duration:
            return base - amp
        else:
            return base
    
    def multisine_profile(self, t: float, frequencies: list = None, amplitudes: list = None) -> float:
        """
        Многочастотный синусоидальный сигнал.
        
        Args:
            t: Время (сек)
            frequencies: Список частот (Гц)
            amplitudes: Список амплитуд (градусы)
        
        Returns:
            Отклонение руля (рад)
        """
        if frequencies is None:
            frequencies = [0.1, 0.3, 0.7, 1.2]  # Гц
        if amplitudes is None:
            amplitudes = [2.0, 1.5, 1.0, 0.5]  # градусы
        
        base = self.delta_e_trim
        signal = 0.0
        
        for freq, amp in zip(frequencies, amplitudes):
            signal += np.deg2rad(amp) * np.sin(2 * np.pi * freq * t)
        
        return base + signal
    
    def step_profile(self, t: float, step_time: float = 5.0, step_amplitude: float = 2.0) -> float:
        """
        Ступенчатое изменение.
        
        Args:
            t: Время (сек)
            step_time: Время начала ступеньки (сек)
            step_amplitude: Амплитуда ступеньки (градусы)
        
        Returns:
            Отклонение руля (рад)
        """
        base = self.delta_e_trim
        
        if t < step_time:
            return base
        else:
            return base + np.deg2rad(step_amplitude)
    
    def chirp_profile(self, t: float, f0: float = 0.1, f1: float = 2.0, duration: float = 30.0) -> float:
        """
        Чирп-сигнал (sweep) - синусоида с изменяющейся частотой.
        
        Args:
            t: Время (сек)
            f0: Начальная частота (Гц)
            f1: Конечная частота (Гц)
            duration: Длительность sweep (сек)
        
        Returns:
            Отклонение руля (рад)
        """
        base = self.delta_e_trim
        
        if t < 2 or t > duration + 2:
            return base
        
        t_rel = t - 2
        # Линейное изменение частоты
        k = (f1 - f0) / duration
        instantaneous_freq = f0 + k * t_rel
        phase = 2 * np.pi * (f0 * t_rel + 0.5 * k * t_rel**2)
        
        amplitude = np.deg2rad(2.0)
        return base + amplitude * np.sin(phase)
    
    def three_two_one_one(self, t: float) -> float:
        """
        Маневр 3-2-1-1 (стандартный идентификационный маневр).
        
        Returns:
            Отклонение руля (рад)
        """
        base = self.delta_e_trim
        
        if t < 2:
            return base
        elif t < 5:  # 3 секунды
            return base + np.deg2rad(3.0)
        elif t < 7:  # 2 секунды
            return base - np.deg2rad(3.0)
        elif t < 8:  # 1 секунда
            return base + np.deg2rad(2.0)
        elif t < 9:  # 1 секунда
            return base - np.deg2rad(2.0)
        else:
            return base
    
    def get_control(self, t: float, profile_type: str = None) -> float:
        """
        Получить управляющий сигнал для заданного времени.
        
        Args:
            t: Время (сек)
            profile_type: Тип профиля ('enhanced', 'doublet', 'multisine', 
                         'step', 'chirp', '3-2-1-1', 'none')
                         Если None, используется текущий профиль
        
        Returns:
            Отклонение руля высоты (рад)
        """
        if profile_type is None:
            profile_type = self.current_profile
        
        profile_map = {
            'enhanced': self.enhanced_profile,
            'doublet': self.doublet_profile,
            'multisine': self.multisine_profile,
            'step': self.step_profile,
            'chirp': self.chirp_profile,
            '3-2-1-1': self.three_two_one_one,
            'none': lambda t: self.delta_e_trim
        }
        
        if profile_type not in profile_map:
            raise ValueError(f"Unknown profile type: {profile_type}")
        
        delta = float(profile_map[profile_type](t))
        # Физическое ограничение руля высоты
        lim = float(AircraftParams.ELEVATOR_LIMIT)
        return float(np.clip(delta, -lim, lim))
    
    def set_profile(self, profile_type: str):
        """
        Устанавливает текущий профиль.
        
        Args:
            profile_type: Тип профиля
        """
        self.current_profile = profile_type
    
    def set_trim(self, delta_e_trim: float):
        """
        Устанавливает балансировочное отклонение.
        
        Args:
            delta_e_trim: Балансировочное отклонение (рад)
        """
        self.delta_e_trim = delta_e_trim
    
    @staticmethod
    def create_custom_profile(func: Callable[[float], float]) -> Callable:
        """
        Создает пользовательский профиль из функции.
        
        Args:
            func: Функция f(t) -> delta_e
        
        Returns:
            Callable профиль
        """
        return func
    
    def get_profile_description(self, profile_type: str = None) -> str:
        """
        Возвращает описание профиля.
        
        Args:
            profile_type: Тип профиля
        
        Returns:
            Текстовое описание
        """
        if profile_type is None:
            profile_type = self.current_profile
        
        descriptions = {
            'enhanced': 'Комбинированный профиль: ступеньки + многочастотная синусоида',
            'doublet': 'Дублеты (положительный и отрицательный импульсы)',
            'multisine': 'Многочастотный синусоидальный сигнал',
            'step': 'Ступенчатое изменение',
            'chirp': 'Чирп (sweep) - синусоида с изменяющейся частотой',
            '3-2-1-1': 'Стандартный маневр 3-2-1-1 секунды',
            'none': 'Без возбуждения (только балансировка)'
        }
        
        return descriptions.get(profile_type, 'Неизвестный профиль')
    
    def print_info(self):
        """Выводит информацию о доступных профилях."""
        print("=" * 80)
        print("ДОСТУПНЫЕ ПРОФИЛИ УПРАВЛЕНИЯ")
        print("=" * 80)
        print(f"Текущий профиль: {self.current_profile}")
        print(f"Балансировочное δe: {np.rad2deg(self.delta_e_trim):.2f}°")
        print("-" * 80)
        
        for profile_type in ['enhanced', 'doublet', 'multisine', 'step', 'chirp', '3-2-1-1', 'none']:
            desc = self.get_profile_description(profile_type)
            print(f"  {profile_type:12} - {desc}")
        
        print("=" * 80)
