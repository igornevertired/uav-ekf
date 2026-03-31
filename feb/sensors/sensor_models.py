"""
Модели датчиков с характерными ошибками и шумами.

Этот модуль содержит модели бортовых датчиков самолета:
- Датчик воздушной скорости (Pitot-static system)
- Датчик угла атаки (AoA vane)
- Гироскопы (угловая скорость)
- Акселерометры (перегрузка)
"""

import numpy as np
from ..config.aircraft_params import AircraftParams


class SensorModels:
    """
    Класс для моделирования измерений бортовых датчиков с шумами и ошибками.
    
    Single Responsibility: только моделирование датчиков.
    """
    
    def __init__(self, params: AircraftParams = None):
        """
        Инициализация моделей датчиков.
        
        Args:
            params: Параметры самолета
        """
        self.params = params if params is not None else AircraftParams
        
        # Уровни шума датчиков (СКО)
        self.noise_levels = {
            'velocity': 0.5,        # м/с
            'aoa': 0.001,          # рад (~0.06°)
            'angular_rate': 0.001, # рад/с (~0.06°/с)
            'attitude': 0.001,     # рад (~0.06°)
            'load_factor': 0.005   # безразмерная (0.005 g)
        }
        
        # Смещения (bias) датчиков
        self.biases = {
            'velocity': 0.0,
            'aoa': 0.0,
            'angular_rate': 0.0,
            'attitude': 0.0,
            'load_factor': 0.0
        }
    
    def generate_sensor_data(
        self,
        V: float,
        alpha: float,
        q: float,
        theta: float,
        delta_e: float,
        identified_params: np.ndarray,
        add_noise: bool = True
    ) -> np.ndarray:
        """
        Генерирует зашумленные измерения всех датчиков.
        
        Args:
            V: Истинная воздушная скорость (м/с)
            alpha: Истинный угол атаки (рад)
            q: Истинная угловая скорость тангажа (рад/с)
            theta: Истинный угол тангажа (рад)
            delta_e: Отклонение руля высоты (рад)
            identified_params: Аэродинамические параметры для расчета nz
            add_noise: Добавлять ли шум
        
        Returns:
            Вектор измерений [V_m, alpha_m, q_m, theta_m, nz_m]
        """
        # 1. Воздушная скорость (с шумом)
        V_measured = self._measure_velocity(V, add_noise)
        
        # 2. Угол атаки (с динамической ошибкой датчика)
        alpha_measured = self._measure_aoa(V, alpha, q, add_noise)
        
        # 3. Угловая скорость тангажа
        q_measured = self._measure_angular_rate(q, add_noise)
        
        # 4. Угол тангажа
        theta_measured = self._measure_attitude(theta, add_noise)
        
        # 5. Нормальная перегрузка
        nz_measured = self._measure_load_factor(V, alpha, q, theta, delta_e, identified_params, add_noise)
        
        return np.array([V_measured, alpha_measured, q_measured, theta_measured, nz_measured])
    
    def _measure_velocity(self, V_true: float, add_noise: bool = True) -> float:
        """
        Модель датчика воздушной скорости (Pitot-static).
        
        Args:
            V_true: Истинная скорость (м/с)
            add_noise: Добавлять ли шум
        
        Returns:
            Измеренная скорость (м/с)
        """
        V_measured = V_true + self.biases['velocity']
        
        if add_noise:
            noise = np.random.randn() * self.noise_levels['velocity']
            V_measured += noise
        
        return V_measured
    
    def _measure_aoa(self, V: float, alpha_true: float, q: float, add_noise: bool = True) -> float:
        """
        Модель датчика угла атаки (AoA vane).
        
        Включает динамическую ошибку: α_measured = α_true + d_aoa * q/V
        где d_aoa - коэффициент динамической ошибки (обусловлен расположением датчика)
        
        Args:
            V: Воздушная скорость (м/с)
            alpha_true: Истинный угол атаки (рад)
            q: Угловая скорость (рад/с)
            add_noise: Добавлять ли шум
        
        Returns:
            Измеренный угол атаки (рад)
        """
        # Динамическая ошибка из-за расположения датчика относительно ЦМ
        d_aoa = self.params.AOA_SENSOR_DYNAMIC_ERROR
        dynamic_error = d_aoa * q / max(V, 1e-2)
        
        alpha_measured = alpha_true + dynamic_error + self.biases['aoa']
        
        if add_noise:
            noise = np.random.randn() * self.noise_levels['aoa']
            alpha_measured += noise
        
        return alpha_measured
    
    def _measure_angular_rate(self, q_true: float, add_noise: bool = True) -> float:
        """
        Модель гироскопа (угловая скорость тангажа).
        
        Args:
            q_true: Истинная угловая скорость (рад/с)
            add_noise: Добавлять ли шум
        
        Returns:
            Измеренная угловая скорость (рад/с)
        """
        q_measured = q_true + self.biases['angular_rate']
        
        if add_noise:
            noise = np.random.randn() * self.noise_levels['angular_rate']
            q_measured += noise
        
        return q_measured
    
    def _measure_attitude(self, theta_true: float, add_noise: bool = True) -> float:
        """
        Модель датчика угла тангажа (из БИНС или AHRS).
        
        Args:
            theta_true: Истинный угол тангажа (рад)
            add_noise: Добавлять ли шум
        
        Returns:
            Измеренный угол тангажа (рад)
        """
        theta_measured = theta_true + self.biases['attitude']
        
        if add_noise:
            noise = np.random.randn() * self.noise_levels['attitude']
            theta_measured += noise
        
        return theta_measured
    
    def _measure_load_factor(
        self,
        V: float,
        alpha: float,
        q: float,
        theta: float,
        delta_e: float,
        identified_params: np.ndarray,
        add_noise: bool = True
    ) -> float:
        """
        Модель акселерометра (нормальная перегрузка).
        
        nz = (L*cos(α) + D*sin(α)) / (m*g)
        
        Args:
            V, alpha, q, theta: Параметры состояния
            delta_e: Управление
            identified_params: Аэродинамические параметры
            add_noise: Добавлять ли шум
        
        Returns:
            Измеренная перегрузка (безразмерная, в единицах g)
        """
        from ..models.aerodynamics import AerodynamicCoefficients
        from ..models.aircraft_dynamics import AircraftDynamics
        
        # Вычисляем истинную перегрузку
        aero = AerodynamicCoefficients(self.params)
        alpha_dot = 0.0
        L, D, _ = aero.compute_forces_moments(V, alpha, alpha_dot, q, delta_e, identified_params)
        
        nz_true = (L * np.cos(alpha) + D * np.sin(alpha)) / (self.params.MASS * self.params.G)
        
        nz_measured = nz_true + self.biases['load_factor']
        
        if add_noise:
            noise = np.random.randn() * self.noise_levels['load_factor']
            nz_measured += noise
        
        return nz_measured
    
    def set_noise_level(self, sensor: str, std: float):
        """
        Устанавливает уровень шума для датчика.
        
        Args:
            sensor: Название датчика ('velocity', 'aoa', 'angular_rate', 'attitude', 'load_factor')
            std: Стандартное отклонение шума
        """
        if sensor in self.noise_levels:
            self.noise_levels[sensor] = std
        else:
            raise ValueError(f"Unknown sensor: {sensor}")
    
    def set_bias(self, sensor: str, bias: float):
        """
        Устанавливает смещение для датчика.
        
        Args:
            sensor: Название датчика
            bias: Величина смещения
        """
        if sensor in self.biases:
            self.biases[sensor] = bias
        else:
            raise ValueError(f"Unknown sensor: {sensor}")
    
    def reset_biases(self):
        """Сбрасывает все смещения в ноль."""
        for key in self.biases:
            self.biases[key] = 0.0
    
    def get_covariance_matrix(self) -> np.ndarray:
        """
        Возвращает ковариационную матрицу шумов измерений.
        
        Returns:
            R - диагональная матрица 5x5
        """
        R = np.diag([
            self.noise_levels['velocity']**2,
            self.noise_levels['aoa']**2,
            self.noise_levels['angular_rate']**2,
            self.noise_levels['attitude']**2,
            self.noise_levels['load_factor']**2
        ])
        return R
    
    def print_sensor_specs(self):
        """Выводит характеристики датчиков."""
        print("=" * 80)
        print("ХАРАКТЕРИСТИКИ ДАТЧИКОВ")
        print("=" * 80)
        print(f"{'Датчик':<25} {'Шум (СКО)':<20} {'Смещение':<15}")
        print("-" * 80)
        
        specs = [
            ('Воздушная скорость', 'velocity', 'м/с', 1.0),
            ('Угол атаки', 'aoa', '°', np.rad2deg(1.0)),
            ('Угловая скорость', 'angular_rate', '°/с', np.rad2deg(1.0)),
            ('Угол тангажа', 'attitude', '°', np.rad2deg(1.0)),
            ('Перегрузка', 'load_factor', 'g', 1.0)
        ]
        
        for name, key, unit, scale in specs:
            noise = self.noise_levels[key] * scale
            bias = self.biases[key] * scale
            print(f"{name:<25} {noise:.6f} {unit:<10} {bias:.6f} {unit}")
        
        print("=" * 80)


class SensorFaultSimulator:
    """
    Симулятор отказов датчиков для тестирования робастности.
    """
    
    def __init__(self, sensor_model: SensorModels):
        """
        Args:
            sensor_model: Базовая модель датчиков
        """
        self.sensor_model = sensor_model
        self.faults = {
            'velocity': False,
            'aoa': False,
            'angular_rate': False,
            'attitude': False,
            'load_factor': False
        }
    
    def inject_fault(self, sensor: str, fault_type: str = 'freeze', **kwargs):
        """
        Вводит отказ датчика.
        
        Args:
            sensor: Название датчика
            fault_type: Тип отказа ('freeze', 'bias', 'noise_increase', 'dropout')
            **kwargs: Параметры отказа
        """
        if sensor not in self.faults:
            raise ValueError(f"Unknown sensor: {sensor}")
        
        self.faults[sensor] = {
            'type': fault_type,
            'params': kwargs
        }
    
    def clear_fault(self, sensor: str):
        """Убирает отказ датчика."""
        self.faults[sensor] = False
    
    def clear_all_faults(self):
        """Убирает все отказы."""
        for sensor in self.faults:
            self.faults[sensor] = False
