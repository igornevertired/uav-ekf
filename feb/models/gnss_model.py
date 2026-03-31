"""
Модель приемника глобальной навигационной спутниковой системы (ГНСС).

Этот модуль моделирует выходные данные ГНСС приемника с характерными ошибками:
- Белый шум измерений позиции и скорости
- Возможность потери сигнала
- Различные классы точности приемников
"""

import numpy as np
from ..config.navigation_params import GNSSParameters, NavigationParams


class GNSSReceiver:
    """
    Модель приемника ГНСС (GPS/ГЛОНАСС/Galileo и т.д.).
    
    Выдает измерения позиции и скорости с шумом согласно классу точности.
    """
    
    def __init__(
        self,
        gnss_params: GNSSParameters = None,
        update_rate: float = None
    ):
        """
        Инициализация модели ГНСС приемника.
        
        Args:
            gnss_params: Параметры точности ГНСС
            update_rate: Частота обновления (Гц). Если None, берется из gnss_params
        """
        if gnss_params is None:
            _, gnss_params = NavigationParams.get_params('commercial')
        
        self.params = gnss_params
        self.update_rate = update_rate if update_rate is not None else gnss_params.update_rate
        self.update_period = 1.0 / self.update_rate
        
        # Внутреннее время для контроля обновлений
        self.time_since_update = 0.0
        self.last_measurement = None
        
        # Статистика для анализа
        self.measurement_count = 0
        self.total_time = 0.0
    
    def get_measurement(
        self,
        true_position: np.ndarray,
        true_velocity: np.ndarray,
        dt: float,
        available: bool = True
    ) -> dict | None:
        """
        Получает зашумленные измерения позиции и скорости.
        
        Args:
            true_position: Истинная позиция [x, y, z] (м)
            true_velocity: Истинная скорость [Vx, Vy, Vz] (м/с)
            dt: Временной шаг с предыдущего вызова (сек)
            available: Доступен ли сигнал ГНСС
        
        Returns:
            dict с ключами {'position', 'velocity', 'std_position', 'std_velocity', 'time'}
            или None если измерение недоступно
        """
        self.time_since_update += dt
        self.total_time += dt
        
        # Проверяем, пора ли выдавать новое измерение
        if self.time_since_update < self.update_period:
            # Возвращаем предыдущее измерение или None
            return self.last_measurement
        
        # Сбрасываем таймер (учитываем возможное превышение)
        self.time_since_update = self.time_since_update % self.update_period
        
        # Проверяем доступность сигнала
        if not available:
            self.last_measurement = None
            return None
        
        # Генерируем шумы
        position_noise = np.random.randn(3) * self.params.position_std
        velocity_noise = np.random.randn(3) * self.params.velocity_std

        # Формируем измерение
        measurement = {
            'position': true_position + position_noise,
            'velocity': true_velocity + velocity_noise,
            'std_position': self.params.position_std,
            'std_velocity': self.params.velocity_std,
            'time': self.total_time,
            'available': True
        }
        
        self.last_measurement = measurement
        self.measurement_count += 1
        
        return measurement
    
    def is_measurement_ready(self, dt: float) -> bool:
        """
        Проверяет, готово ли новое измерение.
        
        Args:
            dt: Прошедшее время с последнего вызова (сек)
        
        Returns:
            True если пора выдавать новое измерение
        """
        self.time_since_update += dt
        return self.time_since_update >= self.update_period
    
    def reset(self):
        """Сбрасывает внутреннее состояние приемника."""
        self.time_since_update = 0.0
        self.last_measurement = None
        self.measurement_count = 0
        self.total_time = 0.0
    
    def get_statistics(self) -> dict:
        """
        Возвращает статистику работы приемника.
        
        Returns:
            dict: Статистическая информация
        """
        avg_rate = self.measurement_count / self.total_time if self.total_time > 0 else 0
        
        return {
            'measurement_count': self.measurement_count,
            'total_time': self.total_time,
            'average_rate': avg_rate,
            'nominal_rate': self.update_rate,
            'position_std': self.params.position_std,
            'velocity_std': self.params.velocity_std
        }
    
    def simulate_signal_loss(self, probability: float = 0.0) -> bool:
        """
        Моделирует потерю сигнала ГНСС.
        
        Args:
            probability: Вероятность потери сигнала на данном измерении [0, 1]
        
        Returns:
            True если сигнал доступен, False если потерян
        """
        return np.random.rand() > probability
    
    def get_covariance_matrix(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Возвращает матрицы ковариации ошибок измерений.
        
        Returns:
            tuple: (R_position, R_velocity) - ковариационные матрицы 3x3
        """
        R_position = np.eye(3) * self.params.position_std ** 2
        R_velocity = np.eye(3) * self.params.velocity_std ** 2
        
        return R_position, R_velocity
    
    def print_info(self):
        """Выводит информацию о параметрах приемника."""
        print("=" * 80)
        print("ПАРАМЕТРЫ ПРИЕМНИКА ГНСС")
        print("=" * 80)
        print(f"Точность позиции (σ): {self.params.position_std:.2f} м")
        print(f"Точность скорости (σ): {self.params.velocity_std:.3f} м/с")
        print(f"Частота обновления: {self.update_rate:.1f} Гц")
        print(f"Период обновления: {self.update_period:.3f} сек")
        
        stats = self.get_statistics()
        if stats['total_time'] > 0:
            print("-" * 80)
            print("СТАТИСТИКА:")
            print(f"  Всего измерений: {stats['measurement_count']}")
            print(f"  Время работы: {stats['total_time']:.1f} сек")
            print(f"  Средняя частота: {stats['average_rate']:.2f} Гц")
        print("=" * 80)


class GNSSSimulator:
    """
    Расширенный симулятор ГНСС с дополнительными эффектами.
    
    Может моделировать:
    - Multipath эффекты
    - Ionospheric delays
    - Satellite geometry (DOP)
    """
    
    def __init__(self, base_receiver: GNSSReceiver):
        """
        Args:
            base_receiver: Базовый приемник ГНСС
        """
        self.receiver = base_receiver
        self.multipath_enabled = False
        self.ionospheric_error_enabled = False
    
    def enable_multipath(self, amplitude: float = 2.0):
        """
        Включает моделирование multipath эффектов.
        
        Args:
            amplitude: Амплитуда multipath ошибки (м)
        """
        self.multipath_enabled = True
        self.multipath_amplitude = amplitude
    
    def enable_ionospheric_error(self, std: float = 5.0):
        """
        Включает моделирование ионосферных задержек.
        
        Args:
            std: Стандартное отклонение ошибки (м)
        """
        self.ionospheric_error_enabled = True
        self.ionospheric_std = std
    
    def get_measurement(
        self,
        true_position: np.ndarray,
        true_velocity: np.ndarray,
        dt: float,
        available: bool = True
    ) -> dict | None:
        """
        Получает измерение с дополнительными эффектами.
        
        Args:
            true_position: Истинная позиция
            true_velocity: Истинная скорость
            dt: Временной шаг
            available: Доступность сигнала
        
        Returns:
            Измерение или None
        """
        # Базовое измерение
        measurement = self.receiver.get_measurement(
            true_position, true_velocity, dt, available
        )
        
        if measurement is None or not measurement['available']:
            return measurement
        
        # Добавляем multipath
        if self.multipath_enabled:
            # Медленно меняющаяся ошибка
            multipath_error = np.random.randn(3) * self.multipath_amplitude * 0.5
            measurement['position'] += multipath_error
        
        # Добавляем ионосферные задержки
        if self.ionospheric_error_enabled:
            # Коррелированная ошибка (влияет в основном на вертикаль)
            iono_error = np.array([0, 0, np.random.randn() * self.ionospheric_std])
            measurement['position'] += iono_error
        
        return measurement
