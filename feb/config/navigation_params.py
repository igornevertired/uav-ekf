"""
Параметры навигационных систем (БИНС и ГНСС).

Этот модуль определяет характеристики точности навигационных датчиков
для различных классов оборудования:
- Commercial (коммерческое)
- Tactical (тактическое)
- Navigation (навигационное)
"""

import numpy as np
from dataclasses import dataclass
from typing import Literal


@dataclass
class INSParameters:
    """
    Параметры инерциальной навигационной системы (БИНС).
    
    Attributes:
        gyro_bias_stability: Стабильность смещения гироскопа (deg/hr)
        gyro_arw: Угловое случайное блуждание (deg/√hr)
        accel_bias_stability: Стабильность смещения акселерометра (mg)
        accel_vrw: Случайное блуждание скорости (m/s/√hr)
        update_rate: Частота обновления (Гц)
    """
    gyro_bias_stability: float  # deg/hr
    gyro_arw: float  # deg/√hr
    accel_bias_stability: float  # mg (миллиграммы)
    accel_vrw: float  # m/s/√hr
    update_rate: float  # Hz
    
    def get_gyro_noise_std(self, dt: float) -> float:
        """
        Вычисляет СКО шума гироскопа для заданного временного шага.
        
        Args:
            dt: Временной шаг (секунды)
            
        Returns:
            СКО в рад/с
        """
        # ARW в рад/√с
        arw_rad_per_sqrt_s = np.deg2rad(self.gyro_arw) / np.sqrt(3600)
        return arw_rad_per_sqrt_s / np.sqrt(dt)
    
    def get_accel_noise_std(self, dt: float) -> float:
        """
        Вычисляет СКО шума акселерометра для заданного временного шага.
        
        Args:
            dt: Временной шаг (секунды)
            
        Returns:
            СКО в м/с²
        """
        # VRW в м/с²/√с
        vrw_per_sqrt_s = self.accel_vrw / np.sqrt(3600)
        return vrw_per_sqrt_s / np.sqrt(dt)
    
    def get_gyro_bias_std(self) -> float:
        """
        Возвращает СКО смещения гироскопа.
        
        Returns:
            СКО в рад/с
        """
        # Конвертация из deg/hr в рад/с
        return np.deg2rad(self.gyro_bias_stability) / 3600
    
    def get_accel_bias_std(self) -> float:
        """
        Возвращает СКО смещения акселерометра.
        
        Returns:
            СКО в м/с²
        """
        # Конвертация из mg в м/с²
        g = 9.81
        return self.accel_bias_stability * g / 1000


@dataclass
class GNSSParameters:
    """
    Параметры приемника ГНСС.
    
    Attributes:
        position_std: СКО ошибки позиции (м)
        velocity_std: СКО ошибки скорости (м/с)
        update_rate: Частота обновления (Гц)
    """
    position_std: float  # м
    velocity_std: float  # м/с
    update_rate: float  # Hz


class NavigationParams:
    """
    Класс, содержащий параметры навигационных систем для разных классов точности.
    
    Single Responsibility: только конфигурация навигационного оборудования.
    """
    
    # === КОММЕРЧЕСКИЙ КЛАСС (COMMERCIAL) ===
    # Низкая стоимость, умеренная точность
    INS_COMMERCIAL = INSParameters(
        gyro_bias_stability=10.0,  # deg/hr
        gyro_arw=0.5,  # deg/√hr
        accel_bias_stability=1.0,  # mg
        accel_vrw=0.05,  # m/s/√hr
        update_rate=100.0  # Hz
    )
    
    GNSS_COMMERCIAL = GNSSParameters(
        position_std=7.5,  # м (среднее между 5-10м)
        velocity_std=0.15,  # м/с
        update_rate=5.0  # Hz
    )
    
    # === ТАКТИЧЕСКИЙ КЛАСС (TACTICAL) ===
    # Средняя стоимость, хорошая точность
    INS_TACTICAL = INSParameters(
        gyro_bias_stability=1.0,  # deg/hr
        gyro_arw=0.1,  # deg/√hr
        accel_bias_stability=0.1,  # mg
        accel_vrw=0.01,  # m/s/√hr
        update_rate=200.0  # Hz
    )
    
    GNSS_TACTICAL = GNSSParameters(
        position_std=2.0,  # м
        velocity_std=0.05,  # м/с
        update_rate=10.0  # Hz
    )
    
    # === НАВИГАЦИОННЫЙ КЛАСС (NAVIGATION) ===
    # Высокая стоимость, высокая точность
    INS_NAVIGATION = INSParameters(
        gyro_bias_stability=0.01,  # deg/hr
        gyro_arw=0.001,  # deg/√hr
        accel_bias_stability=0.01,  # mg
        accel_vrw=0.001,  # m/s/√hr
        update_rate=400.0  # Hz
    )
    
    GNSS_NAVIGATION = GNSSParameters(
        position_std=0.5,  # м
        velocity_std=0.01,  # м/с
        update_rate=20.0  # Hz
    )
    
    @classmethod
    def get_params(cls, accuracy_class: Literal['commercial', 'tactical', 'navigation'] = 'commercial'):
        """
        Возвращает параметры БИНС и ГНСС для указанного класса точности.
        
        Args:
            accuracy_class: Класс точности ('commercial', 'tactical', 'navigation')
            
        Returns:
            tuple: (INSParameters, GNSSParameters)
        """
        if accuracy_class == 'commercial':
            return cls.INS_COMMERCIAL, cls.GNSS_COMMERCIAL
        elif accuracy_class == 'tactical':
            return cls.INS_TACTICAL, cls.GNSS_TACTICAL
        elif accuracy_class == 'navigation':
            return cls.INS_NAVIGATION, cls.GNSS_NAVIGATION
        else:
            raise ValueError(f"Unknown accuracy class: {accuracy_class}")
    
    @classmethod
    def print_comparison(cls):
        """Выводит таблицу сравнения классов точности."""
        print("=" * 100)
        print("СРАВНЕНИЕ КЛАССОВ ТОЧНОСТИ НАВИГАЦИОННОГО ОБОРУДОВАНИЯ")
        print("=" * 100)
        print(f"{'Параметр':<40} {'Commercial':<20} {'Tactical':<20} {'Navigation':<20}")
        print("-" * 100)
        
        # БИНС
        print("БИНС:")
        print(f"  {'Смещение гироскопа (deg/hr)':<38} {cls.INS_COMMERCIAL.gyro_bias_stability:<20} "
              f"{cls.INS_TACTICAL.gyro_bias_stability:<20} {cls.INS_NAVIGATION.gyro_bias_stability:<20}")
        print(f"  {'ARW гироскопа (deg/√hr)':<38} {cls.INS_COMMERCIAL.gyro_arw:<20} "
              f"{cls.INS_TACTICAL.gyro_arw:<20} {cls.INS_NAVIGATION.gyro_arw:<20}")
        print(f"  {'Смещение акселерометра (mg)':<38} {cls.INS_COMMERCIAL.accel_bias_stability:<20} "
              f"{cls.INS_TACTICAL.accel_bias_stability:<20} {cls.INS_NAVIGATION.accel_bias_stability:<20}")
        print(f"  {'VRW акселерометра (m/s/√hr)':<38} {cls.INS_COMMERCIAL.accel_vrw:<20} "
              f"{cls.INS_TACTICAL.accel_vrw:<20} {cls.INS_NAVIGATION.accel_vrw:<20}")
        print(f"  {'Частота обновления (Hz)':<38} {cls.INS_COMMERCIAL.update_rate:<20} "
              f"{cls.INS_TACTICAL.update_rate:<20} {cls.INS_NAVIGATION.update_rate:<20}")
        
        # ГНСС
        print("\nГНСС:")
        print(f"  {'Точность позиции (м)':<38} {cls.GNSS_COMMERCIAL.position_std:<20} "
              f"{cls.GNSS_TACTICAL.position_std:<20} {cls.GNSS_NAVIGATION.position_std:<20}")
        print(f"  {'Точность скорости (м/с)':<38} {cls.GNSS_COMMERCIAL.velocity_std:<20} "
              f"{cls.GNSS_TACTICAL.velocity_std:<20} {cls.GNSS_NAVIGATION.velocity_std:<20}")
        print(f"  {'Частота обновления (Hz)':<38} {cls.GNSS_COMMERCIAL.update_rate:<20} "
              f"{cls.GNSS_TACTICAL.update_rate:<20} {cls.GNSS_NAVIGATION.update_rate:<20}")
        
        print("=" * 100)
