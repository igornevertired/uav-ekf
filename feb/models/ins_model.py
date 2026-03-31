"""
Модель бесплатформенной инерциальной навигационной системы (БИНС).

Этот модуль реализует механизацию БИНС с ошибками гироскопов и акселерометров.
Используется для моделирования навигационных измерений и исследования влияния
точности навигации на идентификацию параметров.

Референции:
- Titterton, D., & Weston, J. (2004). Strapdown inertial navigation technology.
- Savage, P. G. (1998). Strapdown analytics.
"""

import numpy as np
from ..config.navigation_params import INSParameters, NavigationParams
from ..config.aircraft_params import AircraftParams


class StrapdownINS:
    """
    Модель бесплатформенной инерциальной навигационной системы.
    
    Состояние включает:
    - Позицию (x, y, z) в локальной системе координат
    - Скорость (Vx, Vy, Vz)
    - Ориентацию (кватернион q = [q0, q1, q2, q3])
    - Смещения гироскопов (3 компоненты)
    - Смещения акселерометров (3 компоненты)
    """
    
    def __init__(
        self,
        ins_params: INSParameters = None,
        initial_position: np.ndarray = None,
        initial_velocity: np.ndarray = None,
        initial_attitude: np.ndarray = None
    ):
        """
        Инициализация БИНС.
        
        Args:
            ins_params: Параметры точности БИНС
            initial_position: Начальная позиция [x, y, z] (м)
            initial_velocity: Начальная скорость [Vx, Vy, Vz] (м/с)
            initial_attitude: Начальная ориентация [roll, pitch, yaw] (рад)
        """
        if ins_params is None:
            ins_params, _ = NavigationParams.get_params('commercial')
        
        self.params = ins_params
        
        # Инициализация состояния
        self.position = initial_position if initial_position is not None else np.zeros(3)
        self.velocity = initial_velocity if initial_velocity is not None else np.array([150.0, 0.0, 0.0])
        
        # Кватернион ориентации (от тела к локальной СК)
        if initial_attitude is not None:
            self.quaternion = self._euler_to_quaternion(initial_attitude)
        else:
            self.quaternion = np.array([1.0, 0.0, 0.0, 0.0])  # Единичный кватернион
        
        # Смещения датчиков (инициализируем случайно)
        self.gyro_bias = np.random.randn(3) * self.params.get_gyro_bias_std()
        self.accel_bias = np.random.randn(3) * self.params.get_accel_bias_std()
        
        # Для случайного блуждания (random walk)
        self.gyro_rw_state = np.zeros(3)
        self.accel_rw_state = np.zeros(3)
        
        # История для отладки
        self.time = 0.0
    
    def propagate(
        self,
        dt: float,
        specific_force_body: np.ndarray,
        angular_rate_body: np.ndarray,
        include_errors: bool = True
    ) -> dict:
        """
        Распространяет состояние БИНС на один шаг.
        
        Args:
            dt: Временной шаг (сек)
            specific_force_body: Удельная сила в связанной СК [fx, fy, fz] (м/с²)
            angular_rate_body: Угловая скорость в связанной СК [wx, wy, wz] (рад/с)
            include_errors: Включать ли ошибки датчиков
        
        Returns:
            dict: Состояние БИНС {'position', 'velocity', 'attitude_euler', 'quaternion'}
        """
        # Добавляем ошибки датчиков
        if include_errors:
            # Измерения гироскопов с ошибками
            gyro_noise = np.random.randn(3) * self.params.get_gyro_noise_std(dt)
            gyro_measured = angular_rate_body + self.gyro_bias + gyro_noise
            
            # Случайное блуждание смещения гироскопа
            bias_rw = np.random.randn(3) * self.params.get_gyro_bias_std() * np.sqrt(dt) * 0.1
            self.gyro_bias += bias_rw
            
            # Измерения акселерометров с ошибками
            accel_noise = np.random.randn(3) * self.params.get_accel_noise_std(dt)
            accel_measured = specific_force_body + self.accel_bias + accel_noise
            
            # Случайное блуждание смещения акселерометра
            accel_bias_rw = np.random.randn(3) * self.params.get_accel_bias_std() * np.sqrt(dt) * 0.1
            self.accel_bias += accel_bias_rw
        else:
            gyro_measured = angular_rate_body
            accel_measured = specific_force_body
        
        # 1. Обновление ориентации (интегрирование кватерниона)
        self._update_attitude(gyro_measured, dt)
        
        # Сохраняем измеренную угловую скорость для использования в идентификации
        angular_rate_body = gyro_measured.copy()
        
        # 2. Преобразование удельной силы в локальную СК
        rotation_matrix = self._quaternion_to_dcm(self.quaternion)
        specific_force_local = rotation_matrix @ accel_measured
        
        # 3. Добавляем гравитацию (в локальной СК вниз - отрицательное z)
        g_local = np.array([0.0, 0.0, -AircraftParams.G])
        acceleration_local = specific_force_local + g_local
        
        # 4. Обновление скорости
        self.velocity += acceleration_local * dt
        
        # 5. Обновление позиции
        self.position += self.velocity * dt
        
        self.time += dt
        
        # Возвращаем текущее состояние
        return {
            'position': self.position.copy(),
            'velocity': self.velocity.copy(),
            'attitude_euler': self._quaternion_to_euler(self.quaternion),
            'quaternion': self.quaternion.copy(),
            'angular_rate_body': angular_rate_body,
            'time': self.time
        }
    
    def _update_attitude(self, omega_body: np.ndarray, dt: float):
        """
        Обновляет кватернион ориентации.
        
        Args:
            omega_body: Угловая скорость в связанной СК (рад/с)
            dt: Временной шаг (сек)
        """
        # Матрица кватернионной производной
        # q̇ = 0.5 * Ω(ω) * q
        omega_mag = np.linalg.norm(omega_body)
        
        if omega_mag < 1e-8:
            # Нет вращения
            return
        
        # Кватернион приращения вращения
        half_angle = 0.5 * omega_mag * dt
        sin_half = np.sin(half_angle)
        cos_half = np.cos(half_angle)
        
        axis = omega_body / omega_mag
        delta_q = np.array([
            cos_half,
            sin_half * axis[0],
            sin_half * axis[1],
            sin_half * axis[2]
        ])
        
        # Умножение кватернионов: q_new = delta_q ⊗ q_old
        self.quaternion = self._quaternion_multiply(delta_q, self.quaternion)
        
        # Нормализация
        self.quaternion /= np.linalg.norm(self.quaternion)
    
    def _quaternion_multiply(self, q1: np.ndarray, q2: np.ndarray) -> np.ndarray:
        """
        Умножение кватернионов: q1 ⊗ q2.
        
        Args:
            q1, q2: Кватернионы [q0, q1, q2, q3]
        
        Returns:
            Результат умножения
        """
        w1, x1, y1, z1 = q1
        w2, x2, y2, z2 = q2
        
        return np.array([
            w1*w2 - x1*x2 - y1*y2 - z1*z2,
            w1*x2 + x1*w2 + y1*z2 - z1*y2,
            w1*y2 - x1*z2 + y1*w2 + z1*x2,
            w1*z2 + x1*y2 - y1*x2 + z1*w2
        ])
    
    def _quaternion_to_dcm(self, q: np.ndarray) -> np.ndarray:
        """
        Преобразует кватернион в матрицу направляющих косинусов (DCM).
        
        Args:
            q: Кватернион [q0, q1, q2, q3]
        
        Returns:
            Матрица вращения 3x3 (от тела к локальной СК)
        """
        q0, q1, q2, q3 = q
        
        dcm = np.array([
            [q0**2 + q1**2 - q2**2 - q3**2, 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2)],
            [2*(q1*q2 + q0*q3), q0**2 - q1**2 + q2**2 - q3**2, 2*(q2*q3 - q0*q1)],
            [2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), q0**2 - q1**2 - q2**2 + q3**2]
        ])
        
        return dcm
    
    def _euler_to_quaternion(self, euler: np.ndarray) -> np.ndarray:
        """
        Преобразует углы Эйлера в кватернион.
        
        Args:
            euler: Углы Эйлера [roll, pitch, yaw] (рад)
        
        Returns:
            Кватернион [q0, q1, q2, q3]
        """
        roll, pitch, yaw = euler
        
        cy = np.cos(yaw * 0.5)
        sy = np.sin(yaw * 0.5)
        cp = np.cos(pitch * 0.5)
        sp = np.sin(pitch * 0.5)
        cr = np.cos(roll * 0.5)
        sr = np.sin(roll * 0.5)
        
        q0 = cr * cp * cy + sr * sp * sy
        q1 = sr * cp * cy - cr * sp * sy
        q2 = cr * sp * cy + sr * cp * sy
        q3 = cr * cp * sy - sr * sp * cy
        
        return np.array([q0, q1, q2, q3])
    
    def _quaternion_to_euler(self, q: np.ndarray) -> np.ndarray:
        """
        Преобразует кватернион в углы Эйлера.
        
        Args:
            q: Кватернион [q0, q1, q2, q3]
        
        Returns:
            Углы Эйлера [roll, pitch, yaw] (рад)
        """
        q0, q1, q2, q3 = q
        
        # Roll (φ)
        roll = np.arctan2(2*(q0*q1 + q2*q3), 1 - 2*(q1**2 + q2**2))
        
        # Pitch (θ)
        sin_pitch = 2*(q0*q2 - q3*q1)
        sin_pitch = np.clip(sin_pitch, -1.0, 1.0)
        pitch = np.arcsin(sin_pitch)
        
        # Yaw (ψ)
        yaw = np.arctan2(2*(q0*q3 + q1*q2), 1 - 2*(q2**2 + q3**2))
        
        return np.array([roll, pitch, yaw])
    
    def get_state(self) -> dict:
        """
        Возвращает текущее состояние БИНС.
        
        Returns:
            dict: Полное состояние системы
        """
        return {
            'position': self.position.copy(),
            'velocity': self.velocity.copy(),
            'attitude_euler': self._quaternion_to_euler(self.quaternion),
            'quaternion': self.quaternion.copy(),
            'gyro_bias': self.gyro_bias.copy(),
            'accel_bias': self.accel_bias.copy(),
            'time': self.time
        }
    
    def reset(
        self,
        position: np.ndarray = None,
        velocity: np.ndarray = None,
        attitude: np.ndarray = None
    ):
        """
        Сбрасывает состояние БИНС.
        
        Args:
            position: Новая позиция (если None, остается прежней)
            velocity: Новая скорость (если None, остается прежней)
            attitude: Новая ориентация в углах Эйлера (если None, остается прежней)
        """
        if position is not None:
            self.position = position.copy()
        if velocity is not None:
            self.velocity = velocity.copy()
        if attitude is not None:
            self.quaternion = self._euler_to_quaternion(attitude)
        
        self.time = 0.0
