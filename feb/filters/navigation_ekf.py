"""
Расширенный фильтр Калмана для комплексирования БИНС и ГНСС.

Этот модуль реализует Navigation EKF для объединения данных инерциальной
навигационной системы (БИНС) и глобальной навигационной спутниковой системы (ГНСС).

Архитектура: Loose Coupling (слабосвязанное комплексирование)
- БИНС распространяется автономно
- ГНСС предоставляет измерения позиции и скорости
- EKF корректирует ошибки БИНС по измерениям ГНСС

Состояние ошибок (error-state):
[δposition(3), δvelocity(3), attitude_errors(3), gyro_bias(3), accel_bias(3)]
"""

import numpy as np
from .ekf_base import EKFBase
from ..models.ins_model import StrapdownINS
from ..models.gnss_model import GNSSReceiver
from ..config.navigation_params import INSParameters, GNSSParameters
from ..config.aircraft_params import AircraftParams


class NavigationIntegrationEKF(EKFBase):
    """
    EKF для комплексирования БИНС и ГНСС методом слабосвязанной интеграции.
    
    Оценивает ошибки состояния БИНС и корректирует их по измерениям ГНСС.
    """
    
    def __init__(
        self,
        dt: float = 0.01,
        ins_model: StrapdownINS = None,
        gnss_model: GNSSReceiver = None,
        ins_params: INSParameters = None,
        gnss_params: GNSSParameters = None
    ):
        """
        Инициализация навигационного EKF.
        
        Args:
            dt: Временной шаг БИНС (сек)
            ins_model: Модель БИНС
            gnss_model: Модель ГНСС приемника
            ins_params: Параметры точности БИНС
            gnss_params: Параметры точности ГНСС
        """
        self.ins = ins_model
        self.gnss = gnss_model
        self.ins_params = ins_params
        self.gnss_params = gnss_params
        
        # Размерность состояния ошибок: 15
        # [δpos(3), δvel(3), δatt(3), gyro_bias(3), accel_bias(3)]
        state_dim = 15
        
        # Размерность измерений: 6 (позиция + скорость из ГНСС)
        measurement_dim = 6
        
        # Начальное состояние ошибок (все нули)
        initial_state = np.zeros(state_dim)
        
        # Начальная ковариация (увеличена для лучшей сходимости)
        initial_P = np.diag([
            1000.0, 1000.0, 1000.0,  # Ошибки позиции (м²) - увеличено 10x
            100.0, 100.0, 100.0,     # Ошибки скорости (м²/с²) - увеличено 10x
            1.0, 1.0, 1.0,           # Ошибки ориентации (рад²) - увеличено 10x
            (ins_params.get_gyro_bias_std())**2 * 10, (ins_params.get_gyro_bias_std())**2 * 10, (ins_params.get_gyro_bias_std())**2 * 10,  # Смещения гироскопов
            (ins_params.get_accel_bias_std())**2 * 10, (ins_params.get_accel_bias_std())**2 * 10, (ins_params.get_accel_bias_std())**2 * 10  # Смещения акселерометров
        ]) if ins_params is not None else np.eye(state_dim)
        
        # Ковариация шума процесса (увеличена для учета дрейфа БИНС)
        if ins_params is not None:
            gyro_noise_var = (ins_params.get_gyro_noise_std(dt))**2
            accel_noise_var = (ins_params.get_accel_noise_std(dt))**2
            gyro_bias_rw_var = (ins_params.get_gyro_bias_std() * np.sqrt(dt) * 0.1)**2
            accel_bias_rw_var = (ins_params.get_accel_bias_std() * np.sqrt(dt) * 0.1)**2
            
            # Увеличиваем Q для учета быстрого роста ошибок БИНС
            Q = np.diag([
                accel_noise_var * dt**2 * 100,  # Ошибки позиции - увеличено 100x для квадратичного роста
                accel_noise_var * dt**2 * 100, 
                accel_noise_var * dt**2 * 100,
                accel_noise_var * 50,           # Ошибки скорости - увеличено 50x
                accel_noise_var * 50,
                accel_noise_var * 50,
                gyro_noise_var * 10,            # Ошибки ориентации - увеличено 10x
                gyro_noise_var * 10,
                gyro_noise_var * 10,
                gyro_bias_rw_var * 5,           # Блуждание bias гироскопов - увеличено 5x
                gyro_bias_rw_var * 5,
                gyro_bias_rw_var * 5,
                accel_bias_rw_var * 5,          # Блуждание bias акселерометров - увеличено 5x
                accel_bias_rw_var * 5,
                accel_bias_rw_var * 5
            ])
        else:
            Q = np.eye(state_dim) * 1e-6
        
        # Ковариация шума измерений ГНСС
        if gnss_params is not None:
            R = np.diag([
                gnss_params.position_std**2, gnss_params.position_std**2, gnss_params.position_std**2,
                gnss_params.velocity_std**2, gnss_params.velocity_std**2, gnss_params.velocity_std**2
            ])
        else:
            R = np.eye(measurement_dim) * 1e-6
        
        # Инициализация базового класса
        super().__init__(
            state_dim=state_dim,
            measurement_dim=measurement_dim,
            dt=dt,
            initial_state=initial_state,
            initial_P=initial_P,
            Q=Q,
            R=R
        )
        
        # Дополнительная история для навигации
        self.position_history = []
        self.velocity_history = []
        self.attitude_history = []
        
        print("NavigationEKF инициализирован (Loose Coupling)")
    
    def f_dynamics(self, x: np.ndarray, u: any, dt: float) -> np.ndarray:
        """
        Динамика ошибок состояния БИНС.
        
        Для error-state фильтра ошибки распространяются линейно:
        δx_{k+1} = F * δx_k + w_k
        
        Args:
            x: Вектор ошибок состояния [δpos, δvel, δatt, gyro_bias, accel_bias]
            u: Не используется (динамика автономна)
            dt: Временной шаг
        
        Returns:
            Следующий вектор ошибок
        """
        # В error-state EKF ошибки распространяются через якобиан
        # Здесь возвращаем текущее состояние (интеграция будет в predict)
        return x
    
    def h_measurement(self, x: np.ndarray, u: any) -> np.ndarray:
        """
        Функция измерений для ошибок.
        
        Измерения ГНСС: позиция и скорость
        Инновация: y = (GNSS_measurement) - (INS_estimate)
        
        В error-state: y = H * δx
        где H = [I_{6x6}, 0_{6x9}] - измеряем только ошибки позиции и скорости
        
        Args:
            x: Вектор ошибок состояния
            u: Текущие оценки БИНС (для сравнения с ГНСС)
        
        Returns:
            Прогнозируемые ошибки измерений (обычно нули для error-state)
        """
        # В error-state EKF прогнозируемые измерения ошибок - это сами ошибки
        return x[:6]  # Ошибки позиции и скорости
    
    def jacobian_F(self, x: np.ndarray, u: any) -> np.ndarray:
        """
        Якобиан динамики ошибок.
        
        Матрица F описывает, как ошибки распространяются во времени.
        
        Returns:
            F (15x15)
        """
        F = np.eye(15)
        
        # Простая модель: позиция интегрирует скорость
        # δpos_{k+1} = δpos_k + δvel_k * dt
        F[0:3, 3:6] = np.eye(3) * self.dt
        
        # Ориентация интегрирует угловую скорость (через bias гироскопов)
        # δatt_{k+1} = δatt_k + gyro_bias * dt
        F[6:9, 9:12] = -np.eye(3) * self.dt
        
        # Смещения эволюционируют как случайное блуждание (остаются близкими)
        # bias_{k+1} = bias_k + noise
        # Уже есть в единичной матрице
        
        return F
    
    def jacobian_H(self, x: np.ndarray, u: any) -> np.ndarray:
        """
        Якобиан измерений для ошибок.
        
        ГНСС измеряет позицию и скорость, поэтому H = [I_{6x6}, 0_{6x9}]
        
        Returns:
            H (6x15)
        """
        H = np.zeros((6, 15))
        H[0:6, 0:6] = np.eye(6)  # Измеряем ошибки позиции и скорости
        return H
    
    def update_with_gnss(
        self,
        gnss_measurement: dict,
        ins_state: dict
    ):
        """
        Обновление фильтра по измерениям ГНСС.
        
        Args:
            gnss_measurement: Словарь с ключами 'position', 'velocity'
            ins_state: Текущее состояние БИНС с ключами 'position', 'velocity'
        """
        # Инновация: разница между ГНСС и БИНС
        # Это оценка текущих ошибок БИНС
        innovation = np.hstack([
            gnss_measurement['position'] - ins_state['position'],
            gnss_measurement['velocity'] - ins_state['velocity']
        ])
        
        # Обновление error-state через стандартный EKF update
        # innovation уже содержит оценку ошибок, поэтому передаем его напрямую
        y_measured = innovation  # Реальная инновация
        y_predicted = self.h_measurement(self.x, u=None)  # Предсказанные ошибки
        
        # Используем базовый update с правильной инновацией
        actual_innovation = y_measured - y_predicted
        
        H = self.jacobian_H(self.x, u=None)
        S = H @ self.P @ H.T + self.R
        S = (S + S.T) / 2
        
        try:
            K = self.P @ H.T @ np.linalg.inv(S)
            self.x = self.x + K @ actual_innovation
            I_KH = np.eye(self.state_dim) - K @ H
            self.P = I_KH @ self.P @ I_KH.T + K @ self.R @ K.T
        except np.linalg.LinAlgError:
            K = self.P @ H.T @ np.linalg.pinv(S)
            self.x = self.x + K @ actual_innovation
            I_KH = np.eye(self.state_dim) - K @ H
            self.P = I_KH @ self.P @ I_KH.T + K @ self.R @ K.T
        
        self.P = (self.P + self.P.T) / 2
        self._ensure_positive_definite(self.P)
        
        # Сохранение в историю
        self._save_to_history(actual_innovation)
        
        # После обновления применяем коррекцию к БИНС (делается в симуляторе)

        # === Дополнительная подпитка ориентации (pitch) по направлению скорости GNSS ===
        # Loose coupling в базовом виде корректирует только pos/vel, поэтому pitch остаётся почти INS-only.
        # Для 2D/продольной модели добавим псевдо-измерение:
        #   theta_meas ≈ gamma + alpha_trim
        # где gamma — путевой угол по GNSS-скорости в NED, alpha_trim — типичный угол атаки на крейсере.
        try:
            v_gnss = gnss_measurement.get('velocity', None)
            if v_gnss is not None:
                vx, vy, vz = float(v_gnss[0]), float(v_gnss[1]), float(v_gnss[2])
                v_h = np.sqrt(vx * vx + vy * vy)
                if v_h > 1e-3:
                    gamma = np.arctan2(-vz, v_h)
                    theta_meas = gamma + float(AircraftParams.TRIM_ALPHA)
                    theta_ins = float(ins_state['attitude_euler'][1])

                    # Измеряем ошибку ориентации: delta_theta ≈ theta_meas - theta_ins
                    z = theta_meas - theta_ins

                    # Модель измерения: z ≈ δatt_pitch = x[7]
                    Ht = np.zeros((1, self.state_dim))
                    Ht[0, 7] = 1.0

                    # Оценка дисперсии gamma из дисперсии скорости GNSS: sigma_gamma ≈ sigma_v / V
                    sigma_v = float(gnss_measurement.get('std_velocity', self.gnss_params.velocity_std if self.gnss_params else 0.1))
                    V_eff = max(np.linalg.norm(v_gnss), 1.0)
                    sigma_theta = max(sigma_v / V_eff, np.deg2rad(0.5))  # не меньше 0.5°
                    Rt = np.array([[sigma_theta ** 2]])

                    # 1D EKF update (Joseph form)
                    y_pred = float(self.x[7])
                    innov = z - y_pred
                    S_t = Ht @ self.P @ Ht.T + Rt
                    Kt = self.P @ Ht.T @ np.linalg.inv(S_t)
                    self.x = self.x + (Kt[:, 0] * innov)
                    I_KH = np.eye(self.state_dim) - Kt @ Ht
                    self.P = I_KH @ self.P @ I_KH.T + Kt @ Rt @ Kt.T
                    self.P = (self.P + self.P.T) / 2
                    self._ensure_positive_definite(self.P)
        except Exception:
            # Псевдо-измерение не должно ломать навигацию
            pass
    
    def get_corrected_state(self, ins_state: dict) -> dict:
        """
        Применяет оценки ошибок к состоянию БИНС для получения скорректированного состояния.
        
        В error-state EKF:
        δx = истинное - БИНС
        => истинное = БИНС + δx
        => скорректированное = БИНС + δx (NOT БИНС - δx!)
        
        Args:
            ins_state: Текущее состояние БИНС
        
        Returns:
            Скорректированное состояние
        """
        δpos = self.x[0:3]
        δvel = self.x[3:6]
        δatt = self.x[6:9]
        
        corrected = {
            'position': ins_state['position'] + δpos,  # ← ИСПРАВЛЕНО: + вместо -
            'velocity': ins_state['velocity'] + δvel,  # ← ИСПРАВЛЕНО: + вместо -
            'attitude_euler': ins_state['attitude_euler'] + δatt,  # ← ИСПРАВЛЕНО: + вместо -
            'time': ins_state['time']
        }
        
        return corrected
    
    def reset_errors(self):
        """
        Сбрасывает оценки ошибок после коррекции БИНС.
        
        Это важно в error-state EKF: после применения коррекции к БИНС,
        ошибки сбрасываются в ноль.
        """
        self.x[:9] = 0.0  # Сбрасываем ошибки pos, vel, att
        # bias оставляем (они оцениваются накопительно)
    
    def get_navigation_solution(self) -> dict:
        """
        Возвращает текущее навигационное решение с оценками точности.
        
        Returns:
            dict: Навигационное решение с ковариациями
        """
        position_std = np.sqrt(np.diag(self.P)[0:3])
        velocity_std = np.sqrt(np.diag(self.P)[3:6])
        attitude_std = np.sqrt(np.diag(self.P)[6:9])
        
        return {
            'position_errors': self.x[0:3],
            'velocity_errors': self.x[3:6],
            'attitude_errors': self.x[6:9],
            'gyro_bias': self.x[9:12],
            'accel_bias': self.x[12:15],
            'position_std': position_std,
            'velocity_std': velocity_std,
            'attitude_std': attitude_std
        }
    
    def print_navigation_accuracy(self):
        """Выводит информацию о точности навигационного решения."""
        sol = self.get_navigation_solution()
        
        print("=" * 80)
        print("ТОЧНОСТЬ НАВИГАЦИОННОГО РЕШЕНИЯ")
        print("=" * 80)
        print("Позиция (СКО):")
        print(f"  X: {sol['position_std'][0]:.2f} м")
        print(f"  Y: {sol['position_std'][1]:.2f} м")
        print(f"  Z: {sol['position_std'][2]:.2f} м")
        print(f"  Горизонталь: {np.linalg.norm(sol['position_std'][:2]):.2f} м")
        print("\nСкорость (СКО):")
        print(f"  Vx: {sol['velocity_std'][0]:.3f} м/с")
        print(f"  Vy: {sol['velocity_std'][1]:.3f} м/с")
        print(f"  Vz: {sol['velocity_std'][2]:.3f} м/с")
        print("\nОриентация (СКО):")
        print(f"  Roll: {np.rad2deg(sol['attitude_std'][0]):.2f}°")
        print(f"  Pitch: {np.rad2deg(sol['attitude_std'][1]):.2f}°")
        print(f"  Yaw: {np.rad2deg(sol['attitude_std'][2]):.2f}°")
        print("=" * 80)
