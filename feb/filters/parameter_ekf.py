"""
Расширенный фильтр Калмана для идентификации аэродинамических параметров.

Этот модуль реализует EKF для одновременной оценки состояния самолета 
и идентификации аэродинамических параметров.

Состояние: [V, alpha, q, theta, CLα, CLδe, Cmα, Cmq̄, Cmδe]
Измерения: [V, alpha_measured, q, theta, nz]
"""

import numpy as np
from .ekf_base import EKFBase
from ..models.aircraft_dynamics import AircraftDynamics
from ..config.aircraft_params import AircraftParams


class ParameterIdentificationEKF(EKFBase):
    """
    EKF для идентификации аэродинамических параметров самолета.
    
    Расширенное состояние включает:
    - Динамические переменные: V, alpha, q, theta
    - Идентифицируемые параметры: CLα, CLδe, Cmα, Cmq̄, Cmδe
    """
    
    def __init__(
        self,
        dt: float = 0.02,
        aircraft_model: AircraftDynamics = None,
        use_analytic_jacobians: bool = True,
        enable_adaptation: bool = True
    ):
        """
        Инициализация EKF для идентификации параметров.
        
        Args:
            dt: Временной шаг (сек)
            aircraft_model: Модель динамики самолета
            use_analytic_jacobians: Использовать аналитические якобианы (быстрее)
            enable_adaptation: Включить адаптивную фильтрацию Q и R
        """
        self.aircraft = aircraft_model if aircraft_model is not None else AircraftDynamics()
        self.params = self.aircraft.params
        self.use_analytic_jacobians = use_analytic_jacobians
        
        # Балансировка
        V_trim, alpha_trim, _, _ = self.aircraft.compute_trim_conditions()
        
        # Начальное состояние [V, alpha, q, theta, параметры]
        initial_state = np.hstack([
            V_trim, alpha_trim, 0.0, alpha_trim,
            self.params.INITIAL_PARAMS
        ])
        
        # Начальные ковариации (увеличены для параметров, чтобы легче сходиться при навигационных ошибках)
        initial_P = np.diag([
            1.0, 0.0001, 0.0001, 0.0001,   # Состояние
            0.5, 0.05, 0.05, 4.0, 0.05     # Параметры
        ])
        
        # Шум процесса (параметры: допускаем медленную подстройку)
        Q = np.diag([
            0.5, 1e-6, 1e-6, 1e-6,        # Состояние
            5e-7, 5e-7, 1e-6, 2e-6, 1e-6  # Параметры (особенно моментные)
        ])
        
        # Шум измерений [V, alpha_m, q, theta, nz]
        R = np.diag([0.25, 1e-6, 1e-6, 1e-6, 0.000025])
        
        # Инициализация базового класса
        super().__init__(
            state_dim=9,
            measurement_dim=5,
            dt=dt,
            initial_state=initial_state,
            initial_P=initial_P,
            Q=Q,
            R=R
        )
        
        # Специфичные для параметрической идентификации атрибуты
        self.param_history = []
        self.measurement_history = []
        self.control_history = []
        
        # Границы параметров
        self.param_bounds_lower, self.param_bounds_upper = self.params.get_param_bounds_array()
        
        # Адаптивная фильтрация
        if enable_adaptation:
            self.enable_adaptive_filtering(
                adaptation_rate_Q=0.001,
                adaptation_rate_R=0.001,
                window_size=100
            )
            
            # Границы для адаптации
            self.Q_min = np.diag([0.01, 1e-8, 1e-8, 1e-8, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10])
            self.Q_max = np.diag([5.0, 1e-4, 1e-4, 1e-4, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6])
            self.R_min = np.diag([0.01, 1e-8, 1e-8, 1e-8, 1e-6])
            self.R_max = np.diag([2.0, 1e-4, 1e-4, 1e-4, 0.001])
        
        print(f"ParameterEKF инициализирован. Метод якобианов: {'АНАЛИТИЧЕСКИЙ' if use_analytic_jacobians else 'ЧИСЛЕННЫЙ'}")
    
    def f_dynamics(self, x: np.ndarray, u: float, dt: float) -> np.ndarray:
        """
        Функция динамики расширенного состояния.
        
        Args:
            x: Расширенное состояние [state(4), params(5)]
            u: Управление (delta_e)
            dt: Временной шаг
        
        Returns:
            Следующее состояние
        """
        # Разделяем состояние и параметры
        state = x[:4]
        params = x[4:]
        
        # Вычисляем производные состояния
        state_dot = self.aircraft.compute_derivatives(0, state, u, params)
        
        # Интегрируем ТОЛЬКО состояние
        # Параметры ПОСТОЯННЫ (как в оригинале!)
        state_next = state + state_dot * dt
        
        return np.hstack([state_next, params])
    
    def h_measurement(self, x: np.ndarray, u: float) -> np.ndarray:
        """
        Функция измерений.
        
        Измерения: [V, alpha_measured, q, theta, nz]
        где alpha_measured = alpha + d_aoa * q/V (динамическая ошибка датчика АОА)
        
        Args:
            x: Состояние [state(4), params(5)]
            u: Управление (delta_e)
        
        Returns:
            Вектор прогнозируемых измерений
        """
        state = x[:4]
        params = x[4:]
        
        V, alpha, q, theta = state
        
        # Измеряемый угол атаки (с динамической ошибкой датчика)
        alpha_measured = alpha + self.params.AOA_SENSOR_DYNAMIC_ERROR * q / max(V, 1e-2)
        
        # Нормальная перегрузка
        nz = self.aircraft.compute_load_factor(state, u, params)
        
        return np.array([V, alpha_measured, q, theta, nz])
    
    def jacobian_F(self, x: np.ndarray, u: float) -> np.ndarray:
        """
        Вычисляет якобиан функции динамики.
        
        Args:
            x: Состояние
            u: Управление
        
        Returns:
            Матрица Якоби F (9x9)
        """
        if self.use_analytic_jacobians:
            return self._jacobian_F_analytic(x, u)
        else:
            return self._jacobian_F_numeric(x, u)
    
    def jacobian_H(self, x: np.ndarray, u: float) -> np.ndarray:
        """
        Вычисляет якобиан функции измерений.
        
        Args:
            x: Состояние
            u: Управление
        
        Returns:
            Матрица Якоби H (5x9)
        """
        if self.use_analytic_jacobians:
            return self._jacobian_H_analytic(x, u)
        else:
            return self._jacobian_H_numeric(x, u)
    
    def _jacobian_F_analytic(self, x: np.ndarray, u: float) -> np.ndarray:
        """
        Аналитическое вычисление якобиана F для дискретной системы.
        
        Returns:
            F_discrete (9x9)
        """
        V, alpha, q, theta = x[:4]
        CL_alpha, CL_delta_e, Cm_alpha, Cm_qbar, Cm_delta_e = x[4:]
        delta_e = u
        V_safe = max(V, 5.0)
        
        # Параметры
        m = self.params.MASS
        Iyy = self.params.IYY
        S = self.params.WING_AREA
        c = self.params.CHORD
        rho = self.params.RHO
        g = self.params.G
        
        # Балансировочная тяга
        Thrust = self.aircraft.get_trim_thrust()
        
        # Аэродинамические коэффициенты
        qbar = 0.5 * rho * V_safe ** 2
        CL = self.params.CL0 + CL_alpha * alpha + CL_delta_e * delta_e
        CD = self.params.CD0 + self.params.CD_ALPHA * alpha
        CD_alpha = self.params.CD_ALPHA
        
        # Матрица F для НЕПРЕРЫВНОЙ системы (9x9)
        F_cont = np.zeros((9, 9))
        
        # === Строка 1: dV/dt ===
        F_cont[0, 0] = -rho * V_safe * S * CD / m
        F_cont[0, 1] = -Thrust * np.sin(alpha) / m - qbar * S * CD_alpha / m + g * np.cos(theta - alpha)
        F_cont[0, 3] = -g * np.cos(theta - alpha)
        
        # === Строка 2: dα/dt ===
        L = qbar * S * CL
        dL_dV = rho * V_safe * S * CL
        F_cont[1, 0] = L / (m * V_safe ** 2) - dL_dV / (m * V_safe)
        
        dL_dalpha = qbar * S * CL_alpha
        F_cont[1, 1] = -dL_dalpha / (m * V_safe) - g * np.sin(theta - alpha) / V_safe
        F_cont[1, 2] = 1.0
        F_cont[1, 3] = g * np.sin(theta - alpha) / V_safe
        F_cont[1, 4] = -qbar * S * alpha / (m * V_safe)
        F_cont[1, 5] = -qbar * S * delta_e / (m * V_safe)
        
        # === Строка 3: dq/dt ===
        Cm = self.params.CM0 + Cm_alpha * alpha + Cm_qbar * (c * q / (2 * V_safe)) + Cm_delta_e * delta_e
        M = qbar * S * c * Cm
        
        dM_dV = rho * V_safe * S * c * Cm
        F_cont[2, 0] = dM_dV / Iyy
        
        dM_dalpha = qbar * S * c * Cm_alpha
        F_cont[2, 1] = dM_dalpha / Iyy
        
        dM_dq = qbar * S * c * (Cm_qbar * c / (2 * V_safe))
        F_cont[2, 2] = dM_dq / Iyy
        
        F_cont[2, 6] = qbar * S * c * alpha / Iyy
        F_cont[2, 7] = qbar * S * c * (c * q / (2 * V_safe)) / Iyy
        F_cont[2, 8] = qbar * S * c * delta_e / Iyy
        
        # === Строка 4: dθ/dt ===
        F_cont[3, 2] = 1.0
        
        # === Строки 5-9: Параметры ===
        # Параметры считаются постоянными: dp/dt = 0, поэтому в якобиане по ним единичная динамика.
        # (Ранее здесь был искусственный decay, который ухудшал согласованность модели EKF.)
        # F_cont[i,i] уже 0 для i=4..8, этого достаточно.
        
        # Дискретизация: F_discrete = I + F_cont * dt
        F_discrete = np.eye(9) + F_cont * self.dt
        
        return F_discrete
    
    def _jacobian_H_analytic(self, x: np.ndarray, u: float) -> np.ndarray:
        """
        Аналитическое вычисление якобиана H.
        
        Returns:
            H (5x9)
        """
        V, alpha, q, theta = x[:4]
        CL_alpha, CL_delta_e, Cm_alpha, Cm_qbar, Cm_delta_e = x[4:]
        delta_e = u
        
        # Параметры
        m = self.params.MASS
        S = self.params.WING_AREA
        rho = self.params.RHO
        g = self.params.G
        d_aoa = self.params.AOA_SENSOR_DYNAMIC_ERROR
        
        qbar = 0.5 * rho * V ** 2
        CL = self.params.CL0 + CL_alpha * alpha + CL_delta_e * delta_e
        CD = self.params.CD0 + self.params.CD_ALPHA * alpha
        CD_alpha = self.params.CD_ALPHA
        
        H = np.zeros((5, 9))
        
        # 1. V_m = V
        H[0, 0] = 1.0
        
        # 2. α_m = α + d_aoa * q/V
        if V > 0.1:
            H[1, 0] = -d_aoa * q / (V ** 2)
            H[1, 1] = 1.0
            H[1, 2] = d_aoa / V
        
        # 3. q_m = q
        H[2, 2] = 1.0
        
        # 4. θ_m = θ
        H[3, 3] = 1.0
        
        # 5. n_z = (L*cos(α) + D*sin(α))/(m*g)
        cos_alpha = np.cos(alpha)
        sin_alpha = np.sin(alpha)
        
        H[4, 0] = rho * V * S * (CL * cos_alpha + CD * sin_alpha) / (m * g)
        H[4, 1] = (qbar * S / (m * g)) * (
            CL_alpha * cos_alpha - CL * sin_alpha +
            CD_alpha * sin_alpha + CD * cos_alpha
        )
        H[4, 4] = qbar * S * alpha * cos_alpha / (m * g)
        H[4, 5] = qbar * S * delta_e * cos_alpha / (m * g)
        
        return H
    
    def _jacobian_F_numeric(self, x: np.ndarray, u: float) -> np.ndarray:
        """
        Численное вычисление якобиана F.
        
        Returns:
            F_discrete (9x9)
        """
        eps = 1e-5
        
        # Вспомогательная функция для производных непрерывной системы
        def f_extended_continuous(x_ext):
            state = x_ext[:4]
            params = x_ext[4:]
            state_dot = self.aircraft.compute_derivatives(0, state, u, params)
            params_dot = np.zeros(5)
            return np.hstack([state_dot, params_dot])
        
        f0 = f_extended_continuous(x)
        F_cont = np.zeros((9, 9))
        
        for i in range(9):
            x_pert = x.copy()
            x_pert[i] += eps
            f1 = f_extended_continuous(x_pert)
            F_cont[:, i] = (f1 - f0) / eps
        
        # Дискретизация
        F_discrete = np.eye(9) + F_cont * self.dt
        
        # Специальная обработка параметров
        F_discrete[4:, 4:] = 0.999 * np.eye(5)
        
        return F_discrete
    
    def _jacobian_H_numeric(self, x: np.ndarray, u: float) -> np.ndarray:
        """
        Численное вычисление якобиана H.
        
        Returns:
            H (5x9)
        """
        eps = 1e-6
        h0 = self.h_measurement(x, u)
        H = np.zeros((5, 9))
        
        for i in range(9):
            x_pert = x.copy()
            x_pert[i] += eps
            h1 = self.h_measurement(x_pert, u)
            H[:, i] = (h1 - h0) / eps
        
        return H
    
    def update(self, y_measured: np.ndarray, u: float = None):
        """
        Переопределенный метод update с ограничениями на параметры.
        
        Args:
            y_measured: Измерения
            u: Управление
        """
        # Вызываем базовый метод update
        super().update(y_measured, u)
        
        # Ограничения на состояние (V, alpha, theta) — защита при измерениях от навигации
        self.x[0] = np.clip(self.x[0], 10.0, 400.0)   # V (м/с)
        self.x[1] = np.clip(self.x[1], -0.6, 0.6)     # alpha (рад)
        self.x[3] = np.clip(self.x[3], -0.6, 0.6)     # theta (рад)
        # Применяем ограничения на параметры
        self.x[4:] = np.clip(
            self.x[4:],
            self.param_bounds_lower,
            self.param_bounds_upper
        )
        
        # Сохраняем в специальную историю
        self.param_history.append(self.x[4:].copy())
        self.measurement_history.append(y_measured.copy())
        if u is not None:
            self.control_history.append(u)
    
    def get_identified_parameters(self) -> np.ndarray:
        """Возвращает текущие оценки параметров."""
        return self.x[4:].copy()
    
    def get_aircraft_state(self) -> np.ndarray:
        """Возвращает текущее состояние самолета."""
        return self.x[:4].copy()
    
    def get_parameter_errors(self) -> np.ndarray:
        """
        Вычисляет относительные ошибки параметров в процентах.
        
        Returns:
            Массив ошибок в %
        """
        true_params = self.params.TRUE_PARAMS
        estimated_params = self.x[4:]
        errors = 100 * np.abs(estimated_params - true_params) / np.abs(true_params)
        return errors
    
    def print_parameters(self):
        """Выводит текущие оценки параметров."""
        params = self.x[4:]
        errors = self.get_parameter_errors()
        
        print("=" * 80)
        print("ТЕКУЩИЕ ОЦЕНКИ ПАРАМЕТРОВ")
        print("=" * 80)
        print(f"{'Параметр':<10} {'Истинное':<12} {'Оценка':<12} {'Ошибка':<12}")
        print("-" * 80)
        
        for i, name in enumerate(self.params.PARAM_NAMES):
            true_val = self.params.TRUE_PARAMS[i]
            est_val = params[i]
            error = errors[i]
            print(f"{name:<10} {true_val:<12.4f} {est_val:<12.4f} {error:<12.1f}%")
        
        print("=" * 80)
