"""
Модель динамики самолета в продольном движении.

Этот модуль содержит класс для вычисления производных состояния самолета,
балансировочных условий и нагрузочных факторов.

Принцип Single Responsibility: только динамика самолета.
Использует AerodynamicCoefficients через композицию (Dependency Inversion).
"""

import numpy as np
from .aerodynamics import AerodynamicCoefficients
from ..config.aircraft_params import AircraftParams


class AircraftDynamics:
    """
    Класс для моделирования продольной динамики самолета.
    
    Состояние: [V, alpha, q, theta]
    где:
        V - воздушная скорость (м/с)
        alpha - угол атаки (рад)
        q - угловая скорость тангажа (рад/с)
        theta - угол тангажа (рад)
    """
    
    def __init__(
        self,
        params: AircraftParams = None,
        aero_model: AerodynamicCoefficients = None
    ):
        """
        Инициализация модели динамики.
        
        Args:
            params: Параметры самолета
            aero_model: Аэродинамическая модель (Dependency Injection)
        """
        self.params = params if params is not None else AircraftParams
        self.aero = aero_model if aero_model is not None else AerodynamicCoefficients(params)
        
        # Кэшированные балансировочные условия
        self._trim_cache = None
    
    def compute_derivatives(
        self,
        t: float,
        state: np.ndarray,
        control: float,
        identified_params: np.ndarray,
        thrust: float = None
    ) -> np.ndarray:
        """
        Вычисляет производные состояния самолета.
        
        Уравнения продольного движения:
        dV/dt = (T*cos(α) - D - m*g*sin(θ-α)) / m
        dα/dt = (-L + m*g*cos(θ-α)) / (m*V) + q
        dq/dt = M / Iyy
        dθ/dt = q
        
        Args:
            t: Время (сек)
            state: Состояние [V, alpha, q, theta]
            control: Управление (delta_e - отклонение руля высоты, рад)
            identified_params: Идентифицируемые параметры [CLα, CLδe, Cmα, Cmq̄, Cmδe]
            thrust: Тяга двигателя (Н). Если None, используется балансировочная тяга
        
        Returns:
            Производные состояния [dV/dt, dα/dt, dq/dt, dθ/dt]
        """
        V, alpha, q, theta = state
        delta_e = control

        # Защита от неадекватных состояний (для устойчивости интегрирования)
        # В продольной модели деление на V критично: ограничим снизу.
        V_eff = max(float(V), 5.0)
        
        # Используем балансировочную тягу если не указана другая
        if thrust is None:
            if self._trim_cache is None:
                self._trim_cache = self.compute_trim_conditions()
            thrust = self._trim_cache[2]
        
        # Для вычисления производных предполагаем alpha_dot = 0
        # (будет пересчитано на следующем шаге)
        alpha_dot = 0.0
        
        # Получаем аэродинамические силы и момент
        L, D, M = self.aero.compute_forces_moments(
            V_eff, alpha, alpha_dot, q, delta_e, identified_params
        )
        
        # Производные состояния
        dV_dt = (thrust * np.cos(alpha) - D - self.params.MASS * self.params.G * np.sin(theta - alpha)) / self.params.MASS
        
        dalpha_dt = ((-L + self.params.MASS * self.params.G * np.cos(theta - alpha)) / 
                     (self.params.MASS * V_eff) + q)
        
        dq_dt = M / self.params.IYY
        
        dtheta_dt = q
        
        return np.array([dV_dt, dalpha_dt, dq_dt, dtheta_dt])
    
    def compute_trim_conditions(
        self,
        V_target: float = None,
        identified_params: np.ndarray = None
    ) -> tuple[float, float, float, float]:
        """
        Вычисляет балансировочные (trim) условия для горизонтального полета.
        
        В горизонтальном полете:
        - Подъемная сила = Вес
        - Тяга = Сопротивление
        - Момент тангажа = 0
        
        Args:
            V_target: Целевая скорость (м/с). Если None, используется TRIM_VELOCITY
            identified_params: Параметры для расчета. Если None, используются истинные
        
        Returns:
            tuple: (V, alpha_trim, thrust_trim, delta_e_trim)
        """
        if V_target is None:
            V_target = self.params.TRIM_VELOCITY
        
        if identified_params is None:
            identified_params = self.params.TRUE_PARAMS
        
        # Скоростной напор
        q_dyn = 0.5 * self.params.RHO * V_target ** 2
        
        # Требуемый коэффициент подъемной силы для горизонтального полета
        L_required = self.params.MASS * self.params.G
        CL_required = L_required / (q_dyn * self.params.WING_AREA)
        
        # Угол атаки из требования по подъемной силе
        # CL_required = CL0 + CLα * α + ... (остальные члены малы при балансировке)
        CL_alpha = identified_params[0]
        alpha_trim = (CL_required - self.params.CL0) / CL_alpha
        
        # Коэффициент сопротивления
        CD = self.params.CD0 + self.params.CD_ALPHA * alpha_trim
        
        # Требуемая тяга для балансировки сопротивления
        D = q_dyn * self.params.WING_AREA * CD
        thrust_trim = D
        
        # Отклонение руля для нулевого момента
        # Cm = 0 = Cm0 + Cmα * α + Cmδe * δe
        Cm_alpha = identified_params[2]
        Cm_delta_e = identified_params[4]
        delta_e_trim = -(self.params.CM0 + Cm_alpha * alpha_trim) / Cm_delta_e
        
        return V_target, alpha_trim, thrust_trim, delta_e_trim
    
    def compute_load_factor(
        self,
        state: np.ndarray,
        control: float,
        identified_params: np.ndarray
    ) -> float:
        """
        Вычисляет нормальную перегрузку (n_z).
        
        n_z = (L*cos(α) + D*sin(α)) / (m*g)
        
        Args:
            state: Состояние [V, alpha, q, theta]
            control: Управление (delta_e)
            identified_params: Идентифицируемые параметры
        
        Returns:
            Нормальная перегрузка (безразмерная, в единицах g)
        """
        V, alpha, q, theta = state
        delta_e = control
        
        alpha_dot = 0.0
        L, D, _ = self.aero.compute_forces_moments(
            V, alpha, alpha_dot, q, delta_e, identified_params
        )
        
        nz = (L * np.cos(alpha) + D * np.sin(alpha)) / (self.params.MASS * self.params.G)
        
        return nz
    
    def get_trim_state(self) -> np.ndarray:
        """
        Возвращает балансировочное состояние самолета.
        
        Returns:
            Состояние [V_trim, alpha_trim, 0, alpha_trim]
        """
        if self._trim_cache is None:
            self._trim_cache = self.compute_trim_conditions()
        
        V_trim, alpha_trim, _, _ = self._trim_cache
        
        # В балансировке: q = 0, theta = alpha (горизонтальный полет)
        return np.array([V_trim, alpha_trim, 0.0, alpha_trim])
    
    def get_trim_control(self) -> float:
        """
        Возвращает балансировочное управление.
        
        Returns:
            delta_e_trim
        """
        if self._trim_cache is None:
            self._trim_cache = self.compute_trim_conditions()
        
        return self._trim_cache[3]
    
    def get_trim_thrust(self) -> float:
        """
        Возвращает балансировочную тягу.
        
        Returns:
            thrust_trim
        """
        if self._trim_cache is None:
            self._trim_cache = self.compute_trim_conditions()
        
        return self._trim_cache[2]
    
    def print_trim_conditions(self):
        """Выводит информацию о балансировочных условиях."""
        V, alpha, thrust, delta_e = self.compute_trim_conditions()
        
        q_dyn = 0.5 * self.params.RHO * V ** 2
        CL_req = (self.params.MASS * self.params.G) / (q_dyn * self.params.WING_AREA)
        CD = self.params.CD0 + self.params.CD_ALPHA * alpha
        
        print("=" * 80)
        print(f"БАЛАНСИРОВКА ДЛЯ V = {V:.1f} м/с:")
        print("-" * 80)
        print(f"  Требуемый CL = {CL_req:.3f}")
        print(f"  Угол атаки α = {np.rad2deg(alpha):.2f}°")
        print(f"  Балансировочное δe = {np.rad2deg(delta_e):.2f}°")
        print(f"  Сопротивление CD = {CD:.4f}")
        print(f"  Требуемая тяга = {thrust:.0f} Н")
        print("=" * 80)
