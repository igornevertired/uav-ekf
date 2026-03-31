"""
Аэродинамическая модель самолета Cessna 172 Skyhawk.

Этот модуль содержит реалистичную аэродинамическую модель с учётом:
- Нелинейных эффектов (квадратичная зависимость CD от α)
- Сжимаемости воздуха (число Маха)
- Нестационарных эффектов (влияние α̇ и q)
- Эффектов срыва потока при больших углах атаки

ИСТОЧНИКИ:
- NASA CR-172501: "Identification of Light Aircraft from Flight Data"
- Stevens, Lewis "Aircraft Control and Simulation" (2003), Appendix B
- Etkin, Reid "Dynamics of Flight: Stability and Control" (1996), Chapter 4

Принцип Single Responsibility: только аэродинамические расчеты.
"""

import numpy as np
from ..config.aircraft_params import AircraftParams


class AerodynamicCoefficients:
    """
    Класс для вычисления реалистичных аэродинамических коэффициентов Cessna 172.
    
    Модель учитывает:
    1. Линейную область (малые углы атаки)
    2. Нелинейные эффекты (квадратичная зависимость CD)
    3. Нестационарные эффекты (α̇, q)
    4. Эффекты сжимаемости (поправка Прандтля-Глауэрта)
    5. Срыв потока при α > α_stall
    """
    
    def __init__(self, params: AircraftParams = None):
        """
        Инициализация аэродинамической модели.
        
        Args:
            params: Параметры самолета (если None, используются значения по умолчанию)
        """
        self.params = params if params is not None else AircraftParams
        
        # Критические углы для модели срыва
        self.alpha_stall_pos = self.params.ALPHA_STALL  # +16° (сваливание)
        self.alpha_stall_neg = -12.0 * np.pi/180        # -12° (отрицательное сваливание)
        
    def _compressibility_correction(self, V: float) -> float:
        """
        Вычисляет поправку на сжимаемость (формула Прандтля-Глауэрта).
        
        Для дозвуковых скоростей (M < 0.7):
        β = √(1 - M²)
        CL = CL_incomp / β
        
        Args:
            V: Воздушная скорость (м/с)
            
        Returns:
            Коэффициент поправки β
        """
        mach = V / self.params.SPEED_OF_SOUND
        
        if mach < 0.7:
            # Дозвуковая область - поправка Прандтля-Глауэрта
            beta = np.sqrt(1.0 - mach**2)
            return 1.0 / beta if beta > 0.1 else 1.0
        else:
            # Трансзвуковая область (не актуально для Cessna 172)
            return 1.0
    
    def _stall_model(self, alpha: float, CL_linear: float) -> float:
        """
        Модель срыва потока при больших углах атаки.
        
        Использует гладкую интерполяцию для CL при приближении к срыву:
        - Линейная область: |α| < 12°
        - Переходная область: 12° < |α| < 16°
        - Срыв: |α| > 16°
        
        Args:
            alpha: Угол атаки (рад)
            CL_linear: Линейный коэффициент подъемной силы
            
        Returns:
            CL с учетом срыва потока
        """
        # Максимальный CL для Cessna 172
        CL_max = 1.52  # из Cessna 172R POH
        CL_min = -1.20
        
        alpha_deg = np.rad2deg(alpha)
        
        if alpha > self.params.ALPHA_MAX_SAFE:
            # Переходная область к срыву
            if alpha < self.alpha_stall_pos:
                # Плавное насыщение CL
                ratio = (alpha - self.params.ALPHA_MAX_SAFE) / \
                       (self.alpha_stall_pos - self.params.ALPHA_MAX_SAFE)
                CL = CL_linear * (1 - ratio) + CL_max * ratio
            else:
                # Срыв - резкое падение CL
                CL = CL_max * (1.0 - 0.8 * ((alpha - self.alpha_stall_pos) / 
                                            (10 * np.pi/180)))
                CL = max(CL, 0.3)  # Остаточная подъемная сила
        elif alpha < self.alpha_stall_neg:
            # Отрицательный срыв
            CL = CL_min * (1.0 - 0.6 * (abs(alpha - self.alpha_stall_neg) / 
                                        (5 * np.pi/180)))
            CL = max(CL, CL_min * 0.5)
        else:
            # Линейная область
            CL = CL_linear
            
        return CL
    
    def compute_coefficients(
        self,
        V: float,
        alpha: float,
        alpha_dot: float,
        q: float,
        delta_e: float,
        identified_params: np.ndarray
    ) -> tuple[float, float, float]:
        """
        Вычисляет реалистичные аэродинамические коэффициенты.
        
        Args:
            V: Воздушная скорость (м/с)
            alpha: Угол атаки (рад)
            alpha_dot: Скорость изменения угла атаки (рад/с)
            q: Угловая скорость тангажа (рад/с)
            delta_e: Отклонение руля высоты (рад)
            identified_params: Идентифицируемые параметры [CLα, CLδe, Cmα, Cmq̄, Cmδe]
        
        Returns:
            tuple: (CL, CD, Cm) - коэффициенты подъемной силы, сопротивления и момента
        """
        # Распаковка идентифицируемых параметров
        CL_alpha, CL_delta_e, Cm_alpha, Cm_qbar, Cm_delta_e = identified_params
        
        # Защита от деления на ноль
        V_eff = max(V, 1e-3)
        
        # Безразмерные величины (приведённые к характерной длине - хорде)
        alpha_dot_bar = alpha_dot * self.params.CHORD / (2 * V_eff)
        q_bar = q * self.params.CHORD / (2 * V_eff)
        
        # Поправка на сжимаемость
        comp_correction = self._compressibility_correction(V_eff)
        
        # === КОЭФФИЦИЕНТ ПОДЪЕМНОЙ СИЛЫ ===
        # Полная модель с нестационарными членами
        CL_linear = (self.params.CL0 + 
                     CL_alpha * alpha * comp_correction +  # Основной член
                     self.params.CL_ALPHA_DOT * alpha_dot_bar +  # Нестационарность по α̇
                     self.params.CL_Q_BAR * q_bar +  # Влияние угловой скорости
                     CL_delta_e * delta_e)  # Управление рулем высоты
        
        # Применяем модель срыва
        CL = self._stall_model(alpha, CL_linear)
        
        # === КОЭФФИЦИЕНТ СОПРОТИВЛЕНИЯ ===
        # Полная нелинейная модель: CD = CD₀ + CD_α·α + CD_α²·α²
        # Источник: Etkin, Reid "Dynamics of Flight", eq. (4.4.9)
        CD_profile = self.params.CD0  # Паразитное сопротивление
        CD_induced_linear = self.params.CD_ALPHA * abs(alpha)  # Линейная часть
        CD_induced_quad = self.params.CD_ALPHA2 * alpha**2  # Квадратичная часть
        
        CD = CD_profile + CD_induced_linear + CD_induced_quad
        
        # Дополнительное сопротивление при отклонении руля высоты
        # Источник: Stevens, Lewis, Section 2.5
        CD_elevator = 0.003 * (delta_e / (np.pi/180))**2  # Квадратичная зависимость
        CD += CD_elevator
        
        # Дополнительное сопротивление при срыве
        if abs(alpha) > self.params.ALPHA_MAX_SAFE:
            CD_stall = 0.3 * ((abs(alpha) - self.params.ALPHA_MAX_SAFE) / 
                              (self.alpha_stall_pos - self.params.ALPHA_MAX_SAFE))**2
            CD += CD_stall
        
        # === КОЭФФИЦИЕНТ МОМЕНТА ТАНГАЖА ===
        # Полная модель с демпфированием
        Cm = (self.params.CM0 + 
              Cm_alpha * alpha +  # Статическая устойчивость
              self.params.CM_ALPHA_DOT * alpha_dot_bar +  # Демпфирование по α̇
              Cm_qbar * q_bar +  # Демпфирование по q (очень важно!)
              Cm_delta_e * delta_e)  # Управление
        
        # Нелинейность момента при больших углах атаки
        if abs(alpha) > self.params.ALPHA_MAX_SAFE:
            # Момент становится менее эффективным при срыве
            Cm_reduction = 0.7 * (abs(alpha) - self.params.ALPHA_MAX_SAFE) / \
                          (self.alpha_stall_pos - self.params.ALPHA_MAX_SAFE)
            Cm *= (1.0 - Cm_reduction)
        
        return CL, CD, Cm
    
    def compute_forces_moments(
        self,
        V: float,
        alpha: float,
        alpha_dot: float,
        q: float,
        delta_e: float,
        identified_params: np.ndarray
    ) -> tuple[float, float, float]:
        """
        Вычисляет реалистичные аэродинамические силы и моменты.
        
        Args:
            V: Воздушная скорость (м/с)
            alpha: Угол атаки (рад)
            alpha_dot: Скорость изменения угла атаки (рад/с)
            q: Угловая скорость тангажа (рад/с)
            delta_e: Отклонение руля высоты (рад)
            identified_params: Идентифицируемые параметры [CLα, CLδe, Cmα, Cmq̄, Cmδe]
        
        Returns:
            tuple: (L, D, M) - подъемная сила (Н), сила сопротивления (Н), 
                   момент тангажа (Н·м)
        """
        # Получаем коэффициенты
        CL, CD, Cm = self.compute_coefficients(
            V, alpha, alpha_dot, q, delta_e, identified_params
        )
        
        # Скоростной напор q = ½ρV²
        q_dyn = 0.5 * self.params.RHO * V ** 2
        
        # Силы и момент согласно стандартным формулам аэродинамики
        # Источник: Etkin, Reid "Dynamics of Flight", eq. (4.3.1-3)
        L = q_dyn * self.params.WING_AREA * CL  # Подъемная сила [Н]
        D = q_dyn * self.params.WING_AREA * CD  # Сопротивление [Н]
        M = q_dyn * self.params.WING_AREA * self.params.CHORD * Cm  # Момент тангажа [Н·м]
        
        return L, D, M
    
    def compute_dynamic_pressure(self, V: float) -> float:
        """
        Вычисляет скоростной напор.
        
        q_dyn = ½ρV²
        
        Args:
            V: Воздушная скорость (м/с)
        
        Returns:
            Скоростной напор (Па)
        """
        return 0.5 * self.params.RHO * V ** 2
    
    def compute_mach_number(self, V: float) -> float:
        """
        Вычисляет число Маха.
        
        M = V / a, где a - скорость звука
        
        Args:
            V: Воздушная скорость (м/с)
            
        Returns:
            Число Маха
        """
        return V / self.params.SPEED_OF_SOUND
    
    def compute_reynolds_number(self, V: float) -> float:
        """
        Вычисляет число Рейнольдса (для полноты модели).
        
        Re = ρVc̄/μ
        
        Args:
            V: Воздушная скорость (м/с)
            
        Returns:
            Число Рейнольдса
        """
        # Динамическая вязкость воздуха при стандартных условиях
        mu = 1.789e-5  # Па·с (при 15°C)
        
        Re = (self.params.RHO * V * self.params.CHORD) / mu
        return Re
    
    def compute_lift_coefficient_derivative(
        self,
        alpha: float,
        identified_params: np.ndarray
    ) -> float:
        """
        Вычисляет производную коэффициента подъемной силы по углу атаки.
        
        Args:
            alpha: Угол атаки (рад)
            identified_params: Идентифицируемые параметры
        
        Returns:
            dCL/dα [1/рад]
        """
        return identified_params[0]  # CLα
    
    def compute_moment_coefficient_derivative(
        self,
        alpha: float,
        identified_params: np.ndarray
    ) -> float:
        """
        Вычисляет производную коэффициента момента по углу атаки.
        
        Args:
            alpha: Угол атаки (рад)
            identified_params: Идентифицируемые параметры
        
        Returns:
            dCm/dα [1/рад]
        """
        return identified_params[2]  # Cmα
    
    def get_aerodynamic_info(self, V: float, alpha: float) -> dict:
        """
        Возвращает полную аэродинамическую информацию для анализа.
        
        Args:
            V: Воздушная скорость (м/с)
            alpha: Угол атаки (рад)
            
        Returns:
            dict: Словарь с аэродинамическими параметрами
        """
        q_dyn = self.compute_dynamic_pressure(V)
        mach = self.compute_mach_number(V)
        reynolds = self.compute_reynolds_number(V)
        comp_corr = self._compressibility_correction(V)
        
        return {
            'dynamic_pressure_Pa': q_dyn,
            'mach_number': mach,
            'reynolds_number': reynolds,
            'compressibility_correction': comp_corr,
            'alpha_deg': np.rad2deg(alpha),
            'alpha_stall_deg': np.rad2deg(self.alpha_stall_pos),
            'stall_margin_deg': np.rad2deg(self.alpha_stall_pos - alpha),
            'airspeed_kts': V * 1.944  # м/с → узлы
        }
