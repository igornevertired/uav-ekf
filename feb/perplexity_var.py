import numpy as np
from matplotlib.axis import XTick
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib import rcParams
import time

# ========= 1. КОНСТАНТЫ =========
m = 7962.011        # МАССА самолета, кг (~8 тонн - транспортный самолет)
Iyy = 35105.145     # МОМЕНТ ИНЕРАЦИИ вокруг оси тангажа, кг·м²
S = 23.225          # ПЛОЩАДЬ КРЫЛА, м² (~спорт-лайнер)
c = 3.292           # АЭРОДИНАМИЧЕСКАЯ ХОРДА крыла, м
rho = 1.225         # ПЛОТНОСТЬ воздуха (уровень моря), кг/м³
g = 9.81            # УСКОРЕНИЕ СВОБОДНОГО ПАДЕНИЯ, м/с²
d_aoa = 2.0         # ДИНАМИЧЕСКАЯ ОШИБКА датчика АОА (α_m = α + 2·q/V)

CL0 = 0.300     # Коэффициент подъёмной силы при α=0
CLα = 3.450     # dCL/dα - подъемная сила от угла атаки [/рад]
CLα̇ = 1.120    # dCL/dα̇ - от скорости изменения угла атаки
CLq̄ = 0.000    # dCL/dq̄ - от тангажа (0 для простоты)
CLδe = 0.400    # dCL/dδe - от руля высоты

CD0 = 0.038     # Коэффициент сопротивления при α=0
CDα = 0.560     # dCD/dα - сопротивление от угла атаки

Cm0 = 0.000     # Момент тангажа при α=0
Cmα = -0.410    # dCm/dα - стабилизирующий момент [-]
Cmα̇ = -1.650    # dCm/dα̇ - от скорости изменения АОА
Cmq̄ = -4.300    # dCm/dq̄ - демпфирующий момент тангажа
Cmδe = -0.600    # dCm/dδe - эффективность руля высоты


p_true = np.array([3.450, 0.400, -0.410, -4.300, -0.600])
p_init = np.array([3.986, 0.460, -0.472, -3.655, -0.690])

param_names = ['CLα', 'CLδe', 'Cmα', 'Cmq̄', 'Cmδe']
param_descriptions = [
    'Производная коэф. подъемной силы по углу атаки',
    'Производная коэф. подъемной силы по отклонению руля',
    'Производная коэф. момента по углу атаки',
    'Производная коэф. момента по безразмерной скорости тангажа',
    'Производная коэф. момента по отклонению руля'
]

print("=" * 80)
print("ИДЕНТИФИКАЦИЯ ПАРАМЕТРОВ САМОЛЕТА - РЕАЛЬНОЕ ВРЕМЯ С ГРАФИКАМИ")
print("=" * 80)


# ========= 2. БАЛАНСИРОВКА С РАСЧЕТОМ δe_trim =========
def calculate_trim_conditions(V_target=150.0):
    q_dyn = 0.5 * rho * V_target ** 2
    L_required = m * g
    CL_required = L_required / (q_dyn * S)
    alpha_trim = (CL_required - CL0) / CLα
    CD = CD0 + CDα * alpha_trim
    D = q_dyn * S * CD
    thrust_required = D
    delta_e_trim = - (Cm0 + Cmα * alpha_trim) / Cmδe

    print(f"БАЛАНСИРОВКА ДЛЯ V = {V_target:.1f} м/с:")
    print(f"  Требуемый CL = {CL_required:.3f}")
    print(f"  Угол атаки α = {np.rad2deg(alpha_trim):.2f}°")
    print(f"  Балансировочное δe = {np.rad2deg(delta_e_trim):.2f}°")
    print(f"  Сопротивление CD = {CD:.4f}")
    print(f"  Требуемая тяга = {thrust_required:.0f} Н")

    return V_target, alpha_trim, thrust_required, delta_e_trim


V_trim, alpha_trim, Thrust_trim, delta_e_trim = calculate_trim_conditions(150.0)
Thrust = Thrust_trim


# ========= 3. ДИНАМИКА =========
def aero_coeffs(V, alpha, alpha_dot, q, delta_e, p):
    CL_alpha, CL_delta_e, Cm_alpha, Cm_qbar, Cm_delta_e = p
    V_eff = max(V, 1e-3)
    alpha_dot_bar = alpha_dot * c / (2 * V_eff)
    q_bar = q * c / (2 * V_eff)

    CL = CL0 + CL_alpha * alpha + CLα̇ * alpha_dot_bar + CLq̄ * q_bar + CL_delta_e * delta_e
    CD = CD0 + CDα * alpha
    Cm = Cm0 + Cm_alpha * alpha + Cmα̇ * alpha_dot_bar + Cm_qbar * q_bar + Cm_delta_e * delta_e
    return CL, CD, Cm


def forces_moments(V, alpha, alpha_dot, q, delta_e, p):
    CL, CD, Cm = aero_coeffs(V, alpha, alpha_dot, q, delta_e, p)
    qdyn = 0.5 * rho * V ** 2
    return qdyn * S * CL, qdyn * S * CD, qdyn * S * c * Cm


def dynamics_with_control(t, x, delta_e, p, use_controller=False):
    V, alpha, q, theta = x
    current_thrust = Thrust_trim

    alpha_dot = 0.0
    L, D, M = forces_moments(V, alpha, alpha_dot, q, delta_e, p)

    dV = (current_thrust * np.cos(alpha) - D - m * g * np.sin(theta - alpha)) / m
    dalpha = (-L + m * g * np.cos(theta - alpha)) / (m * max(V, 1e-2)) + q
    dq = M / Iyy
    dtheta = q

    return np.array([dV, dalpha, dq, dtheta]), current_thrust


def nz_from_state(x, delta_e, p):
    V, alpha, q, theta = x
    alpha_dot = 0.0
    L, D, _ = forces_moments(V, alpha, alpha_dot, q, delta_e, p)
    return (L * np.cos(alpha) + D * np.sin(alpha)) / (m * g)


# ========= 4. УПРАВЛЕНИЕ С БАЛАНСИРОВКОЙ =========
def delta_e_enhanced_profile(t):
    """Богатое возбуждение с балансировочным отклонением"""
    base = delta_e_trim

    if t < 2:
        return base
    elif t < 4:
        return base + np.deg2rad(-1.0)
    elif t < 6:
        return base + np.deg2rad(1.0)
    elif t < 7:
        return base
    elif t < 8:
        return base + np.deg2rad(-2.0)
    elif t < 10:
        return base + np.deg2rad(2.0)
    elif t < 11:
        return base
    elif t < 12:
        return base + np.deg2rad(-2.0)
    elif t < 14:
        return base + np.deg2rad(2.0)
    else:
        return base + np.deg2rad(
            1.5 * np.sin(0.5 * t) +
            0.8 * np.sin(1.5 * t) +
            0.3 * np.sin(3.0 * t)
        )


# ========= 5. EKF С АНАЛИТИЧЕСКИМ ВЫЧИСЛЕНИЕМ ЯКОБИАНОВ =========
class RealTimeEKF:
    def __init__(self, dt=0.02, use_analytic_jacobians=True):
        self.dt = dt
        self.time = 0.0
        self.use_analytic_jacobians = use_analytic_jacobians  # Флаг выбора метода

        # Начальное состояние [V, alpha, q, theta, CLα, CLδe, Cmα, Cmq̄, Cmδe]
        self.x = np.hstack([V_trim, alpha_trim, 0, alpha_trim, p_init])

        # Начальные ковариационные матрицы
        self.P = np.diag([
            1.0, 0.0001, 0.0001, 0.0001,
            0.1, 0.01, 0.01, 1.0, 0.01
        ])

        self.Q = np.diag([
            0.5, 1e-6, 1e-6, 1e-6,
            1e-8, 1e-8, 1e-8, 1e-6, 1e-8
        ])

        self.R = np.diag([0.25, 1e-6, 1e-6, 1e-6, 0.000025])

        # Адаптивные параметры
        self.innovation_history = []
        self.max_innovation_history = 100
        self.adaptation_enabled = True
        self.adaptation_rate_Q = 0.001
        self.adaptation_rate_R = 0.001

        # Ограничения
        self.Q_min = np.diag([0.01, 1e-8, 1e-8, 1e-8, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10])
        self.Q_max = np.diag([5.0, 1e-4, 1e-4, 1e-4, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6])
        self.R_min = np.diag([0.01, 1e-8, 1e-8, 1e-8, 1e-6])
        self.R_max = np.diag([2.0, 1e-4, 1e-4, 1e-4, 0.001])

        # История
        self.time_history = []
        self.state_history = []
        self.param_history = []
        self.measurement_history = []
        self.control_history = []
        self.Q_history = []
        self.R_history = []
        self.innovation_norm_history = []

        print(f"Используется метод вычисления якобианов: {'АНАЛИТИЧЕСКИЙ' if use_analytic_jacobians else 'ЧИСЛЕННЫЙ'}")

    def f_extended(self, x_ext, u):
        x, p = x_ext[:4], x_ext[4:]
        dx, _ = dynamics_with_control(0, x, u, p, use_controller=False)
        dp = np.zeros(5)  # Параметры постоянны
        return np.hstack([dx, dp])

    def h_meas(self, x_ext, u):
        x = x_ext[:4]
        V, alpha, q, theta = x
        alpha_m = alpha + d_aoa * q / max(V, 1e-2)
        nz = nz_from_state(x, u, x_ext[4:])
        return np.array([V, alpha_m, q, theta, nz])

    def jacobian_F_analytic(self, x_ext, u):
        """АНАЛИТИЧЕСКОЕ вычисление матрицы Якоби F (9x9)"""
        V, alpha, q, theta = x_ext[:4]
        CL_alpha, CL_delta_e, Cm_alpha, Cm_qbar, Cm_delta_e = x_ext[4:]
        delta_e = u

        # Предварительные вычисления
        qbar = 0.5 * rho * V ** 2
        CL = CL0 + CL_alpha * alpha + CL_delta_e * delta_e
        CD = CD0 + CDα * alpha

        # Матрица F (9x9) для НЕПРЕРЫВНОЙ системы
        F = np.zeros((9, 9))

        # === СТРОКА 1: dV/dt = (T*cos(alpha) - D - m*g*sin(theta-alpha)) / m ===
        F[0, 0] = -rho * V * S * CD / m
        F[0, 1] = -Thrust * np.sin(alpha) / m - qbar * S * CDα / m + g * np.cos(theta - alpha)
        F[0, 3] = -g * np.cos(theta - alpha)

        # === СТРОКА 2: dα/dt = (-L + m*g*cos(theta-alpha))/(m*V) + q ===
        # L = qbar * S * CL
        L = qbar * S * CL

        # ∂ᾶ/∂V = [L/(mV²)] - [∂L/∂V/(mV)] где ∂L/∂V = ρVS*CL
        dL_dV = rho * V * S * CL
        F[1, 0] = L / (m * V ** 2) - dL_dV / (m * V)

        # ∂ᾶ/∂α
        dL_dalpha = qbar * S * CL_alpha
        F[1, 1] = -dL_dalpha / (m * V) - g * np.sin(theta - alpha) / V

        F[1, 2] = 1.0  # ∂ᾶ/∂q

        # ∂ᾶ/∂θ
        F[1, 3] = g * np.sin(theta - alpha) / V

        # ∂ᾶ/∂CLα
        F[1, 4] = -qbar * S * alpha / (m * V)

        # ∂ᾶ/∂CLδe
        F[1, 5] = -qbar * S * delta_e / (m * V)

        # === СТРОКА 3: dq/dt = M / Iyy ===
        Cm = Cm0 + Cm_alpha * alpha + Cm_qbar * (c * q / (2 * max(V, 1e-3))) + Cm_delta_e * delta_e
        M = qbar * S * c * Cm

        # ∂q̇/∂V
        dM_dV = rho * V * S * c * Cm
        F[2, 0] = dM_dV / Iyy

        # ∂q̇/∂α
        dM_dalpha = qbar * S * c * Cm_alpha
        F[2, 1] = dM_dalpha / Iyy

        # ∂q̇/∂q
        dM_dq = qbar * S * c * (Cm_qbar * c / (2 * V))
        F[2, 2] = dM_dq / Iyy

        # ∂q̇/∂Cmα
        F[2, 6] = qbar * S * c * alpha / Iyy

        # ∂q̇/∂Cmq̄
        F[2, 7] = qbar * S * c * (c * q / (2 * V)) / Iyy

        # ∂q̇/∂Cmδe
        F[2, 8] = qbar * S * c * delta_e / Iyy

        # === СТРОКА 4: dθ/dt = q ===
        F[3, 2] = 1.0

        # === СТРОКИ 5-9: Параметры ===
        # КЛЮЧЕВОЕ ИСПРАВЛЕНИЕ: нужно поставить те же значения, что в численном методе
        # Но помните: это матрица Якоби для НЕПРЕРЫВНОЙ системы!
        # В численном методе после дискретизации: F_discrete[4:,4:] = 0.999*I
        # Значит для непрерывной: F_continuous[4:,4:] = (0.999 - 1)/dt = -0.001/dt

        param_rate = (0.999 - 1.0) / self.dt  # -0.001 / 0.02 = -0.05
        for i in range(4, 9):
            F[i, i] = param_rate  # -0.05

        return F

    def jacobian_H_analytic(self, x_ext, u):
        """АНАЛИТИЧЕСКОЕ вычисление матрицы Якоби H (5x9)"""
        V, alpha, q, theta = x_ext[:4]
        CL_alpha, CL_delta_e, Cm_alpha, Cm_qbar, Cm_delta_e = x_ext[4:]
        delta_e = u

        # Предварительные вычисления
        qbar = 0.5 * rho * V ** 2
        CL = CL0 + CL_alpha * alpha + CL_delta_e * delta_e
        CD = CD0 + CDα * alpha

        # Матрица H (5x9)
        H = np.zeros((5, 9))

        # 1. V_m = V
        H[0, 0] = 1.0

        # 2. α_m = α + 2*q/V
        if V > 0.1:
            H[1, 0] = -2.0 * q / (V ** 2)  # ∂α_m/∂V
            H[1, 1] = 1.0  # ∂α_m/∂α
            H[1, 2] = 2.0 / V  # ∂α_m/∂q

        # 3. q_m = q
        H[2, 2] = 1.0

        # 4. θ_m = θ
        H[3, 3] = 1.0

        # 5. n_z = (L*cosα + D*sinα)/(m*g)
        cos_alpha = np.cos(alpha)
        sin_alpha = np.sin(alpha)

        # ∂n_z/∂V
        H[4, 0] = rho * V * S * (CL * cos_alpha + CD * sin_alpha) / (m * g)

        # ∂n_z/∂α
        H[4, 1] = (qbar * S / (m * g)) * (
                CL_alpha * cos_alpha - CL * sin_alpha +
                CDα * sin_alpha + CD * cos_alpha
        )

        # ∂n_z/∂CLα
        H[4, 4] = qbar * S * alpha * cos_alpha / (m * g)

        # ∂n_z/∂CLδe
        H[4, 5] = qbar * S * delta_e * cos_alpha / (m * g)

        return H

    def jacobian_F_numeric(self, x_ext, u):
        """ЧИСЛЕННОЕ вычисление матрицы Якоби F"""
        n = 9
        eps = 1e-5

        # f0 = непрерывная производная в точке x
        f0_continuous = self.f_extended(x_ext, u)

        F_continuous = np.zeros((9, 9))

        for i in range(n):
            xh = x_ext.copy()
            xh[i] += eps

            # f1 = непрерывная производная в точке xh
            f1_continuous = self.f_extended(xh, u)

            # Численная производная НЕПРЕРЫВНОЙ системы
            F_continuous[:, i] = (f1_continuous - f0_continuous) / eps

        # Дискретизация: F_discr = I + F_continuous * dt
        F_discrete = np.eye(9) + F_continuous * self.dt

        # Параметры: специальная обработка
        F_discrete[4:, 4:] = 0.999 * np.eye(5)

        return F_discrete

    def jacobian_H_numeric(self, x_ext, u):
        """ЧИСЛЕННОЕ вычисление матрицы Якоби H"""
        m, n = 5, 9
        H = np.zeros((m, n))
        eps = 1e-6

        h0 = self.h_meas(x_ext, u)

        for i in range(n):
            xh = x_ext.copy()
            xh[i] += eps
            H[:, i] = (self.h_meas(xh, u) - h0) / eps

        return H

    def jacobian_F(self, x_ext, u):
        """Выбор метода вычисления матрицы F"""
        if self.use_analytic_jacobians:
            return self.jacobian_F_analytic(x_ext, u)
        else:
            return self.jacobian_F_numeric(x_ext, u)

    def jacobian_H(self, x_ext, u):
        """Выбор метода вычисления матрицы H"""
        if self.use_analytic_jacobians:
            return self.jacobian_H_analytic(x_ext, u)
        else:
            return self.jacobian_H_numeric(x_ext, u)

    # Остальные методы класса остаются без изменений
    def adapt_R_based_on_innovation(self, innovation, H, P):
        """Адаптация R на основе статистики инноваций"""
        if not self.adaptation_enabled or len(self.innovation_history) < 20:
            return

        try:
            innovations_array = np.array(self.innovation_history[-self.max_innovation_history:])
            if len(innovations_array) < 2:
                return

            innovation_cov = np.cov(innovations_array.T)
            theoretical_cov = H @ P @ H.T + self.R
            innovation_error = innovation_cov - theoretical_cov

            R_update = self.adaptation_rate_R * np.diag(np.diag(innovation_error))
            self.R = self.R + R_update

            self.R = np.maximum(self.R, self.R_min)
            self.R = np.minimum(self.R, self.R_max)
            self.R = (self.R + self.R.T) / 2

        except Exception as e:
            print(f"Предупреждение: ошибка адаптации R: {e}")

    def adapt_Q_based_on_residual(self, dx, dt):
        """Адаптация Q на основе несоответствия динамики"""
        if not self.adaptation_enabled:
            return

        try:
            full_dx = np.zeros(9)
            full_dx[:4] = dx

            if len(self.param_history) > 1:
                last_params = self.param_history[-1]
                prev_params = self.param_history[-2] if len(self.param_history) >= 2 else last_params
                param_changes = np.abs(last_params - prev_params) / self.dt
                full_dx[4:] = param_changes

            change_rate = np.abs(full_dx) / dt
            max_rate = 10.0
            change_rate = np.minimum(change_rate, max_rate)

            Q_update_diag = change_rate * self.adaptation_rate_Q
            current_Q_diag = np.diag(self.Q)
            new_Q_diag = current_Q_diag + Q_update_diag

            min_Q_diag = np.diag(self.Q_min)
            max_Q_diag = np.diag(self.Q_max)
            new_Q_diag = np.maximum(new_Q_diag, min_Q_diag)
            new_Q_diag = np.minimum(new_Q_diag, max_Q_diag)

            self.Q = np.diag(new_Q_diag)
            self.Q = (self.Q + self.Q.T) / 2

        except Exception as e:
            print(f"Предупреждение: ошибка адаптации Q: {e}")

    def predict(self, u):
        """Прогноз состояния"""
        x_prev = self.x.copy()

        dx, _ = dynamics_with_control(0, self.x[:4], u, self.x[4:], use_controller=False)
        self.x[:4] += dx * self.dt

        if self.adaptation_enabled:
            self.adapt_Q_based_on_residual(dx, self.dt)

        F = self.jacobian_F(self.x, u)
        self.P = F @ self.P @ F.T + self.Q

        self.P = (self.P + self.P.T) / 2
        for i in range(len(self.x)):
            if self.P[i, i] < 1e-8:
                self.P[i, i] = 1e-8

    def update(self, y_meas, u):
        """Коррекция по измерениям"""
        y_pred = self.h_meas(self.x, u)
        H = self.jacobian_H(self.x, u)
        innovation = y_meas - y_pred

        self.innovation_history.append(innovation.copy())
        self.innovation_norm_history.append(np.linalg.norm(innovation))

        if len(self.innovation_history) > self.max_innovation_history:
            self.innovation_history.pop(0)
            self.innovation_norm_history.pop(0)

        S = H @ self.P @ H.T + self.R

        try:
            K = self.P @ H.T @ np.linalg.pinv(S)
            self.x += K @ innovation
            I_KH = np.eye(len(self.x)) - K @ H
            self.P = I_KH @ self.P @ I_KH.T + K @ self.R @ K.T

            if self.adaptation_enabled and len(self.innovation_history) >= 20:
                self.adapt_R_based_on_innovation(innovation, H, self.P)
                self.R_history.append(np.diag(self.R).copy())

            # Ограничение параметров
            self.x[4] = np.clip(self.x[4], 2.5, 4.5)
            self.x[5] = np.clip(self.x[5], 0.2, 0.6)
            self.x[6] = np.clip(self.x[6], -0.6, -0.2)
            self.x[7] = np.clip(self.x[7], -6.0, -2.5)
            self.x[8] = np.clip(self.x[8], -0.8, -0.4)

        except np.linalg.LinAlgError as e:
            print(f"Предупреждение: проблема с матрицей в EKF: {e}")
            if self.adaptation_enabled:
                self.R = 1.1 * self.R
                self.R = np.minimum(self.R, self.R_max)

        self.time_history.append(self.time)
        self.state_history.append(self.x[:4].copy())
        self.param_history.append(self.x[4:].copy())
        self.measurement_history.append(y_meas.copy())
        self.control_history.append(u)

        if self.adaptation_enabled:
            self.Q_history.append(np.diag(self.Q).copy())

        max_history = 1000
        if len(self.time_history) > max_history:
            self.time_history.pop(0)
            self.state_history.pop(0)
            self.param_history.pop(0)
            self.measurement_history.pop(0)
            self.control_history.pop(0)
            if len(self.Q_history) > max_history:
                self.Q_history.pop(0)
            if len(self.R_history) > max_history:
                self.R_history.pop(0)

        self.time += self.dt
        return self.x

    def get_adaptation_info(self):
        """Получить информацию об адаптации"""
        if not self.adaptation_enabled or len(self.innovation_norm_history) < 10:
            return "Адаптация: неактивна или недостаточно данных"

        avg_innovation = np.mean(self.innovation_norm_history[-50:]) if len(
            self.innovation_norm_history) >= 50 else np.mean(self.innovation_norm_history)

        info = "Адаптация EKF:\n"
        info += f"  Размер окна инноваций: {len(self.innovation_history)}\n"
        info += f"  Средняя норма инноваций: {avg_innovation:.6f}\n"
        if hasattr(self, 'R_history') and len(self.R_history) > 0:
            info += f"  R[0,0] (шум V): {self.R[0, 0]:.6f}\n"
        if hasattr(self, 'Q_history') and len(self.Q_history) > 0:
            info += f"  Q[0,0] (шум V): {self.Q[0, 0]:.6f}"

        return info

    def test_jacobians(self):
        """Тестирование якобианов с правильным сравнением"""
        test_x = self.x.copy()
        test_u = delta_e_trim

        print("\n=== ТЕСТИРОВАНИЕ ЯКОБИАНОВ ===")

        # Аналитический (непрерывный)
        F_analytic_cont = self.jacobian_F_analytic(test_x, test_u)

        # Численный (непрерывный)
        F_numeric_disc = self.jacobian_F_numeric(test_x, test_u)  # уже дискретный
        F_numeric_cont = (F_numeric_disc - np.eye(9)) / self.dt  # преобразуем к непрерывному

        print(f"\nМакс. разница F (непрерывный): {np.max(np.abs(F_analytic_cont - F_numeric_cont)):.6f}")

        # Сравним также дискретные версии
        F_analytic_disc = np.eye(9) + F_analytic_cont * self.dt
        print(f"Макс. разница F (дискретный): {np.max(np.abs(F_analytic_disc - F_numeric_disc)):.6f}")

        print("\nЭлементы для параметров (непрерывные):")
        print("Аналитический | Численный | Разница")
        for i in range(4, 9):
            print(
                f"F[{i},{i}]: {F_analytic_cont[i, i]:.6f} | {F_numeric_cont[i, i]:.6f} | {abs(F_analytic_cont[i, i] - F_numeric_cont[i, i]):.6f}")

        # Тест матрицы H
        H_analytic = self.jacobian_H_analytic(test_x, test_u)
        H_numeric = self.jacobian_H_numeric(test_x, test_u)

        print(f"\nМакс. разница H: {np.max(np.abs(H_analytic - H_numeric)):.6f}")

# ========= 7. КЛАСС ДЛЯ ОТОБРАЖЕНИЯ В РЕАЛЬНОМ ВРЕМЕНИ =========
class RealtimePlotter:
    def __init__(self):
        plt.ion()  # Включаем интерактивный режим
        self.fig = plt.figure(figsize=(16, 10))
        self.fig.suptitle('ИДЕНТИФИКАЦИЯ ПАРАМЕТРОВ САМОЛЕТА - РЕАЛЬНОЕ ВРЕМЯ',
                          fontsize=14, fontweight='bold')

        # История истинных значений
        self.true_v_history = []
        self.true_alpha_history = []
        self.true_q_history = []
        self.true_nz_history = []

        # Создаем сетку графиков
        gs = plt.GridSpec(3, 4, hspace=0.3, wspace=0.3)

        # 1. Скорость
        self.ax1 = self.fig.add_subplot(gs[0, 0])
        self.line_v_true, = self.ax1.plot([], [], 'b-', linewidth=2, label='Истинная')
        self.line_v_meas, = self.ax1.plot([], [], 'g.', markersize=2, alpha=0.5, label='Измерения')
        self.line_v_est, = self.ax1.plot([], [], 'r--', linewidth=1.5, label='EKF')
        self.ax1.set_xlabel('Время, с')
        self.ax1.set_ylabel('V, м/с')
        self.ax1.set_title('СКОРОСТЬ ПОЛЕТА')
        self.ax1.legend(loc='upper right', fontsize=7)
        self.ax1.grid(True, alpha=0.3)

        # 2. Угол атаки
        self.ax2 = self.fig.add_subplot(gs[0, 1])
        self.line_alpha_true, = self.ax2.plot([], [], 'b-', linewidth=2, label='Истинный')
        self.line_alpha_meas, = self.ax2.plot([], [], 'g.', markersize=2, alpha=0.5, label='Измерения')
        self.line_alpha_est, = self.ax2.plot([], [], 'r--', linewidth=1.5, label='EKF')
        self.ax2.set_xlabel('Время, с')
        self.ax2.set_ylabel('α, °')
        self.ax2.set_title('УГОЛ АТАКИ')
        self.ax2.legend(loc='upper right', fontsize=7)
        self.ax2.grid(True, alpha=0.3)

        # 3. Угловая скорость
        self.ax3 = self.fig.add_subplot(gs[0, 2])
        self.line_q_true, = self.ax3.plot([], [], 'b-', linewidth=2, label='Истинная')
        self.line_q_meas, = self.ax3.plot([], [], 'g.', markersize=2, alpha=0.5, label='Измерения')
        self.line_q_est, = self.ax3.plot([], [], 'r--', linewidth=1.5, label='EKF')
        self.ax3.set_xlabel('Время, с')
        self.ax3.set_ylabel('q, °/с')
        self.ax3.set_title('УГЛОВАЯ СКОРОСТЬ ТАНГАЖА')
        self.ax3.legend(loc='upper right', fontsize=7)
        self.ax3.grid(True, alpha=0.3)

        # 4. Управление
        self.ax4 = self.fig.add_subplot(gs[0, 3])
        self.line_delta_e, = self.ax4.plot([], [], 'purple', linewidth=1.5)
        self.ax4.set_xlabel('Время, с')
        self.ax4.set_ylabel('δe, °')
        self.ax4.set_title('УПРАВЛЕНИЕ РУЛЕМ ВЫСОТЫ')
        self.ax4.grid(True, alpha=0.3)

        # 5. CLα - идентификация
        self.ax5 = self.fig.add_subplot(gs[1, 0])
        self.line_cla_true = self.ax5.axhline(y=p_true[0], color='r', linestyle='--',
                                              linewidth=2, label='Истинное', alpha=0.7)
        self.line_cla_est, = self.ax5.plot([], [], 'b-', linewidth=1.5, label='Оценка')
        self.ax5.set_xlabel('Время, с')
        self.ax5.set_ylabel('CLα')
        self.ax5.set_title('CLα - ИДЕНТИФИКАЦИЯ')
        self.ax5.legend(loc='upper right', fontsize=7)
        self.ax5.grid(True, alpha=0.3)

        # 6. Cmα - идентификация
        self.ax6 = self.fig.add_subplot(gs[1, 1])
        self.line_cma_true = self.ax6.axhline(y=p_true[2], color='r', linestyle='--',
                                              linewidth=2, label='Истинное', alpha=0.7)
        self.line_cma_est, = self.ax6.plot([], [], 'b-', linewidth=1.5, label='Оценка')
        self.ax6.set_xlabel('Время, с')
        self.ax6.set_ylabel('Cmα')
        self.ax6.set_title('Cmα - ИДЕНТИФИКАЦИЯ')
        self.ax6.legend(loc='upper right', fontsize=7)
        self.ax6.grid(True, alpha=0.3)

        # 7. Cmq̄ - идентификация
        self.ax7 = self.fig.add_subplot(gs[1, 2])
        self.line_cmq_true = self.ax7.axhline(y=p_true[3], color='r', linestyle='--',
                                              linewidth=2, label='Истинное', alpha=0.7)
        self.line_cmq_est, = self.ax7.plot([], [], 'b-', linewidth=1.5, label='Оценка')
        self.ax7.set_xlabel('Время, с')
        self.ax7.set_ylabel('Cmq̄')
        self.ax7.set_title('Cmq̄ - ИДЕНТИФИКАЦИЯ')
        self.ax7.legend(loc='upper right', fontsize=7)
        self.ax7.grid(True, alpha=0.3)

        # 8. Cmδe - идентификация
        self.ax8 = self.fig.add_subplot(gs[1, 3])
        self.line_cmde_true = self.ax8.axhline(y=p_true[4], color='r', linestyle='--',
                                               linewidth=2, label='Истинное', alpha=0.7)
        self.line_cmde_est, = self.ax8.plot([], [], 'b-', linewidth=1.5, label='Оценка')
        self.ax8.set_xlabel('Время, с')
        self.ax8.set_ylabel('Cmδe')
        self.ax8.set_title('Cmδe - ИДЕНТИФИКАЦИЯ')
        self.ax8.legend(loc='upper right', fontsize=7)
        self.ax8.grid(True, alpha=0.3)

        # 9. Ошибки параметров
        self.ax9 = self.fig.add_subplot(gs[2, 0:2])
        self.lines_error = []
        colors = ['blue', 'green', 'red', 'purple', 'orange']
        for i in range(5):
            line, = self.ax9.plot([], [], color=colors[i], linewidth=1.5, label=param_names[i])
            self.lines_error.append(line)
        self.ax9.set_xlabel('Время, с')
        self.ax9.set_ylabel('Ошибка, %')
        self.ax9.set_title('ОШИБКИ ОЦЕНКИ ПАРАМЕТРОВ')
        self.ax9.legend(loc='upper right', fontsize=7)
        self.ax9.grid(True, alpha=0.3)
        self.ax9.set_ylim([0, 50])

        # 10. Перегрузка
        self.ax10 = self.fig.add_subplot(gs[2, 2])
        self.line_nz_true, = self.ax10.plot([], [], 'b-', linewidth=2, label='Истинная')
        self.line_nz_meas, = self.ax10.plot([], [], 'g.', markersize=2, alpha=0.5, label='Измерения')
        self.ax10.set_xlabel('Время, с')
        self.ax10.set_ylabel('nz, g')
        self.ax10.set_title('НОРМАЛЬНАЯ ПЕРЕГРУЗКА')
        self.ax10.legend(loc='upper right', fontsize=7)
        self.ax10.grid(True, alpha=0.3)

        # 11. Информационная панель
        self.ax11 = self.fig.add_subplot(gs[2, 3])
        self.ax11.axis('off')
        self.info_text = self.ax11.text(0.05, 0.95, '', transform=self.ax11.transAxes,
                                        fontsize=9, verticalalignment='top',
                                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        # Счетчик обновлений
        self.update_counter = 0
        self.last_update_time = 0

    def update_plots(self, ekf, x_true, t_current):
        """Обновление всех графиков"""
        self.update_counter += 1

        # Сохраняем истинные значения
        V_true, alpha_true, q_true, theta_true = x_true
        self.true_v_history.append(V_true)
        self.true_alpha_history.append(np.rad2deg(alpha_true))
        self.true_q_history.append(np.rad2deg(q_true))

        # Обновляем только каждые 10 шагов для производительности
        if self.update_counter % 10 != 0:
            return

        current_time = time.time()
        if current_time - self.last_update_time < 0.05:  # Не чаще 10 Гц
            return

        self.last_update_time = current_time

        # Получаем данные из EKF
        if len(ekf.time_history) == 0:
            return

        time_data = np.array(ekf.time_history)

        # Ограничиваем данные для скользящего окна (последние 20 секунд)
        window_seconds = 20
        window_points = int(window_seconds / ekf.dt)
        if len(time_data) > window_points:
            time_display = time_data[-window_points:]
            true_v_display = self.true_v_history[-window_points:]
            true_alpha_display = self.true_alpha_history[-window_points:]
            true_q_display = self.true_q_history[-window_points:]
        else:
            time_display = time_data
            true_v_display = self.true_v_history
            true_alpha_display = self.true_alpha_history
            true_q_display = self.true_q_history

        # 1. Скорость - ОБНОВЛЕНО: добавляем истинные значения
        if len(true_v_display) == len(time_display):
            self.line_v_true.set_data(time_display, true_v_display)

        if len(ekf.state_history) >= len(time_display):
            v_est = [s[0] for s in ekf.state_history[-len(time_display):]]
            self.line_v_est.set_data(time_display, v_est)

        if len(ekf.measurement_history) >= len(time_display):
            v_meas = [m[0] for m in ekf.measurement_history[-len(time_display):]]
            self.line_v_meas.set_data(time_display, v_meas)

        # 2. Угол атаки - ОБНОВЛЕНО: добавляем истинные значения
        if len(true_alpha_display) == len(time_display):
            self.line_alpha_true.set_data(time_display, true_alpha_display)

        if len(ekf.state_history) >= len(time_display):
            alpha_est = [np.rad2deg(s[1]) for s in ekf.state_history[-len(time_display):]]
            self.line_alpha_est.set_data(time_display, alpha_est)

        if len(ekf.measurement_history) >= len(time_display):
            alpha_meas = [np.rad2deg(m[1]) for m in ekf.measurement_history[-len(time_display):]]
            self.line_alpha_meas.set_data(time_display, alpha_meas)

        # 3. Угловая скорость - ОБНОВЛЕНО: добавляем истинные значения
        if len(true_q_display) == len(time_display):
            self.line_q_true.set_data(time_display, true_q_display)

        if len(ekf.state_history) >= len(time_display):
            q_est = [np.rad2deg(s[2]) for s in ekf.state_history[-len(time_display):]]
            self.line_q_est.set_data(time_display, q_est)

        if len(ekf.measurement_history) >= len(time_display):
            q_meas = [np.rad2deg(m[2]) for m in ekf.measurement_history[-len(time_display):]]
            self.line_q_meas.set_data(time_display, q_meas)

        # 4. Управление
        if len(ekf.control_history) >= len(time_display):
            delta_e_display = [np.rad2deg(u) for u in ekf.control_history[-len(time_display):]]
            self.line_delta_e.set_data(time_display, delta_e_display)

        # 5-8. Параметры
        if len(ekf.param_history) >= len(time_display):
            param_data = np.array(ekf.param_history[-len(time_display):])

            # CLα
            self.line_cla_est.set_data(time_display, param_data[:, 0])
            # Cmα
            self.line_cma_est.set_data(time_display, param_data[:, 2])
            # Cmq̄
            self.line_cmq_est.set_data(time_display, param_data[:, 3])
            # Cmδe
            self.line_cmde_est.set_data(time_display, param_data[:, 4])

        # 9. Ошибки параметров
        if len(ekf.param_history) >= len(time_display):
            param_data = np.array(ekf.param_history[-len(time_display):])

            for i in range(5):
                errors = 100 * np.abs(param_data[:, i] - p_true[i]) / np.abs(p_true[i])
                errors = np.clip(errors, 0, 100)  # Ограничиваем для стабильности
                self.lines_error[i].set_data(time_display, errors)

        # 10. Перегрузка - ОБНОВЛЕНО: добавляем истинные значения
        if len(ekf.measurement_history) >= len(time_display):
            nz_meas = [m[4] for m in ekf.measurement_history[-len(time_display):]]
            self.line_nz_meas.set_data(time_display, nz_meas)

        # 11. Информационная панель
        if len(ekf.param_history) > 0:
            current_params = ekf.param_history[-1]
            errors = 100 * np.abs(current_params - p_true) / np.abs(p_true)

            info_str = f"Время: {t_current:.1f} с\n\n"
            info_str += "Текущие оценки:\n"
            for i in range(5):
                info_str += f"{param_names[i]}: {current_params[i]:.3f}\n"
                info_str += f"  Ошибка: {errors[i]:.1f}%\n"

            self.info_text.set_text(info_str)

        # Автомасштабирование осей X для всех графиков
        for ax in [self.ax1, self.ax2, self.ax3, self.ax4, self.ax5,
                   self.ax6, self.ax7, self.ax8, self.ax9, self.ax10]:
            if len(time_display) > 1:
                ax.set_xlim([time_display[0], time_display[-1]])

        # Автомасштабирование Y для некоторых графиков
        self.autoscale_y_axes(ekf, time_display, true_v_display, true_alpha_display, true_q_display)

        # Перерисовка
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def autoscale_y_axes(self, ekf, time_display, true_v_display, true_alpha_display, true_q_display):
        """Автоматическое масштабирование осей Y с учетом истинных значений"""
        if len(ekf.state_history) >= len(time_display) and len(true_v_display) == len(time_display):
            # Скорость - учитываем и истинные значения и оценки
            v_est = [s[0] for s in ekf.state_history[-len(time_display):]]
            v_min = min(min(v_est), min(true_v_display))
            v_max = max(max(v_est), max(true_v_display))
            self.ax1.set_ylim([v_min - 1, v_max + 1])

            # Угол атаки - учитываем и истинные значения и оценки
            alpha_est = [np.rad2deg(s[1]) for s in ekf.state_history[-len(time_display):]]
            alpha_min = min(min(alpha_est), min(true_alpha_display))
            alpha_max = max(max(alpha_est), max(true_alpha_display))
            self.ax2.set_ylim([alpha_min - 1, alpha_max + 1])

            # Угловая скорость - учитываем и истинные значения и оценки
            q_est = [np.rad2deg(s[2]) for s in ekf.state_history[-len(time_display):]]
            q_min = min(min(q_est), min(true_q_display))
            q_max = max(max(q_est), max(true_q_display))
            self.ax3.set_ylim([q_min - 1, q_max + 1])
        elif len(ekf.state_history) >= len(time_display):
            # Только оценки, если нет истинных значений
            v_data = [s[0] for s in ekf.state_history[-len(time_display):]]
            v_min, v_max = min(v_data), max(v_data)
            self.ax1.set_ylim([v_min - 1, v_max + 1])

            alpha_data = [np.rad2deg(s[1]) for s in ekf.state_history[-len(time_display):]]
            alpha_min, alpha_max = min(alpha_data), max(alpha_data)
            self.ax2.set_ylim([alpha_min - 1, alpha_max + 1])

            q_data = [np.rad2deg(s[2]) for s in ekf.state_history[-len(time_display):]]
            q_min, q_max = min(q_data), max(q_data)
            self.ax3.set_ylim([q_min - 1, q_max + 1])

        # Управление
        if len(ekf.control_history) >= len(time_display):
            delta_e_data = [np.rad2deg(u) for u in ekf.control_history[-len(time_display):]]
            delta_e_min, delta_e_max = min(delta_e_data), max(delta_e_data)
            self.ax4.set_ylim([delta_e_min - 2, delta_e_max + 2])

    def close(self):
        """Закрытие графиков"""
        plt.ioff()
        plt.close(self.fig)


# ========= 8. ГЕНЕРАЦИЯ ДАННЫХ ДАТЧИКОВ =========
def generate_sensor_data(V, alpha, q, theta, delta_e, p_true, add_noise=True):
    """Генерация зашумленных измерений"""
    # Уровни шума
    noise_levels = np.array([0.5, 0.001, 0.001, 0.001, 0.005])

    # Измеряемый угол атаки
    alpha_m = alpha + d_aoa * q / max(V, 1e-2)

    # Перегрузка
    nz = nz_from_state([V, alpha, q, theta], delta_e, p_true)

    measurements = np.array([V, alpha_m, q, theta, nz])

    if add_noise:
        measurements += np.random.randn(5) * noise_levels

    return measurements


# ========= 9. ОСНОВНАЯ СИМУЛЯЦИЯ В РЕАЛЬНОМ ВРЕМЕНИ =========
def simulate_realtime(T=60, dt=0.02):
    """Основная функция симуляции в реальном времени"""
    print(f"\nЗАПУСК СИМУЛЯЦИИ В РЕАЛЬНОМ ВРЕМЕНИ (T={T}с, dt={dt}с):")
    print("=" * 80)

    # Инициализация
    ekf = RealTimeEKF(dt=dt)
    plotter = RealtimePlotter()

    # Начальное состояние
    x_current = np.array([V_trim, alpha_trim, 0, alpha_trim])
    t_current = 0.0

    print("Система готова. Начинаю симуляцию...")
    print("Для остановки нажмите Ctrl+C")
    print("-" * 80)

    try:
        step = 0
        while t_current < T:
            # Управление
            delta_e = delta_e_enhanced_profile(t_current)

            # Интеграция динамики
            sol = solve_ivp(
                lambda t, x: dynamics_with_control(t, x, delta_e, p_true, use_controller=True)[0],
                [t_current, t_current + dt],
                x_current,
                method='RK45',
                max_step=dt
            )
            x_current = sol.y[:, -1]

            V, alpha, q, theta = x_current
            # controller.update_thrust(V)

            # Генерация измерений
            y_meas = generate_sensor_data(V, alpha, q, theta, delta_e, p_true, add_noise=True)

            # EKF шаг
            ekf.predict(delta_e)
            ekf.update(y_meas, delta_e)

            # Обновление графиков
            plotter.update_plots(ekf, x_current, t_current)

            # Прогресс
            step += 1
            if step % 500 == 0:  # Каждые 10 секунд
                print(f"  Время: {t_current:.1f}с")

            t_current += dt

        print("\nСИМУЛЯЦИЯ ЗАВЕРШЕНА")

    except KeyboardInterrupt:
        print("\nСИМУЛЯЦИЯ ПРЕРВАНА ПОЛЬЗОВАТЕЛЕМ")
    except Exception as e:
        print(f"\nОШИБКА В СИМУЛЯЦИИ: {e}")
        import traceback
        traceback.print_exc()
    finally:
        # Закрываем графики
        plotter.close()
        plt.show()  # Показываем финальный график

    return ekf


# ========= 10. АНАЛИЗ РЕЗУЛЬТАТОВ =========
def analyze_results(ekf):
    """Анализ результатов после симуляции"""
    print("\n" + "=" * 80)
    print("АНАЛИЗ РЕЗУЛЬТАТОВ ИДЕНТИФИКАЦИИ")
    print("=" * 80)

    if len(ekf.param_history) == 0:
        print("Нет данных для анализа!")
        return

    # Последние оценки
    final_params = ekf.param_history[-1]

    print(f"\nВРЕМЯ СИМУЛЯЦИИ: {ekf.time_history[-1]:.1f} с")
    print(f"КОЛИЧЕСТВО ШАГОВ: {len(ekf.time_history)}")
    print("-" * 80)

    print(f"\n{'Параметр':<8} {'Истинное':<10} {'Начальное':<10} {'Конечное':<10} "
          f"{'Нач.ошибка':<12} {'Кон.ошибка':<12} {'Улучшение':<12}")
    print("-" * 80)

    for i in range(5):
        param_name = param_names[i]
        true_val = p_true[i]
        init_val = p_init[i]
        final_est = final_params[i]

        final_error = 100 * abs(final_est - true_val) / abs(true_val)
        init_error = 100 * abs(init_val - true_val) / abs(true_val)
        improvement = init_error - final_error

        status = "✓" if improvement > 0 else "✗"

        print(f"{param_name:<8} {true_val:<10.4f} {init_val:<10.4f} {final_est:<10.4f} "
              f"{init_error:<12.1f}% {final_error:<12.1f}% {status}{improvement:<+11.1f}%")

    print("-" * 80)

    # Средняя ошибка
    final_errors = []
    for i in range(5):
        error = 100 * abs(final_params[i] - p_true[i]) / abs(p_true[i])
        final_errors.append(error)

    avg_error = np.mean(final_errors)
    max_error = np.max(final_errors)

    print(f"\nСТАТИСТИКА:")
    print(f"  Средняя конечная ошибка: {avg_error:.1f}%")
    print(f"  Максимальная конечная ошибка: {max_error:.1f}%")

    if avg_error < 5:
        print(f"  ОЦЕНКА: ОТЛИЧНО! ✅")
    elif avg_error < 10:
        print(f"  ОЦЕНКА: ХОРОШО! 👍")
    elif avg_error < 20:
        print(f"  ОЦЕНКА: УДОВЛЕТВОРИТЕЛЬНО 👌")
    else:
        print(f"  ОЦЕНКА: ТРЕБУЕТ УЛУЧШЕНИЯ ⚠️")


# ========= 11. ФИНАЛЬНЫЕ ГРАФИКИ =========
def plot_final_results(ekf):
    """Построение финальных графиков после симуляции"""
    if len(ekf.time_history) == 0:
        print("Нет данных для построения графиков!")
        return

    fig, axes = plt.subplots(3, 2, figsize=(14, 10))

    time_data = np.array(ekf.time_history)
    param_data = np.array(ekf.param_history)

    # 1. Все параметры
    ax1 = axes[0, 0]
    colors = ['blue', 'green', 'red', 'purple', 'orange']
    for i in range(5):
        ax1.plot(time_data, param_data[:, i], color=colors[i],
                 linewidth=1.5, label=param_names[i])
        ax1.axhline(y=p_true[i], color=colors[i], linestyle='--', alpha=0.5)
    ax1.set_xlabel('Время, с')
    ax1.set_ylabel('Значение параметра')
    ax1.set_title('ИСТОРИЯ ИДЕНТИФИКАЦИИ ПАРАМЕТРОВ')
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3)

    # 2. Ошибки
    ax2 = axes[0, 1]
    for i in range(5):
        errors = 100 * np.abs(param_data[:, i] - p_true[i]) / np.abs(p_true[i])
        ax2.plot(time_data, errors, color=colors[i],
                 linewidth=1.5, label=param_names[i])
    ax2.set_xlabel('Время, с')
    ax2.set_ylabel('Ошибка, %')
    ax2.set_title('ИСТОРИЯ ОШИБОК ОЦЕНКИ')
    ax2.legend(loc='best')
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=5, color='red', linestyle=':', alpha=0.7, label='5% порог')

    # 3. Гистограмма конечных ошибок
    ax3 = axes[1, 0]
    final_errors = []
    for i in range(5):
        error = 100 * abs(param_data[-1, i] - p_true[i]) / abs(p_true[i])
        final_errors.append(error)

    bars = ax3.bar(range(5), final_errors, color=colors, alpha=0.7)
    ax3.set_xlabel('Параметр')
    ax3.set_ylabel('Конечная ошибка, %')
    ax3.set_title('КОНЕЧНЫЕ ОШИБКИ ПАРАМЕТРОВ')
    ax3.set_xticks(range(5))
    ax3.set_xticklabels(param_names)
    ax3.grid(True, alpha=0.3, axis='y')
    ax3.axhline(y=5, color='red', linestyle='--', alpha=0.5)

    for bar, error in zip(bars, final_errors):
        height = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width() / 2, height + 0.5,
                 f'{error:.1f}%', ha='center', va='bottom', fontsize=9)

    # 4. Состояния
    ax4 = axes[1, 1]
    state_data = np.array(ekf.state_history)
    ax4.plot(time_data, state_data[:, 0], 'b-', label='V, м/с')
    ax4.set_xlabel('Время, с')
    ax4.set_ylabel('V, м/с', color='b')
    ax4.tick_params(axis='y', labelcolor='b')
    ax4.grid(True, alpha=0.3)

    ax4_twin = ax4.twinx()
    ax4_twin.plot(time_data, np.rad2deg(state_data[:, 1]), 'r-', label='α, °', alpha=0.7)
    ax4_twin.set_ylabel('α, °', color='r')
    ax4_twin.tick_params(axis='y', labelcolor='r')

    ax4.set_title('СКОРОСТЬ И УГОЛ АТАКА (оценки EKF)')

    # 5. Управление
    ax5 = axes[2, 0]
    control_data = np.array(ekf.control_history)
    ax5.plot(time_data, np.rad2deg(control_data), 'purple', linewidth=1)
    ax5.set_xlabel('Время, с')
    ax5.set_ylabel('δe, °')
    ax5.set_title('ИСТОРИЯ УПРАВЛЕНИЯ')
    ax5.grid(True, alpha=0.3)

    # # 6. Перегрузка
    # ax6 = axes[2, 1]
    # if len(ekf.measurement_history) > 0:
    #     meas_data = np.array(ekf.measurement_history)
    #     ax6.plot(time_data, meas_data[:, 4], 'g-', linewidth=1, alpha=0.7, label='Измерения')
    # ax6.set_xlabel('Время, с')
    # ax6.set_ylabel('nz, g')
    # ax6.set_title('НОРМАЛЬНАЯ ПЕРЕГРУЗКА')
    # ax6.legend()
    # ax6.grid(True, alpha=0.3)

    plt.suptitle('ФИНАЛЬНЫЙ АНАЛИЗ РЕЗУЛЬТАТОВ ИДЕНТИФИКАЦИИ',
                 fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.show()



# ========= 12. ОСНОВНАЯ ФУНКЦИЯ =========
def main():
    print("\n" + "=" * 80)
    print("ИДЕНТИФИКАЦИЯ ПАРАМЕТРОВ САМОЛЕТА В РЕАЛЬНОМ ВРЕМЕНИ")
    print("=" * 80)
    use_analytic = False  # True - аналитический, False - численный

    try:
        # Запуск симуляции в реальном времени
        start_time = time.time()
        ekf = RealTimeEKF(dt=0.02, use_analytic_jacobians=use_analytic)
        ekf.test_jacobians()
        plotter = RealtimePlotter()

        # Начальное состояние
        x_current = np.array([V_trim, alpha_trim, 0, alpha_trim])
        t_current = 0.0

        print("Система готова. Начинаю симуляцию...")
        print("Для остановки нажмите Ctrl+C")
        print("-" * 80)

        step = 0
        T = 60  # время симуляции

        while t_current < T:
            # Управление
            delta_e = delta_e_enhanced_profile(t_current)

            # Интеграция динамики
            sol = solve_ivp(
                lambda t, x: dynamics_with_control(t, x, delta_e, p_true, use_controller=False)[0],
                [t_current, t_current + 0.02],
                x_current,
                method='RK45',
                max_step=0.02
            )
            x_current = sol.y[:, -1]

            V, alpha, q, theta = x_current

            # Генерация измерений
            y_meas = generate_sensor_data(V, alpha, q, theta, delta_e, p_true, add_noise=True)

            # EKF шаг
            ekf.predict(delta_e)
            ekf.update(y_meas, delta_e)

            # Обновление графиков
            plotter.update_plots(ekf, x_current, t_current)

            # Прогресс
            step += 1
            if step % 500 == 0:
                print(f"  Время: {t_current:.1f}с")

            t_current += 0.02

        print("\nСИМУЛЯЦИЯ ЗАВЕРШЕНА")

        # Анализ результатов
        analyze_results(ekf)

        # Финальные графики
        print("\nПостроение финальных графиков...")
        plot_final_results(ekf)

        sim_time = time.time() - start_time
        print(f"\nОбщее время выполнения: {sim_time:.1f} секунд")

    except KeyboardInterrupt:
        print("\nСИМУЛЯЦИЯ ПРЕРВАНА ПОЛЬЗОВАТЕЛЕМ")
    except Exception as e:
        print(f"\nКРИТИЧЕСКАЯ ОШИБКА: {e}")
        import traceback
        traceback.print_exc()
    finally:
        plt.ioff()
        plt.show()

    print("\n" + "=" * 80)
    print("ПРОГРАММА ЗАВЕРШЕНА")
    print("=" * 80)


if __name__ == "__main__":
    main()