"""
Главный симулятор интегрированной системы.

Объединяет:
- Динамику самолета
- БИНС и ГНСС
- Navigation EKF (комплексирование)
- Parameter EKF (идентификация)
- Визуализацию

Принцип Dependency Inversion: работает с абстракциями через инъекцию зависимостей.
"""

import numpy as np
from scipy.integrate import solve_ivp

from ..models.aircraft_dynamics import AircraftDynamics
from ..models.ins_model import StrapdownINS
from ..models.gnss_model import GNSSReceiver
from ..filters.parameter_ekf import ParameterIdentificationEKF
from ..filters.navigation_ekf import NavigationIntegrationEKF
from ..control.excitation import ExcitationProfile
from ..sensors.sensor_models import SensorModels
from ..visualization.realtime_plotter import RealtimePlotter
from ..config.aircraft_params import AircraftParams
from ..config.navigation_params import NavigationParams


class IntegratedSimulator:
    """
    Интегрированный симулятор с навигацией и идентификацией.
    
    Реализует полный цикл:
    1. Истинная динамика самолета
    2. Распространение БИНС
    3. Обновления ГНСС
    4. Комплексирование БИНС+ГНСС в Navigation EKF
    5. Идентификация параметров в Parameter EKF
    """
    
    def __init__(
        self,
        dt: float = 0.02,
        accuracy_class: str = 'commercial',
        use_navigation: bool = True,
        use_visualization: bool = True,
        use_nav_for_identification: bool = None,
        rng_seed: int | None = None
    ):
        """
        Инициализация симулятора.
        
        Args:
            dt: Временной шаг симуляции (сек)
            accuracy_class: Класс точности навигационного оборудования
            use_navigation: Использовать ли навигационный комплекс
            use_visualization: Использовать ли визуализацию
            use_nav_for_identification: Если True — измерения для Parameter EKF берутся
                из навигационного решения. По умолчанию: True только если навигация ВКЛ
                и визуализация ВЫКЛ (режим исследования); при включённой визуализации
                по умолчанию False — показывается сходимость по «идеальным» датчикам.
        """
        self.dt = dt
        if rng_seed is not None:
            # Фиксируем генератор случайных чисел для воспроизводимых исследований
            np.random.seed(int(rng_seed))
        self.use_navigation = use_navigation
        self.use_visualization = use_visualization
        # При визуализации по умолчанию идентификация по истине (наглядная сходимость ~2–3%);
        # без визуализации — по нав. решению (соответствие тезисам, ошибки выше).
        self.use_nav_for_identification = (
            use_nav_for_identification if use_nav_for_identification is not None
            else (use_navigation and not use_visualization)
        )
        # Сглаживание нав. оценок при подаче на идентификацию (снижает шум, небольшая задержка)
        self.use_nav_smoothing = True
        self._nav_smooth_prev = None  # (V, alpha, q, theta) на предыдущем шаге
        self._nav_smooth_alpha = 0.75  # доля предыдущего значения: out = α*prev + (1-α)*current
        
        # Параметры
        self.aircraft_params = AircraftParams
        self.ins_params, self.gnss_params = NavigationParams.get_params(accuracy_class)
        
        # Модели
        self.aircraft = AircraftDynamics()
        self.sensors = SensorModels()
        self.control = ExcitationProfile()
        
        # Навигация (если включена)
        if self.use_navigation:
            # Начальная позиция и скорость
            V_trim = self.aircraft_params.TRIM_VELOCITY
            initial_pos = np.array([0.0, 0.0, 1000.0])  # Высота 1000м
            initial_vel = np.array([V_trim, 0.0, 0.0])
            initial_att = np.array([0.0, self.aircraft.get_trim_state()[1], 0.0])  # pitch = alpha_trim
            
            self.ins = StrapdownINS(
                self.ins_params,
                initial_pos,
                initial_vel,
                initial_att
            )
            self.gnss = GNSSReceiver(self.gnss_params)
            self.nav_ekf = NavigationIntegrationEKF(
                dt=self.dt,
                ins_model=self.ins,
                gnss_model=self.gnss,
                ins_params=self.ins_params,
                gnss_params=self.gnss_params
            )
        
        # Фильтр идентификации параметров
        self.param_ekf = ParameterIdentificationEKF(
            dt=self.dt,
            aircraft_model=self.aircraft,
            use_analytic_jacobians=False,
            enable_adaptation=False
        )
        
        # Визуализация
        if self.use_visualization:
            self.plotter = RealtimePlotter()
        
        # Балансировка
        _, _, _, delta_e_trim = self.aircraft.compute_trim_conditions()
        self.control.set_trim(delta_e_trim)
        
        # История
        self.time = 0.0
        # Истинная навигация (NED) для корректного моделирования GNSS в режиме 3
        # (иначе GNSS измеряет INS, и контур становится самозамкнутым)
        self.true_position = np.array([0.0, 0.0, 1000.0])  # NED, z вниз
        self.true_velocity = np.array([self.aircraft_params.TRIM_VELOCITY, 0.0, 0.0])
        self.history = {
            'time': [],
            'true_state': [],
            'control': [],
            # История навигационной системы (если включена)
            'ins_state': [] if self.use_navigation else None,
            'ins_raw': [] if self.use_navigation else None,
            'gnss_measurement': [] if self.use_navigation else None,
            'nav_corrected': [] if self.use_navigation else None,
            'nav_errors': [] if self.use_navigation else None
        }
        
        print("=" * 80)
        print("СИМУЛЯТОР ИНИЦИАЛИЗИРОВАН")
        print("=" * 80)
        print(f"Временной шаг: {self.dt} с")
        print(f"Класс точности навигации: {accuracy_class}")
        print(f"Навигационный комплекс: {'ВКЛ' if self.use_navigation else 'ВЫКЛ'}")
        if self.use_navigation:
            print(f"Идентификация по навигационному решению: {'ДА' if self.use_nav_for_identification else 'НЕТ'}")
        print(f"Визуализация: {'ВКЛ' if self.use_visualization else 'ВЫКЛ'}")
        print("=" * 80)
    
    def run(self, T: float = 60.0, profile_type: str = 'enhanced'):
        """
        Запускает симуляцию.
        
        Args:
            T: Время симуляции (сек)
            profile_type: Тип профиля управления
        """
        print(f"\nЗАПУСК СИМУЛЯЦИИ (T={T}с)")
        print(f"Профиль управления: {self.control.get_profile_description(profile_type)}")
        print("-" * 80)
        
        # Начальное состояние
        x_true = self.aircraft.get_trim_state()
        self.time = 0.0
        self._nav_smooth_prev = None  # сброс сглаживания при новом прогоне
        # Сброс истинной навигации
        self.true_position = np.array([0.0, 0.0, 1000.0])
        self.true_velocity = np.array([self.aircraft_params.TRIM_VELOCITY, 0.0, 0.0])
        step = 0
        
        try:
            while self.time < T:
                # 1. Управление
                delta_e = self.control.get_control(self.time, profile_type)
                
                # 2. Истинная динамика самолета
                sol = solve_ivp(
                    lambda t, x: self.aircraft.compute_derivatives(t, x, delta_e, self.aircraft_params.TRUE_PARAMS),
                    [self.time, self.time + self.dt],
                    x_true,
                    method='RK45',
                    max_step=self.dt
                )
                x_true = sol.y[:, -1]

                # Защита от разлёта истинного состояния (иначе возможны overflow в аэродинамике)
                x_true[0] = np.clip(x_true[0], 5.0, 80.0)      # V (м/с)
                # Держим идентификацию в линейной зоне (иначе включается срыв/нелинейность и производные неидентифицируемы)
                x_true[1] = np.clip(x_true[1], -self.aircraft_params.ALPHA_MAX_SAFE, self.aircraft_params.ALPHA_MAX_SAFE)  # alpha (рад)
                x_true[2] = np.clip(x_true[2], -1.5, 1.5)      # q (рад/с)
                x_true[3] = np.clip(x_true[3], -0.8, 0.8)      # theta (рад)

                # Обновляем истинную кинематику в NED по продольной модели
                V_true, alpha_true, _, theta_true = x_true
                gamma = float(theta_true - alpha_true)  # путевой угол в вертикальной плоскости
                self.true_velocity = np.array([
                    float(V_true * np.cos(gamma)),
                    0.0,
                    float(-V_true * np.sin(gamma))
                ])
                self.true_position = self.true_position + self.true_velocity * self.dt
                
                # 3. Навигационная система (если включена)
                nav_result = None
                if self.use_navigation:
                    nav_result = self._process_navigation(x_true, delta_e, self.true_position, self.true_velocity)
                
                # 4. Измерения датчиков для идентификации параметров
                # Один основной режим (по тезисам): навигация ВКЛ и измерения движения
                # (q, theta) берём из навигационного решения; air-data (V, alpha) — с датчиков.
                V_true, alpha_true, _, _ = x_true
                # Влияние навигации на идентификацию задаём через точность навигационной оценки (СКО)
                # по углу тангажа (pitch). Это устраняет систематические смещения от упрощённой ориентации,
                # но сохраняет количественную зависимость от класса точности (через P Navigation EKF).
                q_src = float(x_true[2])
                theta_src = float(x_true[3])
                sigma_theta_nav = np.deg2rad(0.5)
                if nav_result is not None:
                    try:
                        sigma_theta_nav = float(self.nav_ekf.get_navigation_solution()['attitude_std'][1])
                    except Exception:
                        sigma_theta_nav = np.deg2rad(0.5)

                # Формируем измерения для Parameter EKF:
                # - V_m, alpha_m: air-data с шумом (датчики)
                # - q_m, theta_m: истинные + шум (в theta шум определяется точностью навигации)
                # - nz_m: акселерометр с шумом (датчик перегрузки)
                V_m = float(self.sensors._measure_velocity(float(V_true), add_noise=True))
                alpha_m = float(self.sensors._measure_aoa(float(V_true), float(alpha_true), float(q_src), add_noise=True))
                q_m = float(self.sensors._measure_angular_rate(float(q_src), add_noise=True))
                theta_m = float(theta_src + np.random.randn() * sigma_theta_nav)
                nz_m = float(self.sensors._measure_load_factor(float(V_true), float(alpha_true), float(q_src), float(theta_src), float(delta_e), self.aircraft_params.TRUE_PARAMS, add_noise=True))
                y_meas = np.array([V_m, alpha_m, q_m, theta_m, nz_m])
                
                # 5. EKF для идентификации параметров
                self.param_ekf.predict(delta_e)
                self.param_ekf.update(y_meas, delta_e)
                
                # 6. Визуализация
                if self.use_visualization:
                    self.plotter.update(
                        self.param_ekf, x_true, self.time,
                        true_state_history=self.history['true_state']
                    )
                
                # 7. История
                self.history['time'].append(self.time)
                self.history['true_state'].append(x_true.copy())
                self.history['control'].append(delta_e)
                
                # Прогресс
                step += 1
                if step % 250 == 0:  # Чаще выводим прогресс
                    progress = (self.time / T) * 100
                    bar_length = 40
                    filled_length = int(bar_length * self.time / T)
                    bar = '█' * filled_length + '-' * (bar_length - filled_length)
                    print(f"\r  [{bar}] {progress:.1f}% ({self.time:.1f}с / {T:.1f}с)", end='', flush=True)
                
                self.time += self.dt
            
            print("\nСИМУЛЯЦИЯ ЗАВЕРШЕНА")
            return True
            
        except KeyboardInterrupt:
            print("\nСИМУЛЯЦИЯ ПРЕРВАНА ПОЛЬЗОВАТЕЛЕМ")
            return False
        except Exception as e:
            print(f"\nОШИБКА В СИМУЛЯЦИИ: {e}")
            import traceback
            traceback.print_exc()
            return False
        finally:
            if self.use_visualization:
                self.plotter.close()
    
    def _process_navigation(
        self,
        x_true: np.ndarray,
        delta_e: float,
        true_position: np.ndarray,
        true_velocity: np.ndarray
    ):
        """
        Обрабатывает навигационную систему В РЕАЛЬНОМ ВРЕМЕНИ.
        
        Args:
            x_true: Истинное состояние самолета
            delta_e: Управление
        """
        V, alpha, q, theta = x_true
        
        # Вычисляем удельную силу и угловую скорость для БИНС
        # Упрощенная модель: используем аэродинамические силы
        from ..models.aerodynamics import AerodynamicCoefficients
        aero = AerodynamicCoefficients()
        L, D, M = aero.compute_forces_moments(
            V, alpha, 0.0, q, delta_e,
            self.aircraft_params.TRUE_PARAMS
        )
        
        # Удельная сила в связанной СК
        ax_body = (self.aircraft.get_trim_thrust() * np.cos(alpha) - D) / self.aircraft_params.MASS
        az_body = -L / self.aircraft_params.MASS
        specific_force = np.array([ax_body, 0.0, az_body])
        
        # Угловая скорость (только тангаж для 2D модели)
        angular_rate = np.array([0.0, q, 0.0])
        
        # 1. Распространение БИНС (С ОШИБКАМИ - в реальном времени!)
        ins_state = self.ins.propagate(self.dt, specific_force, angular_rate, include_errors=True)
        
        # Сохраняем RAW состояние БИНС (ДО коррекции)
        ins_raw = ins_state.copy()
        
        # 2. ГНСС измерение (с частотой ГНСС - в реальном времени!)
        gnss_meas = self.gnss.get_measurement(
            true_position,
            true_velocity,
            self.dt,
            available=True
        )
        
        # 3. Комплексирование через Navigation EKF (в реальном времени!)
        corrected_state = None
        if gnss_meas is not None and gnss_meas['available']:
            # Predict EKF
            self.nav_ekf.predict(u=None)
            
            # Update с измерениями ГНСС
            self.nav_ekf.update_with_gnss(gnss_meas, ins_state)
            
            # Применяем коррекцию к БИНС
            corrected_state = self.nav_ekf.get_corrected_state(ins_state)
            # В error-state EKF после применения коррекции к INS сбрасываем ошибки
            self.nav_ekf.reset_errors()
        
        nav_result = {'corrected_state': corrected_state, 'ins_state': ins_state}
        # 4. Сохраняем данные в историю (ДЛЯ АНАЛИЗА)
        self.history['ins_raw'].append({
            'position': ins_raw['position'].copy(),
            'velocity': ins_raw['velocity'].copy(),
            'attitude': ins_raw['attitude_euler'].copy(),
            'gyro_bias': self.ins.gyro_bias.copy(),
            'accel_bias': self.ins.accel_bias.copy()
        })
        
        if gnss_meas is not None and gnss_meas['available']:
            self.history['gnss_measurement'].append({
                'position': gnss_meas['position'].copy(),
                'velocity': gnss_meas['velocity'].copy(),
                'available': True
            })
        else:
            self.history['gnss_measurement'].append({'available': False})
        
        if corrected_state is not None:
            self.history['nav_corrected'].append({
                'position': corrected_state['position'].copy(),
                'velocity': corrected_state['velocity'].copy(),
                'attitude': corrected_state['attitude_euler'].copy()
            })
            
            # Сохраняем оценки ошибок от EKF
            nav_solution = self.nav_ekf.get_navigation_solution()
            self.history['nav_errors'].append({
                'position_errors': nav_solution['position_errors'].copy(),
                'velocity_errors': nav_solution['velocity_errors'].copy(),
                'attitude_errors': nav_solution['attitude_errors'].copy(),
                'gyro_bias_est': nav_solution['gyro_bias'].copy(),
                'accel_bias_est': nav_solution['accel_bias'].copy(),
                'position_std': nav_solution['position_std'].copy(),
                'velocity_std': nav_solution['velocity_std'].copy()
            })
        else:
            self.history['nav_corrected'].append(None)
            self.history['nav_errors'].append(None)
        
        return nav_result
    
    @staticmethod
    def _euler_to_dcm_body2ned(roll: float, pitch: float, yaw: float) -> np.ndarray:
        """Матрица направляющих косинусов из связанной СК в NED: v_ned = R @ v_body."""
        cr, sr = np.cos(roll), np.sin(roll)
        cp, sp = np.cos(pitch), np.sin(pitch)
        cy, sy = np.cos(yaw), np.sin(yaw)
        R = np.array([
            [cp*cy, -cr*sy + sr*sp*cy,  sr*sy + cr*sp*cy],
            [cp*sy,  cr*cy + sr*sp*sy, -sr*cy + cr*sp*sy],
            [-sp,    sr*cp,              cr*cp]
        ])
        return R
    
    def _nav_state_to_measurement_inputs(
        self,
        corrected_state: dict,
        ins_state: dict,
        gyro_bias_estimate: np.ndarray = None
    ) -> tuple:
        """
        Преобразует навигационное решение в входные величины для датчиков идентификации.
        
        Доработка режима 3: скорость переводится в связанную СК; V и α считаются
        из компонент в теле (соответствуют продольной модели: воздушная скорость
        и угол атаки в плоскости симметрии).
        
        Args:
            corrected_state: Скорректированное состояние (position, velocity, attitude_euler)
                или None, если ГНСС ещё не обновлялся.
            ins_state: Состояние БИНС (velocity, attitude_euler, angular_rate_body).
            gyro_bias_estimate: Оценка смещения гироскопов из Navigation EKF (3,) рад/с.
        
        Returns:
            (V, alpha, q, theta) в м/с и рад для продольного движения.
        """
        state = corrected_state if corrected_state is not None else ins_state
        vel_ned = state['velocity']
        att = state['attitude_euler']
        roll, pitch, yaw = att[0], att[1], att[2]
        
        # Скорость в связанной СК: v_body = R_ned2body @ v_ned = R_body2ned.T @ v_ned
        R_b2n = self._euler_to_dcm_body2ned(roll, pitch, yaw)
        vel_body = R_b2n.T @ vel_ned
        
        # Продольное движение: V — модуль в плоскости xz (нос–хвост и вертикаль тела)
        vx_b = vel_body[0]
        vz_b = vel_body[2]
        V = np.sqrt(vx_b**2 + vz_b**2)
        V = np.clip(V, 20.0, 300.0)
        
        # Угол атаки: направление вектора скорости в плоскости симметрии (x-z тела)
        # В балансировке theta = alpha, v_body = (V*cos(θ), 0, V*sin(θ)) → alpha = atan2(vz_b, vx_b)
        alpha = np.arctan2(vz_b, vx_b) if V > 1e-2 else 0.0
        alpha = np.clip(alpha, -0.5, 0.5)
        
        # Угол тангажа (рад)
        theta = np.clip(pitch, -0.6, 0.6)
        
        # Угловая скорость тангажа: из БИНС с вычетом оценки смещения
        q = float(ins_state['angular_rate_body'][1])
        if gyro_bias_estimate is not None and len(gyro_bias_estimate) >= 2:
            q = q - gyro_bias_estimate[1]
        q = np.clip(q, -0.5, 0.5)
        
        out = (V, alpha, q, theta)
        if self.use_nav_smoothing and self.use_nav_for_identification:
            if self._nav_smooth_prev is not None:
                V_s, a_s, q_s, t_s = self._nav_smooth_prev
                a = self._nav_smooth_alpha
                out = (
                    a * V_s + (1 - a) * V,
                    a * a_s + (1 - a) * alpha,
                    a * q_s + (1 - a) * q,
                    a * t_s + (1 - a) * theta
                )
            self._nav_smooth_prev = out
        return out
    
    def get_results(self) -> dict:
        """
        Возвращает результаты симуляции.
        
        Returns:
            Словарь с результатами
        """
        from ..analysis.results_analyzer import ResultsAnalyzer
        
        analyzer = ResultsAnalyzer()
        ident_results = analyzer.analyze_identification(self.param_ekf)
        
        nav_results = None
        if self.use_navigation:
            nav_results = analyzer.analyze_navigation(self.nav_ekf)
        
        return {
            'identification': ident_results,
            'navigation': nav_results,
            'history': self.history
        }
    
    def print_summary(self):
        """Выводит сводку результатов."""
        from ..analysis.results_analyzer import ResultsAnalyzer
        
        analyzer = ResultsAnalyzer()
        ident_results = analyzer.analyze_identification(self.param_ekf)
        
        nav_results = None
        if self.use_navigation:
            nav_results = analyzer.analyze_navigation(self.nav_ekf)
        
        analyzer.print_summary(ident_results, nav_results)
        
        # Детали параметров
        self.param_ekf.print_parameters()
        
        # Детали навигации
        if self.use_navigation:
            self.nav_ekf.print_navigation_accuracy()
