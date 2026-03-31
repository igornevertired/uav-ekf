# Справочник реализации: что и как сделано в коде

Документ описывает векторы состояния, матрицы, модели БИНС/ГНСС, шаги предсказания/обновления и общий поток данных.

**Простое пошаговое объяснение «от А до Я»** — в файле **[WALKTHROUGH_FROM_A_TO_Z.md](WALKTHROUGH_FROM_A_TO_Z.md)** (что происходит в коде на каждом шаге, простым языком и со всей нужной информацией).

---

## 1. Общая схема: что делается в коде

### 1.1 Главный цикл симуляции (`simulation/simulator.py`)

На каждом шаге по времени (dt = 0.02 с) выполняется:

```
1. Задать управление δe = control.get_control(t, profile_type)

2. Проинтегрировать истинную динамику самолёта (RK45):
   x_true = [V, α, q, θ]  ← обновляется из уравнений движения

3. Если включена навигация:
   a) По x_true вычислить удельную силу и угловую скорость в связанной СК
   b) БИНС: ins.propagate(dt, specific_force, angular_rate) → ins_state
   c) ГНСС: gnss.get_measurement(true_pos, true_vel, dt) → gnss_meas (если наступил период ГНСС)
   d) Navigation EKF: predict → update_with_gnss(gnss_meas, ins_state)
   e) Скорректированное состояние: nav_ekf.get_corrected_state(ins_state)

4. Датчики для идентификации: при включённой навигации и флаге `use_nav_for_identification` (по умолчанию True) входы для датчиков (V, α, q, θ) берутся из навигационного решения: `_nav_state_to_measurement_inputs(corrected_state, ins_state)` → затем `y_meas = sensors.generate_sensor_data(V_nav, alpha_nav, q_nav, theta_nav, δe, TRUE_PARAMS)`. Иначе — из истинного состояния: `y_meas = sensors.generate_sensor_data(V, α, q, θ, δe, TRUE_PARAMS)`.

5. Parameter EKF: predict(δe) → update(y_meas, δe)

6. Визуализация (если включена): plotter.update(param_ekf, x_true, t)

7. Запись в history
```

Истинная динамика считается по эталонным параметрам `TRUE_PARAMS`. Оценка параметров ведётся в Parameter EKF; навигация (БИНС+ГНСС+Navigation EKF) при включённом флаге считается параллельно и в результаты идентификации не подмешивается (влияние — через качество полётной информации в общем эксперименте).

---

## 2. Динамика самолёта и идентификация (Parameter EKF)

### 2.1 Вектор состояния (Parameter EKF)

Расширенное состояние **x** размерности **9**:

| Индекс | Переменная | Размерность | Описание |
|--------|------------|-------------|----------|
| 0 | V | 1 | Воздушная скорость (м/с) |
| 1 | α | 1 | Угол атаки (рад) |
| 2 | q | 1 | Угловая скорость тангажа (рад/с) |
| 3 | θ | 1 | Угол тангажа (рад) |
| 4 | CLα | 1 | Коэффициент подъёмной силы по α (1/рад) |
| 5 | CLδe | 1 | Коэффициент подъёмной силы по δe (1/рад) |
| 6 | Cmα | 1 | Момент тангажа по α (1/рад) |
| 7 | Cmq̄ | 1 | Момент тангажа по безразмерной q (1/рад) |
| 8 | Cmδe | 1 | Момент тангажа по δe (1/рад) |

Итого: **x = [V, α, q, θ, CLα, CLδe, Cmα, Cmq̄, Cmδe]ᵀ**, размер **9×1**.

### 2.2 Динамика (непрерывная модель)

Управление: **u = δe** (отклонение руля высоты, рад).

Производные состояния (используются в `f_dynamics` через `aircraft.compute_derivatives`):

- **dV/dt** = (T cos α − D − mg sin(θ−α)) / m  
- **dα/dt** = (−L + mg cos(θ−α)) / (mV) + q  
- **dq/dt** = M / Iyy  
- **dθ/dt** = q  

L, D, M считаются в аэродинамической модели по V, α, q, δe и **текущим оценкам параметров** (x[4:9]). Параметры в непрерывной модели **не меняются**: d(params)/dt = 0.

В коде (`parameter_ekf.py`):

- `f_dynamics`: по текущему **x** и **u** вызывается `compute_derivatives` только для первых четырёх переменных; затем состояние интегрируется методом Эйлера: `state_next = state + state_dot * dt`, а параметры просто копируются: `return np.hstack([state_next, params])`.

### 2.3 Дискретный якобиан F (Parameter EKF)

Используется дискретизация: **F_disc = I + F_cont · dt**.

**F_cont** (9×9) — якобиан непрерывной правой части по **x**. Ниже — какие блоки ненулевые (в коде — `_jacobian_F_analytic`).

- **Строка 0 (dV/dt):**  
  ∂(dV/dt)/∂V, ∂(dV/dt)/∂α, ∂(dV/dt)/∂θ (через CD, тягу, гравитацию).

- **Строка 1 (dα/dt):**  
  ∂(dα/dt)/∂V, ∂(dα/dt)/∂α, ∂(dα/dt)/∂q, ∂(dα/dt)/∂θ, ∂(dα/dt)/∂CLα, ∂(dα/dt)/∂CLδe.

- **Строка 2 (dq/dt):**  
  ∂(dq/dt)/∂V, ∂(dq/dt)/∂α, ∂(dq/dt)/∂q, ∂(dq/dt)/∂Cmα, ∂(dq/dt)/∂Cmq̄, ∂(dq/dt)/∂Cmδe (через M = q̄ S c̄ Cm).

- **Строка 3 (dθ/dt):**  
  ∂(dθ/dt)/∂q = 1.

- **Строки 4–8 (параметры):**  
  В коде задаётся затухание: F_cont[i,i] = (0.999−1)/dt для i=4..8. Тогда в дискретной форме параметры слегка «затухают». В самой функции `f_dynamics` параметры не меняются — фактическая эволюция параметров в предсказании постоянна; несоответствие можно убрать, обнулив F_cont[4:9, 4:9].

Итог: **F** — матрица 9×9, дискретная: **F = I + F_cont*dt**.

**Примечание:** В аналитическом F для параметров (строки 4–8) задано затухание F_cont[i,i] = (0.999−1)/dt, что даёт в дискретной форме множитель ~0.999. В `f_dynamics` параметры при этом не меняются (возвращаются как есть). То есть эволюция состояния параметров в predict корректна (они постоянны), а ковариация P для параметров умножается на 0.999² на каждом шаге. Для строгой согласованности можно обнулить F_cont[4:9, 4:9], тогда F[4:9, 4:9] = I.

### 2.4 Вектор измерений и функция h (Parameter EKF)

Измерения **y** размерности **5**:

| Индекс | Величина | Описание |
|--------|----------|----------|
| 0 | V_m | Воздушная скорость (м/с) |
| 1 | α_m | Угол атаки с динамической ошибкой (рад) |
| 2 | q_m | Угловая скорость тангажа (рад/с) |
| 3 | θ_m | Угол тангажа (рад) |
| 4 | nz_m | Нормальная перегрузка (g) |

Модель измерений:

- V_m = V (в коде шум добавляется в датчике, в h — просто V).
- α_m = α + d_aoa · q/V (d_aoa — параметр датчика из конфига).
- q_m = q, θ_m = θ.
- nz = (L cos α + D sin α) / (mg); L, D считаются по текущему состоянию и параметрам.

В коде: `h_measurement` возвращает **[V, α_m, q, θ, nz]** по текущему **x** и **u**.

### 2.5 Якобиан H (Parameter EKF), 5×9

- Строка 0 (V_m): ∂V_m/∂V = 1, остальное 0.
- Строка 1 (α_m): ∂α_m/∂V = −d_aoa·q/V², ∂α_m/∂α = 1, ∂α_m/∂q = d_aoa/V.
- Строка 2 (q_m): ∂q_m/∂q = 1.
- Строка 3 (θ_m): ∂θ_m/∂θ = 1.
- Строка 4 (nz): производные по V, α, CLα, CLδe (через L, D и выражение для nz).

Точные формулы — в `_jacobian_H_analytic`. Итог: **H** — матрица **5×9**.

### 2.6 Ковариации Parameter EKF

- **initial_P**: диагональ (1, 1e-4, 1e-4, 1e-4, 0.1, 0.01, 0.01, 1, 0.01).
- **Q**: диагональ (0.5, 1e-6, 1e-6, 1e-6 для состояния; 1e-8, 1e-8, 1e-8, 1e-6, 1e-8 для параметров).
- **R**: диагональ (0.25, 1e-6, 1e-6, 1e-6, 0.000025) для [V, α_m, q, θ, nz].

После `update` оценки параметров ограничиваются: `x[4:9] = np.clip(x[4:9], param_bounds_lower, param_bounds_upper)`.

### 2.7 Предсказание и обновление (Parameter EKF)

- **predict(u):**  
  x := f_dynamics(x, u, dt)  
  F = jacobian_F(x_prev, u)  
  P := F P Fᵀ + Q  

- **update(y_meas, u):**  
  y_pred = h_measurement(x, u)  
  innovation = y_meas - y_pred  
  H = jacobian_H(x, u)  
  S = H P Hᵀ + R  
  K = P Hᵀ S⁻¹  
  x := x + K·innovation  
  P := (I − K H) P (I − K H)ᵀ + K R Kᵀ (форма Иосифа)  
  затем ограничение параметров и запись в param_history, measurement_history.

---

## 3. Навигация: БИНС

### 3.1 Модель БИНС (`models/ins_model.py`)

Внутреннее состояние объекта БИНС (не вектор EKF):

- **position** (3) — позиция в локальной НПЗ (NED), м.
- **velocity** (3) — скорость в NED, м/с.
- **quaternion** (4) — ориентация (тело → NED).
- **gyro_bias** (3), **accel_bias** (3) — смещения гироскопов и акселерометров (рад/с и м/с²).

### 3.2 Шаг распространения БИНС: `propagate(dt, specific_force_body, angular_rate_body, include_errors)`

Входы:

- **specific_force_body** — удельная сила в связанной СК (м/с²), без g.
- **angular_rate_body** — угловая скорость в связанной СК (рад/с).

При **include_errors=True**:

1. **Гироскопы:**  
   ω_meas = ω_true + gyro_bias + N(0, σ²_gyro_noise).  
   gyro_bias обновляется случайным блужданием: += N(0, (σ_bias·√dt·0.1)²).

2. **Акселерометры:**  
   f_meas = f_true + accel_bias + N(0, σ²_accel_noise).  
   accel_bias — случайное блуждание аналогично.

3. **Ориентация:**  
   Кватернион обновляется по ω_meas:  
   q̇ = (1/2) Ω(ω) q, дискретизация через кватернион приращения вращения (угол |ω|·dt, ось ω/|ω|), затем умножение кватернионов и нормализация.

4. **Скорость:**  
   DCM = quaternion_to_dcm(quaternion).  
   a_local = DCM · f_meas + g_local, g_local = [0, 0, −g].  
   velocity += a_local · dt.

5. **Позиция:**  
   position += velocity · dt.

Возвращается словарь: `position`, `velocity`, `attitude_euler`, `quaternion`, `time`.

Параметры σ задаются из `INSParameters` (класс точности): гироскопы (bias, ARW), акселерометры (bias, VRW).

---

## 4. Навигация: ГНСС

### 4.1 Модель ГНСС (`models/gnss_model.py`)

Приёмник выдаёт измерения с периодом **update_period = 1/update_rate** (сек).

**get_measurement(true_position, true_velocity, dt, available)**:

- Накапливает `time_since_update += dt`.
- Если `time_since_update < update_period`, возвращает последнее измерение (или None).
- Иначе сбрасывает таймер и, если `available`:
  - position_noise = N(0, position_std²), velocity_noise = N(0, velocity_std²);
  - возвращает:
    - **position** = true_position + position_noise,
    - **velocity** = true_velocity + velocity_noise,
    - а также std и флаги.

То есть модель: **идеальное измерение = истинная позиция/скорость + белый шум**. Класс точности задаёт `position_std` и `velocity_std` в `GNSSParameters`.

**Важно:** в текущем симуляторе в `_process_navigation` в `gnss.get_measurement(...)` передаются не истинные траектории в NED, а текущее состояние БИНС: `ins_state['position']`, `ins_state['velocity']`. То есть «истина» для ГНСС подставлена как БИНС; для строгой модели нужно было бы считать истинную траекторию в NED и передавать её в `get_measurement`.

---

## 5. Navigation EKF (комплексирование БИНС и ГНСС)

### 5.1 Вектор состояния (error-state)

Оценивается **вектор ошибок** размерности **15**:

| Индекс | Блок | Размер | Описание |
|--------|------|--------|----------|
| 0–2 | δpos | 3 | Ошибка позиции (м), NED |
| 3–5 | δvel | 3 | Ошибка скорости (м/с), NED |
| 6–8 | δatt | 3 | Ошибки ориентации (малые углы, рад) |
| 9–11 | gyro_bias | 3 | Смещение гироскопов (рад/с) |
| 12–14 | accel_bias | 3 | Смещение акселерометров (м/с²) |

Итого: **x = [δposᵀ, δvelᵀ, δattᵀ, gyro_biasᵀ, accel_biasᵀ]ᵀ**, **15×1**.

Смысл: истинное состояние ≈ состояние БИНС + оценка ошибки (для позиции/скорости/ориентации).

### 5.2 Динамика ошибок (Navigation EKF)

В коде `f_dynamics` для Navigation EKF просто возвращает **x** без изменения (интеграция ошибок выполняется в predict через матрицу F).

Линейная модель ошибок (по смыслу error-state):

- δpos_{k+1} = δpos_k + δvel_k · dt  
- δvel зависит от δatt и accel_bias (через удельную силу и DCM); в текущей реализации F упрощён.
- δatt_{k+1} ≈ δatt_k − gyro_bias · dt (вклад угловой скорости ошибки).
- gyro_bias и accel_bias — случайное блуждание (диагональ в F сохраняется 1).

### 5.3 Якобиан F (Navigation EKF), 15×15

В коде (`jacobian_F`):

- F = I₁₅.
- F[0:3, 3:6] = I·dt (позиция интегрирует скорость).
- F[6:9, 9:12] = −I·dt (ошибка ориентации под действием смещения гироскопов).

Остальные блоки (связь δvel с δatt и accel_bias, нелинейные члены) не заносятся — упрощённая модель.

### 5.4 Измерения и якобиан H (Navigation EKF)

Измерение ГНСС: **z = [position_GNSSᵀ, velocity_GNSSᵀ]ᵀ** (6×1).

Инновация в error-state формулировке:

**innovation = (position_GNSS − position_INS, velocity_GNSS − velocity_INS)**.

Интерпретация: измеряем «разницу ГНСС − БИНС», которая аппроксимирует ошибку БИНС по позиции и скорости.

Прогноз измерения в фильтре: **h(x) = x[0:6]** (оценки ошибок позиции и скорости). Тогда **H** — матрица 6×15: первые 6 столбцов — единичная матрица 6×6, остальные нули:

**H = [ I₆ | 0₆ₓ₉ ]**.

### 5.5 Обновление по ГНСС: `update_with_gnss(gnss_measurement, ins_state)`

- innovation = [gnss_position − ins_position, gnss_velocity − ins_velocity].
- y_measured = innovation.
- y_predicted = h(x) = x[0:6].
- actual_innovation = y_measured − y_predicted.
- Далее стандартный шаг EKF: H, S = H P Hᵀ + R, K = P Hᵀ S⁻¹, x := x + K·actual_innovation, P по форме Иосифа.

R строится из `gnss_params.position_std` и `velocity_std` (диагональ 6×6).

### 5.6 Скорректированное навигационное решение

**get_corrected_state(ins_state)** возвращает:

- position_corrected = ins_position **+** δpos  
- velocity_corrected = ins_velocity **+** δvel  
- attitude_corrected = ins_attitude **+** δatt  

То есть скорректированное состояние = БИНС + оценка ошибки (знак «+» соответствует определению ошибки как «истина − БИНС»).

### 5.7 Ковариации Navigation EKF

- **initial_P**: диагональ (1000 для позиции, 100 для скорости, 1 для ориентации, затем дисперсии по bias из ins_params).
- **Q**: диагональ, зависит от dt и параметров БИНС (усилена для учёта роста ошибок).
- **R**: diag(σ²_pos, σ²_pos, σ²_pos, σ²_vel, σ²_vel, σ²_vel) из gnss_params.

---

## 6. Базовый EKF (`filters/ekf_base.py`)

- **predict(u):**  
  x := f_dynamics(x, u, dt)  
  F = jacobian_F(x_prev, u)  
  P := F P Fᵀ + Q  
  (и симметризация, принудительная положительная определённость P).

- **update(y_measured, u):**  
  y_pred = h_measurement(x, u)  
  innovation = y_measured - y_pred  
  H = jacobian_H(x, u)  
  S = H P Hᵀ + R  
  K = P Hᵀ S⁻¹  
  x := x + K·innovation  
  P := (I − K H) P (I − K H)ᵀ + K R Kᵀ  
  плюс сохранение в историю (time, state, covariance, innovation).

Абстрактные методы: `f_dynamics`, `h_measurement`, `jacobian_F`, `jacobian_H`.

---

## 7. Датчики для идентификации (`sensors/sensor_models.py`)

**generate_sensor_data** возвращает вектор **[V_m, α_m, q_m, θ_m, nz_m]**:

- V_m: V_true + bias + N(0, σ²_V).
- α_m: α_true + d_aoa·q/V + bias + N(0, σ²_aoa).
- q_m, θ_m: истина + bias + шум.
- nz_m: значение nz по модели + bias + шум.

Эти зашумлённые измерения подаются в Parameter EKF как **y_meas** в **update**.

---

## 8. Краткая сводка размерностей и матриц

| Объект | Вектор/матрица | Размер |
|--------|-----------------|--------|
| Parameter EKF состояние | x | 9×1 |
| Parameter EKF измерение | y | 5×1 |
| Parameter EKF F | F | 9×9 |
| Parameter EKF H | H | 5×9 |
| Parameter EKF Q | Q | 9×9 (диагональ) |
| Parameter EKF R | R | 5×5 (диагональ) |
| Navigation EKF состояние | x (ошибки) | 15×1 |
| Navigation EKF измерение | z (инновация) | 6×1 |
| Navigation EKF F | F | 15×15 |
| Navigation EKF H | H | 6×15 |
| Navigation EKF Q | Q | 15×15 (диагональ) |
| Navigation EKF R | R | 6×6 (диагональ) |
| Состояние самолёта (истина) | [V, α, q, θ] | 4×1 |
| Состояние БИНС (выход) | position(3), velocity(3), attitude(3) | — |
| Измерение ГНСС | position(3), velocity(3) | 6×1 |

---

## 9. Файлы и ответственность

| Файл | Назначение |
|------|------------|
| `filters/ekf_base.py` | Базовый EKF: predict, update, история, адаптация |
| `filters/parameter_ekf.py` | Расширенное состояние 9, f_dynamics, h_measurement, F и H (аналитически/численно), ограничение параметров |
| `filters/navigation_ekf.py` | Error-state 15, F и H для ошибок, update_with_gnss, get_corrected_state |
| `models/aircraft_dynamics.py` | compute_derivatives, trim, load_factor |
| `models/aerodynamics.py` | L, D, M по V, α, q, δe и параметрам |
| `models/ins_model.py` | StrapdownINS: propagate (кватернион, DCM, скорость, позиция, ошибки датчиков) |
| `models/gnss_model.py` | GNSSReceiver: get_measurement (период, шум позиции/скорости) |
| `sensors/sensor_models.py` | Генерация зашумлённых [V_m, α_m, q_m, θ_m, nz_m] |
| `simulation/simulator.py` | Главный цикл: истинная динамика, БИНС, ГНСС, Navigation EKF, датчики, Parameter EKF, визуализация, history |
| `config/aircraft_params.py` | Параметры самолёта, TRUE_PARAMS, INITIAL_PARAMS, границы параметров |
| `config/navigation_params.py` | INSParameters, GNSSParameters, классы точности (commercial, tactical, navigation) |

Этот файл можно использовать как единый справочник по тому, что и как реализовано в коде: векторы состояния, матрицы якобианов, модели БИНС и ГНСС, предсказание и обновление в обоих фильтрах Калмана.
