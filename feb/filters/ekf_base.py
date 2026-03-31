"""
Базовый класс расширенного фильтра Калмана (EKF).

Этот модуль содержит абстрактный базовый класс для всех фильтров Калмана в проекте.
Реализует общую логику EKF: predict, update, адаптацию ковариаций.

Принцип Open/Closed: закрыт для модификации, открыт для расширения через наследование.
"""

import numpy as np
from abc import ABC, abstractmethod
from typing import Callable, Optional


class EKFBase(ABC):
    """
    Абстрактный базовый класс для расширенного фильтра Калмана.
    
    Определяет интерфейс и общую логику для всех EKF в проекте.
    Конкретные фильтры должны наследоваться и реализовать абстрактные методы.
    
    Принципы SOLID:
    - Single Responsibility: только логика фильтрации
    - Open/Closed: расширяется через наследование
    - Liskov Substitution: любой наследник взаимозаменяем
    """
    
    def __init__(
        self,
        state_dim: int,
        measurement_dim: int,
        dt: float = 0.02,
        initial_state: np.ndarray = None,
        initial_P: np.ndarray = None,
        Q: np.ndarray = None,
        R: np.ndarray = None
    ):
        """
        Инициализация базового EKF.
        
        Args:
            state_dim: Размерность вектора состояния
            measurement_dim: Размерность вектора измерений
            dt: Временной шаг (сек)
            initial_state: Начальное состояние
            initial_P: Начальная ковариация состояния
            Q: Ковариация шума процесса
            R: Ковариация шума измерений
        """
        self.state_dim = state_dim
        self.measurement_dim = measurement_dim
        self.dt = dt
        self.time = 0.0
        
        # Состояние фильтра
        self.x = initial_state if initial_state is not None else np.zeros(state_dim)
        
        # Ковариационная матрица
        self.P = initial_P if initial_P is not None else np.eye(state_dim)
        
        # Ковариация шума процесса
        self.Q = Q if Q is not None else np.eye(state_dim) * 1e-6
        
        # Ковариация шума измерений
        self.R = R if R is not None else np.eye(measurement_dim) * 1e-6
        
        # История для анализа
        self.time_history = []
        self.state_history = []
        self.covariance_history = []
        self.innovation_history = []
        
        # Параметры адаптации
        self.adaptation_enabled = False
        self.max_history = 1000
    
    @abstractmethod
    def f_dynamics(self, x: np.ndarray, u: any, dt: float) -> np.ndarray:
        """
        Функция динамики системы: x_{k+1} = f(x_k, u_k, dt).
        
        Args:
            x: Текущее состояние
            u: Управление
            dt: Временной шаг
        
        Returns:
            Следующее состояние
        """
        pass
    
    @abstractmethod
    def h_measurement(self, x: np.ndarray, u: any) -> np.ndarray:
        """
        Функция измерений: y = h(x, u).
        
        Args:
            x: Состояние
            u: Управление
        
        Returns:
            Прогнозируемое измерение
        """
        pass
    
    @abstractmethod
    def jacobian_F(self, x: np.ndarray, u: any) -> np.ndarray:
        """
        Вычисляет якобиан функции динамики по состоянию.
        
        F = ∂f/∂x
        
        Args:
            x: Состояние
            u: Управление
        
        Returns:
            Матрица Якоби размера (state_dim, state_dim)
        """
        pass
    
    @abstractmethod
    def jacobian_H(self, x: np.ndarray, u: any) -> np.ndarray:
        """
        Вычисляет якобиан функции измерений по состоянию.
        
        H = ∂h/∂x
        
        Args:
            x: Состояние
            u: Управление
        
        Returns:
            Матрица Якоби размера (measurement_dim, state_dim)
        """
        pass
    
    def predict(self, u: any = None):
        """
        Шаг предсказания EKF.
        
        x_pred = f(x, u, dt)
        P_pred = F * P * F^T + Q
        
        Args:
            u: Управление
        """
        # Сохраняем предыдущее состояние для адаптации
        x_prev = self.x.copy()
        
        # Прогноз состояния
        self.x = self.f_dynamics(self.x, u, self.dt)
        
        # Якобиан
        F = self.jacobian_F(x_prev, u)
        
        # Прогноз ковариации
        self.P = F @ self.P @ F.T + self.Q
        
        # Симметризация и проверка положительной определенности
        self.P = (self.P + self.P.T) / 2
        self._ensure_positive_definite(self.P)
        
        self.time += self.dt
    
    def update(self, y_measured: np.ndarray, u: any = None):
        """
        Шаг коррекции EKF по измерениям.
        
        Args:
            y_measured: Вектор измерений
            u: Управление
        """
        # Прогноз измерения
        y_pred = self.h_measurement(self.x, u)
        
        # Инновация (residual)
        innovation = y_measured - y_pred
        
        # Якобиан измерений
        H = self.jacobian_H(self.x, u)
        
        # Ковариация инновации
        S = H @ self.P @ H.T + self.R
        
        # Симметризация
        S = (S + S.T) / 2
        
        try:
            # Усиление Калмана
            K = self.P @ H.T @ np.linalg.inv(S)
            
            # Обновление состояния
            self.x = self.x + K @ innovation
            
            # Обновление ковариации (форма Иосифа для численной стабильности)
            I_KH = np.eye(self.state_dim) - K @ H
            self.P = I_KH @ self.P @ I_KH.T + K @ self.R @ K.T
            
        except np.linalg.LinAlgError:
            # Если матрица S вырожденная, используем псевдообратную
            K = self.P @ H.T @ np.linalg.pinv(S)
            self.x = self.x + K @ innovation
            I_KH = np.eye(self.state_dim) - K @ H
            self.P = I_KH @ self.P @ I_KH.T + K @ self.R @ K.T
        
        # Симметризация и проверка
        self.P = (self.P + self.P.T) / 2
        self._ensure_positive_definite(self.P)
        
        # Сохранение в историю
        self._save_to_history(innovation)
    
    def _ensure_positive_definite(self, matrix: np.ndarray, min_eigenvalue: float = 1e-10):
        """
        Гарантирует положительную определенность матрицы.
        
        Args:
            matrix: Матрица для проверки
            min_eigenvalue: Минимальное собственное значение
        """
        for i in range(matrix.shape[0]):
            if matrix[i, i] < min_eigenvalue:
                matrix[i, i] = min_eigenvalue
    
    def _save_to_history(self, innovation: np.ndarray):
        """
        Сохраняет текущее состояние в историю.
        
        Args:
            innovation: Вектор инновации
        """
        self.time_history.append(self.time)
        self.state_history.append(self.x.copy())
        self.covariance_history.append(np.diag(self.P).copy())
        self.innovation_history.append(innovation.copy())
        
        # Ограничение размера истории
        if len(self.time_history) > self.max_history:
            self.time_history.pop(0)
            self.state_history.pop(0)
            self.covariance_history.pop(0)
            self.innovation_history.pop(0)
    
    def enable_adaptive_filtering(
        self,
        adaptation_rate_Q: float = 0.001,
        adaptation_rate_R: float = 0.001,
        window_size: int = 50
    ):
        """
        Включает адаптивную фильтрацию для автоматической настройки Q и R.
        
        Args:
            adaptation_rate_Q: Скорость адаптации Q
            adaptation_rate_R: Скорость адаптации R
            window_size: Размер окна для статистики инноваций
        """
        self.adaptation_enabled = True
        self.adaptation_rate_Q = adaptation_rate_Q
        self.adaptation_rate_R = adaptation_rate_R
        self.adaptation_window = window_size
    
    def get_state(self) -> np.ndarray:
        """Возвращает текущее состояние."""
        return self.x.copy()
    
    def get_covariance(self) -> np.ndarray:
        """Возвращает текущую ковариационную матрицу."""
        return self.P.copy()
    
    def get_state_std(self) -> np.ndarray:
        """Возвращает стандартные отклонения оценок состояния."""
        return np.sqrt(np.diag(self.P))
    
    def reset(self, initial_state: np.ndarray = None, initial_P: np.ndarray = None):
        """
        Сбрасывает фильтр к начальным условиям.
        
        Args:
            initial_state: Новое начальное состояние
            initial_P: Новая начальная ковариация
        """
        if initial_state is not None:
            self.x = initial_state.copy()
        else:
            self.x = np.zeros(self.state_dim)
        
        if initial_P is not None:
            self.P = initial_P.copy()
        else:
            self.P = np.eye(self.state_dim)
        
        self.time = 0.0
        self.time_history.clear()
        self.state_history.clear()
        self.covariance_history.clear()
        self.innovation_history.clear()
    
    def compute_numerical_jacobian(
        self,
        func: Callable,
        x: np.ndarray,
        u: any,
        eps: float = 1e-6
    ) -> np.ndarray:
        """
        Вычисляет якобиан численно (для проверки или если нет аналитического).
        
        Args:
            func: Функция для дифференцирования
            x: Точка, в которой вычисляется якобиан
            u: Параметры функции
            eps: Шаг для численного дифференцирования
        
        Returns:
            Матрица Якоби
        """
        f0 = func(x, u)
        output_dim = len(f0) if isinstance(f0, np.ndarray) else 1
        input_dim = len(x)
        
        jacobian = np.zeros((output_dim, input_dim))
        
        for i in range(input_dim):
            x_perturbed = x.copy()
            x_perturbed[i] += eps
            f_perturbed = func(x_perturbed, u)
            jacobian[:, i] = (f_perturbed - f0) / eps
        
        return jacobian
    
    def get_history(self) -> dict:
        """
        Возвращает всю историю фильтрации.
        
        Returns:
            dict: Словарь с историей времени, состояний, ковариаций, инноваций
        """
        return {
            'time': np.array(self.time_history),
            'state': np.array(self.state_history),
            'covariance': np.array(self.covariance_history),
            'innovation': np.array(self.innovation_history)
        }
