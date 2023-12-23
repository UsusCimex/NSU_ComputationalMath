import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

# Параметры сетки и уравнения
L = 5.0  # длина отрезка
T = 1.0  # общее время
nx = 100  # число шагов по пространству
nt = 100  # число шагов по времени
dx = L / nx
dt = T / nt
c = 1.0   # скорость волны

# Функция для начального условия
def initial_condition(x):
    return np.piecewise(x, [x < 1, (1 <= x) & (x < 4), x >= 4],
                        [0, lambda x: np.sin(np.pi * (x - 1) / 3), 0])

# Инициализация сетки
x = np.linspace(-1, L, nx)
u0 = initial_condition(x)

# Явная схема "Крест"
def cross_scheme(u, c, dx, dt, nx):
    u_new = np.copy(u)
    for i in range(1, nx - 1):
        u_new[i] = u[i] - c * dt / (2 * dx) * (u[i+1] - u[i-1])
    return u_new

# Неявная схема с центральной разностью
def implicit_scheme(u, dx, dt):
    nx = len(u)
    A = np.zeros((3, nx))
    A[0, 1:] = -dt / (2 * dx)  # верхняя диагональ
    A[1, :] = 1                # главная диагональ
    A[2, :-1] = dt / (2 * dx)  # нижняя диагональ
    b = np.copy(u)
    return solve_banded((1, 1), A, b)

# Расчеты
u_cross = np.zeros((nt, nx))
u_implicit = np.zeros((nt, nx))
u_cross[0, :] = u0
u_implicit[0, :] = u0

for n in range(1, nt):
    u_cross[n, :] = cross_scheme(u_cross[n-1, :], c, dx, dt, nx)
    u_implicit[n, :] = implicit_scheme(u_implicit[n-1, :], dx, dt)

# Визуализация
plt.figure(figsize=(12, 6))
plt.plot(x, u_cross[-1, :], 'o-', label='Явная схема "Крест"')
plt.plot(x, u_implicit[-1, :], 'x-', label='Неявная схема с центральной разностью')
plt.title('Решение задачи переноса на последнем временном слое')
plt.xlabel('x')
plt.ylabel('u')
plt.legend()
plt.show()
