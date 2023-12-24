import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded
from matplotlib.animation import FuncAnimation

# Параметры сетки и уравнения
L = 5.0  # длина отрезка
T = 1.0  # общее время
nx = 500  # число шагов по пространству
nt = 500  # число шагов по времени
dx = L / nx
dt = T / nt
a = 1.0

# Функция для начального условия
def initial_condition(x):
    return np.piecewise(x, [x < 1, (1 <= x) & (x < 4), x >= 4],
                        [0, lambda x: np.sin(np.pi * (x - 1) / 3), 0])

c = 1
def f(u):
    return c * u
    # return u ** 2 / 2 

# Инициализация сетки
x = np.linspace(-1, L, nx)
u0 = initial_condition(x)

# Явная схема "Крест"
def cross_scheme(u, a, dx, dt, nx):
    u_new = np.copy(u)
    for i in range(1, nx - 1):
        u_new[i] = u[i] - a * dt / (2 * dx) * (f(u[i+1]) - f(u[i-1]))
    return u_new

# Неявная схема с центральной разностью
def implicit_scheme(u, a, dx, dt):
    nx = len(u)
    A = np.zeros((3, nx))
    A[0, 1:] = -dt * a / (2 * dx)  # верхняя диагональ
    A[1, :] = 1                # главная диагональ
    A[2, :-1] = dt * a / (2 * dx)  # нижняя диагональ
    b = np.copy(f(u))
    return solve_banded((1, 1), A, b)

# Расчеты
u_cross = np.zeros((nt, nx))
u_implicit = np.zeros((nt, nx))
u_cross[0, :] = u0
u_implicit[0, :] = u0

for n in range(1, nt):
    u_cross[n, :] = cross_scheme(u_cross[n-1, :], a, dx, dt, nx)
    u_implicit[n, :] = implicit_scheme(u_implicit[n-1, :], a, dx, dt)

# Визуализация на последнем временном слое
# plt.figure(figsize=(12, 6))
# plt.plot(x, u_cross[-1, :], 'o-', label='Явная схема "Крест"')
# plt.plot(x, u_implicit[-1, :], 'x-', label='Неявная схема с центральной разностью')
# plt.title('Решение задачи переноса на последнем временном слое')
# plt.xlabel('x')
# plt.ylabel('u')
# plt.legend()
# plt.show()

fig, ax = plt.subplots()
line1, = ax.plot(x, u_cross[0, :], 'o-', label='Явная схема "Крест"')
line2, = ax.plot(x, u_implicit[0, :], 'x-', label='Неявная схема с центральной разностью')
ax.set_xlim(x[0], x[-1])
ax.set_ylim(-1.1, 1.1)
ax.set_xlabel('x')
ax.set_ylabel('u')
ax.set_title('Решение задачи переноса')
ax.legend()

# Функция анимации
def animate(i):
    line1.set_ydata(u_cross[i, :])
    line2.set_ydata(u_implicit[i, :])
    return line1, line2

# Создание анимации
ani = FuncAnimation(fig, animate, frames=nt, interval=20, blit=True)

# Показать анимацию
plt.show()