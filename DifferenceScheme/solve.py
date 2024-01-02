import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded
from matplotlib.animation import FuncAnimation

left = 2.0
right = 1.0
# Функция для начального условия
def initial_condition(x):
    # return np.piecewise(x, [x < 1,         (1 <= x) & (x < 4)           , x >= 4],
    #                        [0    , lambda x: np.sin(np.pi * (x - 1) / 3), 0     ])
    return np.piecewise(x, [x < 0, x >= 0],
                           [left , right ])

c = 1
def f(u):
    # return c * u
    return (u ** 2) / 2 

# Инициализация сетки
# L = 5.0  # длина отрезка
# T = 1.0  # общее время
# nx = 500  # число шагов по пространству
# nt = 100  # число шагов по времени
# dx = L / nx # h
# dt = T / nt # tau
# x = np.linspace(0, 5, nx)

L = 2.0  # длина отрезка
T = 1.0  # общее время
nx = 200  # число шагов по пространству
nt = 100  # число шагов по времени
dx = L / nx # h
dt = T / nt # tau
x = np.linspace(-1, 1, nx)

# Явная схема "Крест"
def cross_scheme(u_old, u, dx, dt, nx):
    u_new = np.copy(u)
    
    dt = dx / max(u)
    # a = c
    a = 1
    for i in range(nx):
        if (i == 0):
            # u_new[i] = u_old[i] - a * dt / dx * (f(u[i+1]) - f(0))
            u_new[i] = u_old[i] - a * dt / dx * (f(u[i+1]) - f(left))
        elif (i == nx - 1):
            # u_new[i] = u_old[i] - a * dt / dx * (f(0) - f(u[i-1]))
            u_new[i] = u_old[i] - a * dt / dx * (f(right) - f(u[i-1]))
        else:
            u_new[i] = u_old[i] - a * dt / dx * (f(u[i+1]) - f(u[i-1]))
    print(f"Courant number: {a * dt / dx}")
    return u_new

# Явная схема Годунова
def godunov_scheme(u, dx, dt, nx):
    u_new = np.copy(u)
    dt = dx / max(u)
    for i in range(nx):
        # a = c
        # a = u[i]
        a = 1
        if i == 0:
            # u_new[i] = u[i] - a * dt / dx * f(u[i])
            u_new[i] = u[i] + a * dt / dx * (f(left) - f(u[i]))
        else:
            u_new[i] = u[i] + a * dt / dx * (f(u[i-1]) - f(u[i]))
    print(f"Courant number: {a * dt / dx}")
    return u_new

def persecution(A, b):
    n = len(b)
    c = np.zeros(n-1)
    d = np.zeros(n)
    
    # Прямая прогонка
    c[0] = A[0, 1] / A[1, 0]
    d[0] = b[0] / A[1, 0]
    for i in range(1, n-1):
        temp = A[1, i] - A[2, i-1] * c[i-1]
        c[i] = A[0, i+1] / temp
        d[i] = (b[i] - A[2, i-1] * d[i-1]) / temp

    d[n-1] = (b[n-1] - A[2, n-2] * d[n-2]) / (A[1, n-1] - A[2, n-2] * c[n-2])

    # Обратная прогонка
    x = np.zeros(n)
    x[-1] = d[-1]
    for i in range(n-2, -1, -1):
        x[i] = d[i] - c[i] * x[i+1]
    return x

def implicit_scheme(u, dx, dt):
    nx = len(u)
    A = np.zeros((3, nx))
    dt = dx / max(u) / 2
    # a = c
    for i in range(nx):
        a = u[i]
        A[2, i] = -dt * a / (2 * dx)  # Нижняя диагональ
        A[1, i] = 1                   # Главная диагональ
        A[0, i] = dt * a / (2 * dx)   # Верхняя диагональ
        print(f"Courant number: {a * dt / dx}")
    nb = np.copy(u)
    return persecution(A, nb)

# Расчеты
u_cross = np.zeros((nt, nx))
u_implicit = np.zeros((nt, nx))
u0 = initial_condition(x)
u_cross[0, :] = u0
u_implicit[0, :] = u0

for n in range(1, nt):
    print(f"timeline: {n}")
    print("Cross scheme")
    if (n == 1):
        u_cross[n, :] = godunov_scheme(u_cross[n-1, :], dx, dt, nx)
    else:
        u_cross[n, :] = cross_scheme(u_cross[n-2, :], u_cross[n-1, :], dx, dt, nx)
    print("Implicit scheme")
    u_implicit[n, :] = implicit_scheme(u_implicit[n-1, :], dx, dt)
    u_implicit[n, :5] = u_implicit[n-1, :5]
    u_implicit[n, -5:] = u_implicit[n-1, -5:]

# print(u_implicit[0, :])
# print(u_implicit[1, :])
# print(u_implicit[2, :])
# print(u_implicit[3, :])

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
ax.set_ylim(0, 4)
ax.set_xlabel('x')
ax.set_ylabel('u')
ax.set_title('Решение задачи переноса')
ax.legend()

def animate(i):
    line1.set_ydata(u_cross[i, :])
    line2.set_ydata(u_implicit[i, :])
    return line1, line2
ani = FuncAnimation(fig, animate, frames=nt, interval=50, blit=True)

plt.show()