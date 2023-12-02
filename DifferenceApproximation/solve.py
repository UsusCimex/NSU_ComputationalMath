import numpy as np
import matplotlib.pyplot as plt

# Метод Эйлера (k=1)
def euler_method(f, x0, y0, h, xn):
    n = int((xn - x0) / h) + 1
    x = np.linspace(x0, xn, n)
    y = np.zeros(n)
    y[0] = y0
    for i in range(1, n):
        y[i] = y[i-1] + h * f(x[i-1], y[i-1])
    return x, y

# Улучшенный метод Эйлера или метод средней точки (k=2)
def improved_euler_method(f, x0, y0, h, xn):
    n = int((xn - x0) / h) + 1
    x = np.linspace(x0, xn, n)
    y = np.zeros(n)
    y[0] = y0
    for i in range(1, n):
        k1 = f(x[i-1], y[i-1])
        k2 = f(x[i-1] + h/2, y[i-1] + h/2 * k1)
        y[i] = y[i-1] + h * k2
    return x, y

# Метод Рунге-Кутты(k=4)
def runge_kutta_method(f, x0, y0, h, xn):
    n = int((xn - x0) / h) + 1
    x = np.linspace(x0, xn, n)
    y = np.zeros(n)
    y[0] = y0
    for i in range(1, n):
        k1 = f(x[i-1], y[i-1])
        k2 = f(x[i-1] + h/2, y[i-1] + h/2 * k1)
        k3 = f(x[i-1] + h/2, y[i-1] + h/2 * k2)
        k4 = f(x[i-1] + h, y[i-1] + h * k3)
        y[i] = y[i-1] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    return x, y


# Функция f(x) = y'(x) = -y(x)
# Аналитическое решение: y=\frac{1}{e^{x}}
def f(x, y):
    return -y

# Зададим начальные условия
y0 = 1
x0 = 0
xn = 5 # Предел интегрирования
h = 0.1 # Шаг интегрирования

x_euler, y_euler = euler_method(f, x0, y0, h, xn)
x_improved_euler, y_improved_euler = improved_euler_method(f, x0, y0, h, xn)
x_runge_kutta, y_runge_kutta = runge_kutta_method(f, x0, y0, h, xn)

# Построение графика
plt.figure(figsize=(12, 4))

plt.subplot(1, 3, 1)
plt.plot(x_euler, y_euler, label='Euler (k=1)', marker='o')
plt.title('Euler Method')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.legend()

plt.subplot(1, 3, 2)
plt.plot(x_improved_euler, y_improved_euler, label='Improved Euler (k=2)', marker='x')
plt.title('Improved Euler Method')
plt.xlabel('x')
plt.grid(True)
plt.legend()

plt.subplot(1, 3, 3)
plt.plot(x_runge_kutta, y_runge_kutta, label='Runge-Kutta (k=4)', marker='s')
plt.title('Runge-Kutta Method')
plt.xlabel('x')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()

# Функция g(x) = y^(k) = e^x * cos(x)
# Аналитическое решение: \frac{e^{x}\sin x}{2}+\frac{e^{x}\cos x}{2}-\frac{1}{2}
def g(x, y):
    return np.exp(x) * np.cos(x)

# Зададим начальные условия для g(x)
x0 = 0
xn = 3
h = 0.1
y0_g = 1  # g(0) = 1

# Применим методы к функции g(x)
x_euler, y_euler = euler_method(g, x0, y0_g, h, xn)
x_improved_euler, y_improved_euler = improved_euler_method(g, x0, y0_g, h, xn)
x_runge_kutta, y_runge_kutta = runge_kutta_method(g, x0, y0_g, h, xn)

# Построим графики
plt.figure(figsize=(12, 4))

plt.subplot(1, 3, 1)
plt.plot(x_euler, y_euler, label='Euler (k=1)', marker='o')
plt.title('Euler Method')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.legend()

plt.subplot(1, 3, 2)
plt.plot(x_improved_euler, y_improved_euler, label='Improved Euler (k=2)', marker='x')
plt.title('Improved Euler Method')
plt.xlabel('x')
plt.grid(True)
plt.legend()

plt.subplot(1, 3, 3)
plt.plot(x_runge_kutta, y_runge_kutta, label='Runge-Kutta (k=4)', marker='s')
plt.title('Runge-Kutta Method')
plt.xlabel('x')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()
