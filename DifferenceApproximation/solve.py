import numpy as np
import matplotlib.pyplot as plt
import math

def y_first_order_approximation(x0, y0, h, xn):
    n = int((xn - x0) / h) + 1
    x = np.linspace(x0, xn, n)
    y = np.zeros(n)
    y[0] = y0
    for i in range(n - 1):
        y[i+1] = y[i]*(1 - h)
    return x,y

def y_second_order_approximation(x0, y0, h, xn):
    n = int((xn - x0) / h) + 1
    x = np.linspace(x0, xn, n)
    y = np.zeros(n)
    y[0] = y0
    for i in range(1, n):
        y[i] = y[i-1] * (1 - h/2) / (1 + h/2)
    return x, y

def y_fourth_order_approximation(x0, y0, h, xn):
    n = int((xn - x0) / h) + 1
    x = np.linspace(x0, xn, n)
    y = np.zeros(n)
    y[0] = y0
    y[1] = y[0] * (2 - h) / (2 + h)
    for i in range(1, n - 1):
        y[i+1] = (y[i-1] - h/3*(y[i-1] + 4*y[i])) / (1 + h/3)
    return x, y

def first_order_approximation(f, x0, y0, h, xn):
    n = int((xn - x0) / h) + 1
    x = np.linspace(x0, xn, n)
    y = np.zeros(n)
    y[0] = y0
    for i in range(1, n):
        y[i] = y[i-1] + h * f(x[i-1])
    return x, y

def second_order_approximation(f, x0, y0, h, xn):
    n = int((xn - x0) / h) + 1
    x = np.linspace(x0, xn, n)
    y = np.zeros(n)
    y[0] = y0
    for i in range(1, n):
        y[i] = y[i-1] + (h/2) * (f(x[i-1]) + f(x[i]))
    return x, y

def fourth_order_approximation(f, x0, y0, h, xn): 
    n = int((xn - x0) / h) + 1
    x = np.linspace(x0, xn, n)
    y = np.zeros(n)
    y[0] = y0
    y[1] = y[0] + h * (f(x[0]) + f(x[1])) / 2
    for i in range(1, n - 1):
        y[i+1] = y[i-1] + h/3 * (f(x[i-1]) + 4*f(x[i]) + f(x[i+1]))
    return x, y

# Функция f(x) = y'(x) = -y(x)
# Аналитическое решение: y=\frac{1}{e^{x}}

# Зададим начальные условия
y0 = 1
x0 = 0
xn = 5 # Предел интегрирования
h = 0.1 # Шаг интегрирования

x_1, y_1 = y_first_order_approximation(x0, y0, h, xn)
x_2, y_2 = y_second_order_approximation(x0, y0, h, xn)
x_4, y_4 = y_fourth_order_approximation(x0, y0, h, xn)

# Построение графика
plt.figure(figsize=(12, 4))

plt.subplot(1, 3, 1)
plt.plot(x_1, y_1, label='first_order_approximation', marker='o')
plt.title('first_order_approximation')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.legend()

plt.subplot(1, 3, 2)
plt.plot(x_2, y_2, label='second_order_approximation', marker='x')
plt.title('second_order_approximation')
plt.xlabel('x')
plt.grid(True)
plt.legend()

plt.subplot(1, 3, 3)
plt.plot(x_4, y_4, label='fourth_order_approximation', marker='o')
plt.title('fourth_order_approximation')
plt.xlabel('x')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()

# Функция g(x) = y^(k) = e^x * cos(x)
# Аналитическое решение: \frac{e^{x}\sin x}{2}+\frac{e^{x}\cos x}{2}-\frac{1}{2}
def g(x):
    return np.exp(x) * np.cos(x)

# Зададим начальные условия для g(x)
x0 = 0
xn = 5
h = 0.1
y0_g = 1  # g(0) = 1

# Правило Рунге
# p = log_2{(I2h - I) / (Ih - I)}

x4h, y4h = fourth_order_approximation(g, x0, y0_g, (xn - x0) / 100, xn)
x2h, y2h = fourth_order_approximation(g, x0, y0_g, (xn - x0) / 40, xn)
xh, yh =   fourth_order_approximation(g, x0, y0_g, (xn - x0) / 20, xn)
p = np.zeros(len(xh))

for i in range(len(xh)):
    p[i] = np.log2(abs( (yh[i] - y4h[5*i]) / (y2h[2*i] - y4h[5*i]) ))

print(p)

# Применим методы к функции g(x)
x_1, y_1 = first_order_approximation(g, x0, y0_g, h, xn)
x_2, y_2 = second_order_approximation(g, x0, y0_g, h, xn)
x_4, y_4 = fourth_order_approximation(g, x0, y0_g, h, xn)

# Построим графики
plt.figure(figsize=(12, 4))

plt.subplot(1, 3, 1)
plt.plot(x_1, y_1, label='first_order_approximation', marker='o')
plt.title('first_order_approximation')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.legend()

plt.subplot(1, 3, 2)
plt.plot(x_2, y_2, label='second_order_approximation', marker='x')
plt.title('second_order_approximation')
plt.xlabel('x')
plt.grid(True)
plt.legend()

plt.subplot(1, 3, 3)
plt.plot(x_4, y_4, label='fourth_order_approximation', marker='s')
plt.title('fourth_order_approximation')
plt.xlabel('x')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()
