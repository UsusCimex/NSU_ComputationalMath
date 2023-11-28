from sympy import symbols, log, diff, lambdify, exp, sin, printing
from math import pi
import numpy as np
import matplotlib.pyplot as plt

NUM_ITERATIONS = 10

def simple_iteration(func, a_val, x0, iterations=NUM_ITERATIONS):
    x_vals = [x0]
    for i in range(iterations):
        x_new = func(x_vals[-1], a_val)
        x_vals.append(x_new)
    return x_vals

def newton_method(func, func_prime, a_val, x0, iterations=NUM_ITERATIONS):
    x_vals = [x0]
    for i in range(iterations):
        f_x = func(x_vals[-1], a_val)
        f_prime_x = func_prime(x_vals[-1], a_val)
        x_new = x_vals[-1] - f_x/f_prime_x
        x_vals.append(x_new)
    return x_vals

# Пользователь определяет свою функцию и её производную
x = symbols('x')
a = symbols('a')

# x - a = ln(1 + x)
sym_function_simple = a + log(1 + x)

sym_function_newton = a + log(1 + x) - x
sym_function_newton_prime = diff(sym_function_newton, x)  # Производная функции

# x = e^(x - a) - 1
# sym_function_simple = exp(x - a) - 1

# sym_function_newton = exp(x - a) - 1 - x
# sym_function_newton_prime = diff(sym_function_newton, x)  # Производная функции

# x = a sin(x)
# sym_function_simple = a * sin(x);

# sym_function_newton = a * sin(x) - x
# sym_function_newton_prime = diff(sym_function_newton, x)  # Производная функции

# Превращение символьных функций в численные функции
func_simple = lambdify([x, a], sym_function_simple)

func_newton = lambdify([x, a], sym_function_newton)
func_newton_prime = lambdify([x, a], sym_function_newton_prime)

a_val = 5
x0 = 0.1

simple_iter_values = simple_iteration(func_simple, a_val, x0)
print("simple: " + str(simple_iter_values[-1]))

# for i in range(0, 10, 1):
#     x0 = pi / 4 + (pi * i)
#     newtons_values = newton_method(func_newton, func_newton_prime, a_val, x0)
#     if (abs(newtons_values[-1] - newtons_values[-2]) > 0.1):
#         break
#     print(f"x={round(x0, 2)}: newton: {newtons_values[-1]}")
#     print(f"x={round(-x0, 2)}: newton: {-newtons_values[-1]}")

newtons_values = newton_method(func_newton, func_newton_prime, a_val, x0)
print("newton: " + str(newtons_values[-1]))

func_simple_str = printing.latex(sym_function_simple)
func_newton_str = printing.latex(sym_function_newton)

fig, ax = plt.subplots()
plt.rc('text', usetex=True)
ax.plot(simple_iter_values, label='Simple Iteration', marker='o')
ax.plot(newtons_values, label='Newton\'s Method', marker='x')
ax.set_xlabel('Iteration')
ax.set_ylabel('Value of $x$')
ax.set_title(f'Convergence of Iteration Methods for $x = {func_simple_str}$')
ax.legend()

iteration_convergence_variant = 'iteration_convergence_variant.png'
plt.savefig(iteration_convergence_variant)
plt.close()

x_vals = np.linspace(-10, 10, 400)
f_x_vals = x_vals
g_x_vals = func_simple(x_vals, a_val)

plt.figure(figsize=(10, 6))
plt.plot(x_vals, f_x_vals, label='$f(x) = x$')
plt.plot(x_vals, g_x_vals, label=f'$g(x) = {func_simple_str}$', linestyle='--')
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'Graphs of $f(x) = x$ and $g(x) = {func_simple_str}$')
plt.legend()
plt.grid(True)

for i in range(1, len(simple_iter_values)):
    prev_x, new_x = simple_iter_values[i-1], simple_iter_values[i]
    plt.plot([prev_x, prev_x], [prev_x, new_x], 'r:', alpha=0.5)
    plt.plot([prev_x, new_x], [new_x, new_x], 'r:', alpha=0.5)
    plt.plot(new_x, new_x, 'ro', alpha=0.5)

##for i in range(1, len(newtons_values)):
##    prev_x, new_x = newtons_values[i-1], newtons_values[i]
##    plt.plot([prev_x, prev_x], [prev_x, new_x], 'b:', alpha=0.5)
##    plt.plot([prev_x, new_x], [new_x, new_x], 'b:', alpha=0.5)
##    plt.plot(new_x, new_x, 'bo', alpha=0.5)

iterative_process_plot_v2_path = 'iterative_process_plot_v2.png'
plt.savefig(iterative_process_plot_v2_path)
plt.show()
plt.close()
