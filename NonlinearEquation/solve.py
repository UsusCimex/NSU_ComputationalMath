from sympy import symbols, Eq, solve, log
from sympy.plotting import plot
import matplotlib.pyplot as plt
import numpy as np

def simple_iteration(a_val, x0, iterations=10):
    x_vals = [x0]
    for i in range(iterations):
        x_new = a_val + log(1 + x_vals[-1])
        x_vals.append(x_new)
    return x_vals

def newtons_method(a_val, x0, iterations=10):
    x_vals = [x0]
    for i in range(iterations):
        f_x = x_vals[-1] - a_val - log(1 + x_vals[-1])
        f_prime_x = 1 - 1 / (1 + x_vals[-1])
        x_new = x_vals[-1] - f_x/f_prime_x
        x_vals.append(x_new)
    return x_vals

a_val = 5
x0 = 0.5

simple_iter_values = simple_iteration(a_val, x0)
newtons_values = newtons_method(a_val, x0)

print("simple: " + str(simple_iter_values[-1]))
print("newton: " + str(newtons_values[-1]))

fig, ax = plt.subplots()
ax.plot(simple_iter_values, label='Simple Iteration', marker='o')
ax.plot(newtons_values, label='Newton\'s Method', marker='x')
ax.set_xlabel('Iteration')
ax.set_ylabel('Value of x')
ax.set_title('Convergence of Iteration Methods for x = l*sin(x)')
ax.legend()

iteration_convergence_variant = 'iteration_convergence_variant.png'
plt.savefig(iteration_convergence_variant)
plt.close()

x_vals = np.linspace(0, 10, 400)
f_x_vals = x_vals
g_x_vals = a_val + np.log1p(x_vals)

plt.figure(figsize=(10, 6))
plt.plot(x_vals, f_x_vals, label='f(x) = x')
plt.plot(x_vals, g_x_vals, label='g(x) = a + ln(1+x)', linestyle='--')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Graphs of f(x) = x and g(x) = a + ln(1+x)')
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
