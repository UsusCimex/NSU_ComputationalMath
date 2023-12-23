# Исправим ошибку и сформируем правильно матрицу для solve_banded
def implicit_scheme_corrected(u, dx, dt):
    nx = len(u)
    # Матрица системы A и вектор b
    A = np.zeros((3, nx))  # Трехдиагональная матрица для solve_banded
    b = np.zeros(nx)
    
    # Заполнение диагоналей матрицы A
    A[1, :] = 1  # Главная диагональ
    A[0, 1:] = -dt / (2 * dx)  # Верхняя диагональ
    A[2, :-1] = -dt / (2 * dx)  # Нижняя диагональ

    # Вектор правой части b
    b[1:-1] = u[1:-1] - dt / (2 * dx) * (0.5 * u[2:]**2 - 0.5 * u[:-2]**2)

    # Граничные условия
    b[0] = u[0]
    b[-1] = u[-1]

    # Решаем систему линейных уравнений с помощью solve_banded
    u_new = solve_banded((1, 1), A, b)

    return u_new

# Выполнение расчетов для всех временных слоев
u_cross_solutions = np.zeros((nt, nx))
u_implicit_solutions = np.zeros((nt, nx))

# Задаем начальные условия для обеих схем
u_cross_solutions[0, :] = u0
u_implicit_solutions[0, :] = u0

# Основной цикл по времени для решения уравнений
for t in range(1, nt):
    # Явная схема "Крест"
    u_cross_solutions[t, :] = explicit_cross_scheme(u_cross_solutions[t - 1, :], c, dx, dt, nx)
    # Неявная схема с центральной разностью
    u_implicit_solutions[t, :] = implicit_scheme_corrected(u_implicit_solutions[t - 1, :], dx, dt)

# Визуализация результатов для последнего временного слоя
plt.figure(figsize=(12, 6))

plt.plot(x, u_cross_solutions[-1, :], label='Явная схема "Крест"', marker='o')
plt.plot(x, u_implicit_solutions[-1, :], label='Неявная схема с центральной разностью', marker='x')

plt.title('Сравнение явной и неявной схем на последнем временном слое')
plt.xlabel('x')
plt.ylabel('u')
plt.legend()
plt.grid(True)
plt.show()
