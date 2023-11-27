import numpy as np

def find_eigenvalues(matrix):
    eigenvalues, _ = np.linalg.eig(matrix)
    return eigenvalues

print(" Матрица Гилиберта")
matrix = np.array([[1, 1/2, 1/3], # Матрица Гильберта3 ~ 524, для 4 ~ 15513
                  [1/2, 1/3, 1/4],
                  [1/3, 1/4, 1/5]])

eigenvalues = find_eigenvalues(matrix)
print("Собственные числа матрицы:")
print(eigenvalues)
print("Норма матрицы A2")
print(abs(max(eigenvalues)) / abs(min(eigenvalues)))