#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

void LU_Decomposition(vector<vector<double>>& A, vector<vector<double>>& L, vector<vector<double>>& U) {
    int n = A.size();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i < j) {
                L[i][j] = 0;
                U[j][i] = 0;
            }
            else if (j == i) {
                L[j][i] = 1;
                U[j][i] = A[j][i];
            }
            else {
                L[j][i] = A[j][i] / U[i][i];
                U[j][i] = 0;
                for (int k = 0; k < i; ++k) {
                    U[j][i] += -L[j][k] * U[k][i];
                }
            }
        }
    }
}

vector<double> solveUsingLU(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));
    
    LU_Decomposition(A, L, U);
    
    // Решение системы Ly = b
    vector<double> y(n, 0);
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < i; ++j) {
            sum += L[i][j] * y[j];
        }
        y[i] = (b[i] - sum) / L[i][i];
    }

    // Решение системы Ux = y
    vector<double> x(n, 0);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; ++j) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }

    return x;
}

void QR_Decomposition(vector<vector<double>>& A, vector<vector<double>>& Q, vector<vector<double>>& R) {
    int n = A.size();

    for (int j = 0; j < n; ++j) {
        vector<double> v(n, 0);
        for (int i = 0; i < j; ++i) {
            double dot_product = 0;
            for (int k = 0; k < n; ++k) {
                dot_product += Q[k][i] * A[k][j];
            }
            for (int k = 0; k < n; ++k) {
                v[k] += dot_product * Q[k][i];
            }
        }

        double norm_v = 0;
        for (int k = 0; k < n; ++k) {
            v[k] = A[k][j] - v[k];
            norm_v += v[k] * v[k];
        }

        norm_v = sqrt(norm_v);

        for (int k = 0; k < n; ++k) {
            Q[k][j] = v[k] / norm_v;
            R[j][k] = 0;
            for (int i = 0; i < n; ++i) {
                R[j][k] += Q[i][j] * A[i][k];
            }
        }
    }
}

// Решение системы методом QR-разложения
vector<double> solveUsingQR(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    vector<vector<double>> Q(n, vector<double>(n, 0));
    vector<vector<double>> R(n, vector<double>(n, 0));
    
    QR_Decomposition(A, Q, R);

    // Решение системы Q^T * Q * x = Q^T * b
    vector<double> Qt_b(n, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Qt_b[i] += Q[j][i] * b[j];
        }
    }

    // Решение системы R * x = Qt_b методом обратной подстановки
    vector<double> x(n, 0);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; ++j) {
            sum += R[i][j] * x[j];
        }
        x[i] = (Qt_b[i] - sum) / R[i][i];
    }

    return x;
}

// Функция для решения системы методом Гаусса-Зейделя
vector<double> solveUsingGaussSeidel(const vector<vector<double>> A, const vector<double> b, double tol = 1e-6, int max_iter = 100) {
    int n = A.size();
    vector<double> x(n, 0);
    vector<double> x_new(n, 0);
    
    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < n; ++i) {
            double sum1 = 0;
            double sum2 = 0;
            for (int j = 0; j < i; ++j) {
                sum1 += A[i][j] * x_new[j];
            }
            for (int j = i + 1; j < n; ++j) {
                sum2 += A[i][j] * x[j];
            }
            x_new[i] = (b[i] - sum1 - sum2) / A[i][i];
        }
        
        // Проверка сходимости
        bool converged = true;
        for (int i = 0; i < n; ++i) {
            if (fabs(x_new[i] - x[i]) > tol) {
                converged = false;
                break;
            }
        }
        
        if (converged) {
            return x_new;
        }
        
        x = x_new;
    }
    
    return x_new;
}

// Функция для решения системы методом Якоби
vector<double> solveUsingJacobi(const vector<vector<double>> A, const vector<double> b, double tol = 1e-6, int max_iter = 100) {
    int n = A.size();
    vector<double> x(n, 0);
    vector<double> x_new(n, 0);
    
    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < n; ++i) {
            double sum = 0;
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    sum += A[i][j] * x[j];
                }
            }
            x_new[i] = (b[i] - sum) / A[i][i];
        }
        
        // Проверка сходимости
        bool converged = true;
        for (int i = 0; i < n; ++i) {
            if (fabs(x_new[i] - x[i]) > tol) {
                converged = false;
                break;
            }
        }
        
        if (converged) {
            return x_new;
        }
        
        x = x_new;
    }
    
    return x_new;
}

// Функция для решения системы методом прогонки для трехдиагональной матрицы
vector<double> solveUsingTridiagonal(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    
    // Прямой ход
    for (int i = 1; i < n; ++i) {
        double m = A[i][i - 1] / A[i - 1][i - 1];
        A[i][i] -= m * A[i - 1][i];
        b[i] -= m * b[i - 1];
    }
    
    // Обратный ход
    vector<double> x(n);
    x[n - 1] = b[n - 1] / A[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = (b[i] - A[i][i + 1] * x[i + 1]) / A[i][i];
    }
    
    return x;
}

void print_matrix(vector<vector<double>> matrix) {
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < matrix.at(0).size(); ++j) {
            cout << matrix.at(i).at(j) << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void print_vector(vector<double> vec) {
    for (int i = 0; i < vec.size(); ++i) {
        cout << "x[" << i << "] = " << vec.at(i) << endl;
    }
    cout << endl;
}

int main(int argc, char* argv[]) {
    // Пример трехдиагональной матрицы
    vector<vector<double>> A = {{2, 1, 0},
                                {1, 2, 1},
                                {0, 1, 2}};
    vector<double> b = {1, 0, 1};

    vector<double> solution = solveUsingTridiagonal(A, b);
    cout << "System solve sweep method:" << endl;
    print_vector(solution);
    vector<double> solution2 = solveUsingGaussSeidel(A, b);
    cout << "System solve Gauss-Seidel method:" << endl;
    print_vector(solution2);
    vector<double> solution3 = solveUsingJacobi(A, b);
    cout << "System solve Jacobi method:" << endl;
    print_vector(solution3);
    vector<double> solution4 = solveUsingLU(A, b);
    cout << "System solve LU method:" << endl;
    print_vector(solution4);
    vector<double> solution5 = solveUsingQR(A, b);
    cout << "System solve QR method:" << endl;
    print_vector(solution5);

    return 0;
}