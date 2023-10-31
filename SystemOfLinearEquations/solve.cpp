#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>

using namespace std;

void PrintMatrix(vector<vector<double>>& matrix) {
    for (int i = 0; i < matrix.size(); ++i) {
        cout << "|";
        for (int j = 0; j < matrix.at(0).size(); ++j) {
            cout << std::setw(8) << std::setprecision(3) << matrix.at(i).at(j) << " ";
        }
        cout << "|" << endl;
    }
    cout << endl;
}

void PrintVector(vector<double>& vec) {
    for (int i = 0; i < vec.size(); ++i) {
        cout << "x[" << i << "] = " << vec.at(i) << endl;
    }
    cout << endl;
}

void TransponseMatrix(vector<vector<double>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();
    
    vector<vector<double>> transposedMatrix(cols, vector<double>(rows));
    
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            transposedMatrix[j][i] = matrix[i][j];
        }
    }
    
    matrix = transposedMatrix;
}

// LU-разложение матрицы A
void LU_Decomposition(vector<vector<double>>& A, vector<vector<double>>& L, vector<vector<double>>& U) {
    int n = A.size();

    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double sum = 0;
            for (int k = 0; k < i; ++k) {
                sum += (L[i][k] * U[k][j]);
            }
            U[i][j] = A[i][j] - sum;
        }
        for (int j = i; j < n; ++j) {
            if (i == j) {
                L[i][i] = 1;
            } else {
                double sum = 0;
                for (int k = 0; k < i; ++k) {
                    sum += (L[j][k] * U[k][i]);
                }
                L[j][i] = (A[j][i] - sum) / U[i][i];
            }
        }
    }

    cout << "L Matrix:" << endl;
    PrintMatrix(L);
    cout << "U Matrix:" << endl;
    PrintMatrix(U);
}

// Решение СЛАУ при помощи LU-разложения
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

// QR-разложение матрицы A (метод Холецкого)
bool QR_Decomposition_Cholesky(vector<vector<double>>& A, vector<vector<double>>& Q, vector<vector<double>>& R) {
    int n = A.size();
    
    for (int i = 0; i < n; ++i) {
        //Сначала вычисляем значения элементов слева от диагонального элемента
        for (int j = 0; j < i; ++j)
        {
            double sum = 0;
            for (int k = 0; k < j; ++k)
            {
                sum += Q[i][k] * Q[j][k];
            }
            Q[i][j] = (A[i][j] - sum) / Q[j][j];
        }

        //Находим значение диагонального элемента
        double temp = A[i][i];
        for (int k = 0; k < i; ++k)
        {
            temp -= Q[i][k] * Q[i][k];
        }
        Q[i][i] = sqrtl(temp);
    }

    R = Q;
    TransponseMatrix(R);
    
    cout << "Q Matrix:" << endl;
    PrintMatrix(Q);
    cout << "R Matrix:" << endl;
    PrintMatrix(R);

    return true;
}

// Решение системы методом QR-разложения
vector<double> solveUsingQR(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    vector<vector<double>> Q(n, vector<double>(n, 0));
    vector<vector<double>> R(n, vector<double>(n, 0));
    
    if (!QR_Decomposition_Cholesky(A, Q, R)) {
        cout << "Choletsky decomposition failed (the matrix was not positively definite)..." << endl;
        return vector<double>(n, 0);
    }

    // Решение системы Qy = b
    vector<double> y(n, 0);
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < i; ++j) {
            sum += Q[i][j] * y[j];
        }
        y[i] = (b[i] - sum) / Q[i][i];
    }

    // Решение системы Rx = y
    vector<double> x(n, 0);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; ++j) {
            sum += R[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / R[i][i];
    }

    return x;
}

// Функция для решения системы методом Гаусса-Зейделя
vector<double> solveUsingGaussSeidel(const vector<vector<double>> A, const vector<double> b, double eps = 1e-8, int max_iter = 100) {
    int n = A.size();
    vector<double> x(n, 0);
    vector<double> x_new(n, 0);
    
    for (int iter = 0; iter < max_iter; ++iter) {
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
            if (fabs(x_new[i] - x[i]) > eps) {
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
vector<double> solveUsingJacobi(const vector<vector<double>> A, const vector<double> b, double eps = 1e-8, int max_iter = 100) {
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
            if (fabs(x_new[i] - x[i]) > eps) {
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

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (abs(i - j) > 1) {
                if (A[i][j] != 0) {
                    cout << "The matrix is not tridiagonal..." << endl;
                    return vector<double>(n, 0);
                }
            }
        }
    }
    
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

void TestMatrix(vector<vector<double>>& A, vector<double> b) {
    cout << "Matrix A:" << endl;
    PrintMatrix(A);
    cout << "Vector b:" << endl;
    PrintVector(b);

    cout << "System solve sweep method:" << endl;
    vector<double> solution = solveUsingTridiagonal(A, b);
    PrintVector(solution);

    cout << "System solve Gauss-Seidel method:" << endl;
    vector<double> solution2 = solveUsingGaussSeidel(A, b);
    PrintVector(solution2);

    cout << "System solve Jacobi method:" << endl;
    vector<double> solution3 = solveUsingJacobi(A, b);
    PrintVector(solution3);

    cout << "System solve LU method:" << endl;
    vector<double> solution4 = solveUsingLU(A, b);
    PrintVector(solution4);

    cout << "System solve QR method (Cholesky):" << endl;
    vector<double> solution5 = solveUsingQR(A, b);
    PrintVector(solution5);
}

int main(int argc, char* argv[]) {
    vector<vector<double>> A = {{81, -45, 45},
                                {-45, 50, -15},
                                {45, -15, 38}};
    vector<double> b = {531, -460, 193};
    TestMatrix(A, b);
    return 0;
}