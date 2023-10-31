#include <iostream>
#include <math.h>
#include <vector>
#include <string>

double cubicEquation(double a, double b, double c, double d, double x) {
    return (a * x * x * x) + (b * x * x) + (c * x) + (d);
}

double quadraticEquation(double a, double b, double c, double x) {
    return (a * x * x) + (b * x) + (c);
}

std::vector<double> solveQuadraticEquation(double a, double b, double c) { //Method : DISCRIMINANT
    double discriminant = b * b - 4 * a * c;
    std::vector<double> result;
    if (discriminant > 0) {
        result.push_back((-b - sqrtf64(discriminant)) / (2 * a));
        result.push_back((-b + sqrtf64(discriminant)) / (2 * a));
    } else if (discriminant == 0) {
        result.push_back((-b) / (2 * a));
    }
    return result;
}

double findSolve(double a, double b, double c, double d, double eps, double left, double right) {
    if (left > right) {
        std::swap(left, right);
    }
    while (right - left > eps) {
        double middle = left + (right - left) / 2;
        double valueOfMiddle = cubicEquation(a, b, c, d, middle);
        if (valueOfMiddle < 0) {
            if (cubicEquation(a, b, c, d, left) > 0) {
                right = middle;
            } else {
                left = middle;
            }
        } else if (valueOfMiddle > 0) {
            if (cubicEquation(a, b, c, d, left) > 0) {
                left = middle;
            } else {
                right = middle;
            }
        } else {
            return middle;
        }
    }
    return left + (right - left) / 2;
}

std::vector<double> solveCubicEquation(double a, double b, double c, double d, double delta, double eps) { //Method : DICHOTOMY
    std::vector<double> extremum = solveQuadraticEquation(3 * a, 2 * b, c);
    std::vector<double> result;
    if (extremum.size() == 2) {
        if (quadraticEquation(3 * a, 2 * b, c, extremum[0] - delta) * cubicEquation(a, b, c, d, extremum[0]) > 0) {
            result.push_back(findSolve(a, b, c, d, eps, extremum[0] - delta, extremum[0]));
        }
        if (cubicEquation(a, b, c, d, extremum[0]) * cubicEquation(a, b, c, d, extremum[1]) < 0) {
            result.push_back(findSolve(a, b, c, d, eps, extremum[0], extremum[1]));
        }
        if (quadraticEquation(3 * a, 2 * b, c, extremum[1] + delta) * cubicEquation(a, b, c, d, extremum[1]) < 0) {
            result.push_back(findSolve(a, b, c, d, eps, extremum[1], extremum[1] + delta));
        }
        if (cubicEquation(a, b, c, d, extremum[0]) == 0) {
            result.push_back(extremum[0]);
        }
        if (cubicEquation(a, b, c, d, extremum[1]) == 0) {
            result.push_back(extremum[1]);
        }
    } else if (extremum.size() == 1) {
        if (quadraticEquation(3 * a, 2 * b, c, extremum[0] - delta) * cubicEquation(a, b, c, d, extremum[0]) > 0) {
            result.push_back(findSolve(a, b, c, d, eps, extremum[0] - delta, extremum[0]));
        }
        if (quadraticEquation(3 * a, 2 * b, c, extremum[0] + delta) * cubicEquation(a, b, c, d, extremum[0]) < 0) {
            result.push_back(findSolve(a, b, c, d, eps, extremum[0] + delta, extremum[0]));
        }
        if (cubicEquation(a, b, c, d, extremum[0]) == 0) {
            result.push_back(extremum[0]);
        }
    } else {
        result.push_back(findSolve(a, b, c, d, eps, -delta, delta));
    }
    return result;
}

void printArray(std::vector<double> array) {
    for (int i = 0; i < array.size(); ++i) {
        std::cout << array[i] << " ";
    }
    std::cout << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cout << "Usage: " << argv[0] << " a b c d" << std::endl;
        return 1;
    }
    double a = std::stod(argv[1]);
    double b = std::stod(argv[2]);
    double c = std::stod(argv[3]);
    double d = std::stod(argv[4]); 
    std::vector<double> result = solveCubicEquation(a, b, c, d, 1e5, 1e-6);
    printArray(result);
    return 0;
}