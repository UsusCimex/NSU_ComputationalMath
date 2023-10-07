#include <iostream>
#include <math.h>

double func1(double x) {
    return logf(1 + x);
}

double func2(double x) {
    return expf(x)*cosf(x);
}

typedef double (*integralFunction)(double);
// Rectangle formula
double integral(double from, double to, integralFunction func, int numNodes) {
    double area = 0;
    double step = (to - from) / numNodes;
    for (int i = 0; i < numNodes; ++i) {
        double a = from + i * step;
        double b = from + (i + 1) * step;
        area += func((a + b) / 2) * step;
    }
    return area;
}

// Trapezoidal formula
double integral2(double from, double to, integralFunction func) {
    return (to - from) / 2 * (func(from) + func(to)); 
}

// Parabola formula (Sympson)
double integral3(double from, double to, integralFunction func) {
    return (to - from) / 6 * (func(from) + 4 * func((from + to) / 2) + func(to)); 
}

// Three-eighths formula
double integral4(double from, double to, integralFunction func) {
    return (to - from) / 8 * (func(from) + 3 * func(from + (to - from) / 3) + 3 * func(to - (to - from) / 3) + func(to));
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cout << "Usage: " << argv[0] << " from to numNodes" << std::endl;
        return 1;
    }
    double from = std::stoi(argv[1]);
    double to = std::stoi(argv[2]);
    int numNodes = std::stoi(argv[3]);
    int p; //Order of accuracy

    std::cout << "Function = ln(1+x)" << std::endl;
    std::cout << "Real integral result: " << integral(from, to, func1, (to - from) * 1000) << std::endl;
    std::cout << "Rectangle integral result: " << integral(from, to, func1, numNodes) << std::endl;
    std::cout << "Trapezoidal integral result: " << integral2(from, to, func1) << std::endl;
    std::cout << "Parabola integral result: " << integral3(from, to, func1) << std::endl;
    std::cout << "Three-eighths integral result: " << integral4(from, to, func1) << std::endl;

    std::cout << "Rectangle integral error: " << abs(integral(from, to, func1, (to - from) * 1000) -
                                                     integral(from, to, func1, numNodes)) << std::endl;
    p = 2;
    std::cout << "Runge integral error: " << abs(integral(from, to, func1, numNodes / 2) -
                                                 integral(from, to, func1, numNodes)) / (pow(2, p) - 1) << std::endl << std::endl;
    
    std::cout << "\nFunction = e^x cos(x)" << std::endl;
    std::cout << "Real integral result: " << integral(from, to, func2, (to - from) * 1000) << std::endl;
    std::cout << "Rectangle integral result: " << integral(from, to, func2, numNodes) << std::endl;
    std::cout << "Trapezoidal integral result: " << integral2(from, to, func2) << std::endl;
    std::cout << "Parabola integral result: " << integral3(from, to, func2) << std::endl;
    std::cout << "Three-eighths integral result: " << integral4(from, to, func2) << std::endl;

    std::cout << "Rectangle integral error: " << abs(integral(from, to, func2, (to - from) * 1000) -
                                                     integral(from, to, func2, numNodes)) << std::endl;
    p = 2;
    std::cout << "Runge integral error: " << abs(integral(from, to, func2, numNodes / 2) -
                                                 integral(from, to, func2, numNodes)) / (pow(2, p) - 1) << std::endl << std::endl;
    return 0;
}