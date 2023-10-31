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
    double result = 0;
    double step = (to - from) / numNodes;
    for (int i = 0; i < numNodes; ++i) {
        double a = from + i * step;
        double b = from + (i + 1) * step;
        result += func((a + b) / 2) * step;
    }
    return result;
}

// Trapezoidal formula
double integral2(double from, double to, integralFunction func, int numNodes) {
    double result = 0;
    double step = (to - from) / numNodes;
    for (int i = 0; i < numNodes; ++i) {
        double a = from + i * step;
        double b = from + (i + 1) * step;
        result += (b - a) / 2 * (func(a) + func(b)); ;
    }
    return result;
}

// Parabola formula (Sympson)
double integral3(double from, double to, integralFunction func, int numNodes) {
    double result = 0;
    double step = (to - from) / numNodes;
    for (int i = 0; i < numNodes; ++i) {
        double a = from + i * step;
        double b = from + (i + 1) * step;
        result += (b - a) / 6 * (func(a) + 4 * func((a + b) / 2) + func(b));
    }
    return result;
}

// Three-eighths formula
double integral4(double from, double to, integralFunction func, int numNodes) {
    double result = 0;
    double step = (to - from) / numNodes;
    for (int i = 0; i < numNodes; ++i) {
        double a = from + i * step;
        double b = from + (i + 1) * step;
        result += (b - a) / 8 * (func(a) + 3 * func(a + (b - a) / 3) + 3 * func(b - (b - a) / 3) + func(b));
    }
    return result;
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
    double realIntegral, calculation;

    std::cout << "Function = ln(1+x)" << std::endl;
    realIntegral = integral(from, to, func1, (to - from) * 1000);
    std::cout << "Real integral result: " << realIntegral << std::endl;
    calculation = integral(from, to, func1, numNodes);
    std::cout << "Rectangle integral result: " << calculation << " (" << 
                                                  abs(realIntegral - calculation) << ")" << std::endl;
    calculation = integral2(from, to, func1, numNodes);
    std::cout << "Trapezoidal integral result: " << calculation << " (" <<
                                                    abs(realIntegral - calculation) << ")" << std::endl;
    calculation = integral3(from, to, func1, numNodes);
    std::cout << "Parabola(Simpson) integral result: " << calculation << " (" <<
                                                 abs(realIntegral - calculation) << ")" << std::endl;
    calculation = integral4(from, to, func1, numNodes);
    std::cout << "Three-eighths integral result: " << calculation << " (" <<
                                                      abs(realIntegral - calculation) << ")" << std::endl;

    p = 2;
    std::cout << "Runge integral error (rectangle): " << abs(integral(from, to, func1, numNodes / 2) -
                                                 integral(from, to, func1, numNodes)) / (pow(2, p) - 1) << std::endl;
    p = 2;
    std::cout << "Runge integral error (trapezoidal): " << abs(integral2(from, to, func1, numNodes / 2) -
                                                 integral2(from, to, func1, numNodes)) / (pow(2, p) - 1) << std::endl;
    p = 4;
    std::cout << "Runge integral error (Simpson): " << abs(integral3(from, to, func1, numNodes / 2) -
                                                 integral3(from, to, func1, numNodes)) / (pow(2, p) - 1) << std::endl << std::endl;
    
    std::cout << "Function = e^x * cos(x)" << std::endl;
    realIntegral = integral(from, to, func2, (to - from) * 1000);
    std::cout << "Real integral result: " << realIntegral << std::endl;
    calculation = integral(from, to, func2, numNodes);
    std::cout << "Rectangle integral result: " << calculation << " (" << 
                                                  abs(realIntegral - calculation) << ")" << std::endl;
    calculation = integral2(from, to, func2, numNodes);
    std::cout << "Trapezoidal integral result: " << calculation << " (" <<
                                                    abs(realIntegral - calculation) << ")" << std::endl;
    calculation = integral3(from, to, func2, numNodes);
    std::cout << "Parabola(Simpson) integral result: " << calculation << " (" <<
                                                 abs(realIntegral - calculation) << ")" << std::endl;
    calculation = integral4(from, to, func2, numNodes);
    std::cout << "Three-eighths integral result: " << calculation << " (" <<
                                                      abs(realIntegral - calculation) << ")" << std::endl;

    p = 2;
    std::cout << "Runge integral error (rectangle): " << abs(integral(from, to, func2, numNodes / 2) -
                                                 integral(from, to, func2, numNodes)) / (pow(2, p) - 1) << std::endl;
    p = 2;
    std::cout << "Runge integral error (trapezoidal): " << abs(integral2(from, to, func2, numNodes / 2) -
                                                 integral2(from, to, func2, numNodes)) / (pow(2, p) - 1) << std::endl;
    p = 4;
    std::cout << "Runge integral error (Simpson): " << abs(integral3(from, to, func2, numNodes / 2) -
                                                 integral3(from, to, func2, numNodes)) / (pow(2, p) - 1) << std::endl << std::endl;
    return 0;
}