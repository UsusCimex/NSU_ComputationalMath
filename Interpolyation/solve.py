import matplotlib.pyplot as plt
import numpy as np
import sys

def func(x):
    return np.sqrt(x) + np.cos(x)

def lagrange_interpolation(X, Y):
    def interpolated_function(x):
        result = 0
        for i in range(len(X)):
            temp = Y[i]
            for j in range(len(X)):
                if i != j:
                    temp *= (x - X[j]) / (X[i] - X[j])
            result += temp
        return result
    return interpolated_function

print("Enter 4 param: FROM TO NUMBERNODES ERRORDOT")
inp = input().split()
if len(inp) == 4:
    a = float(inp[0])
    if a < 0:
        sys.exit()
    b = float(inp[1])
    if b < a:
        sys.exit()

    numDotOfFunc = int((b - a) * 100)
    X = np.linspace(a + (b-a)*0.25, b - (b-a)*0.25, numDotOfFunc) #0.25 because there are jumps at the edges of the graph
    Y = func(X)
    plt.plot(X, Y, label="Основная функция", linewidth=4)

    numNodes = int(inp[2])
    if numNodes < 0:
        sys.exit()
    errorDot = float(inp[3])
    if errorDot < a or errorDot > b:
        sys.exit()
    plt.scatter(errorDot, func(errorDot), s=50)

    lagX = np.linspace(a, b, numNodes)
    lagY = func(lagX)
    print(lagX)
    print(lagY)
    li = lagrange_interpolation(lagX, lagY)

    lagX2 = np.linspace(a, b, numNodes * 2)
    lagY2 = func(lagX2)
    print(lagX2)
    print(lagY2)
    li2 = lagrange_interpolation(lagX2, lagY2)

    Y2 = li(X)
    plt.plot(X, Y2, label=f"Лагранж на {numNodes} точках")
    plt.scatter(errorDot, li(errorDot), s=40)
    Y3 = li2(X)
    plt.plot(X, Y3, label=f"Лагранж на {numNodes * 2} точках")
    plt.scatter(errorDot, li2(errorDot), s=40)

    error = abs(func(errorDot) - li(errorDot))
    print(f"Погрешность между основной функцией и полученной по методу Лагранжа на {numNodes} точках: {error}.")

plt.legend()
plt.show()
