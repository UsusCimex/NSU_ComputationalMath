import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib.pyplot as plt
import math
import sys

# Определение основной функции
def mainFunction(x):
    return math.sqrt(x) + math.cos(x)

# Функция интерполяции Лагранжа
def lagrangeInterpolation(X):
    def interpolatedFunction(x):
        result = 0
        for i in range(len(X)):
            temp = mainFunction(X[i])
            for j in range(len(X)):
                if i != j:
                    temp *= (x - X[j]) / (X[i] - X[j])
            result += temp
        return result
    return interpolatedFunction

# Вывод приветствия и получение параметров от пользователя
print("Enter 4 params: FROM TO NUMBERNODES ERRORDOT")
params = input().split()
if len(params) != 4:
    print("Invalid number of parameters")
    sys.exit()

# Обработка параметров
a = float(params[0])
if a < 0:
    print("Invalid parameter value: FROM")
    sys.exit()
b = float(params[1])
if b < a:
    print("Invalid parameter value: TO")
    sys.exit()
numNodes = int(params[2])
if numNodes < 0:
    print("Invalid parameter value: NUMBERNODES")
    sys.exit()
errorDot = float(params[3])
if errorDot < a or errorDot > b:
    print("Invalid parameter value: ERRORDOT")
    sys.exit()
plt.scatter(errorDot, mainFunction(errorDot), s=50)

# Вычисление точек основной функции
numDotOfFunc = int((b - a) * 100)
X = []
Y = []
mainStep = (b - a) / numDotOfFunc
for i in range(numDotOfFunc):
    X.append(a + mainStep * i)
    Y.append(mainFunction(X[i]))
plt.plot(X, Y, label="Main Function", linewidth=4)

# Создание наборов точек для интерполяции Лагранжа
lagX = []
step = (b - a) / numNodes
for i in range(numNodes):
    lagX.append(a + step * i)
interpolatedFunc1 = lagrangeInterpolation(lagX)

lagX2 = []
step = (b - a) / (numNodes * 2)
for i in range(numNodes * 2):
    lagX2.append(a + step * i)
interpolatedFunc2 = lagrangeInterpolation(lagX2)

# Вычисление и отображение интерполированных функций
Y2 = []
for i in range(numDotOfFunc):
    Y2.append(interpolatedFunc1(a + mainStep * i))
plt.plot(X, Y2, label=f"Lagrange Interpolation with {numNodes} Nodes")
plt.scatter(errorDot, interpolatedFunc1(errorDot), s=40)
Y3 = []
for i in range(numDotOfFunc):
    Y3.append(interpolatedFunc2(a + mainStep * i))
plt.plot(X, Y3, label=f"Lagrange Interpolation with {numNodes * 2} Nodes")
plt.scatter(errorDot, interpolatedFunc2(errorDot), s=40)

# Вычисление и вывод погрешности
error = abs(interpolatedFunc1(errorDot) - mainFunction(errorDot))
print(f"Error between the main function and the Lagrange interpolation with {numNodes} nodes: {error}.")

# Вывод легенды и отображение графика
plt.legend()
plt.show()
