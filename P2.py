# Jacobi, Gauss-Seidl, LU
# Jakub Kilianczyk 184301
# 03.05.2022

# a1 = 5 + 3
# a2 = a3 = -1
# N = 901
# b[n] = sin(n*(4 + 1))

import matplotlib.pyplot as plt

from funkcje import *
from metody import *


N = 901
t = [1]

A = initializeA(N, 8, -1, -1)
b = initializeB(N)

jacobi(A, b, t)
gaus_seidel(A, b, t)

A = initializeA(N, 3, -1, -1)

gaus_seidel(A, b, t)
jacobi(A, b, t)

# exit()

lu(A, b, t)

# exit()
arrOfN = [100, 500, 1000, 2000, 3000, 4000, 5000]#, 7500, 10000, 15000]
times = [[], [], []]
for n in arrOfN:
    N = n
    A = initializeA(N, 8, -1, -1)
    b = initializeB(N)

    jacobi(A, b, t)
    times[0].append(t[0])
    gaus_seidel(A, b, t)
    times[1].append(t[0])
    lu(A, b, t)
    times[2].append(t[0])

names = ["Jacobi", "Gauss-Seidel", "LU"]
for i in range(len(names)):
    plt.figure()
    plt.plot(arrOfN, times[i])
    plt.title("T(N): {0}".format(names[i]))
    plt.xlabel("N")
    plt.ylabel("czas wykonywania")
    # plt.show(block = False)
    plt.savefig("{0}Times".format(names[i]))