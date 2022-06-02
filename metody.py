import time
import matplotlib.pyplot as plt
import copy

from funkcje import *

def jacobi(A, b, t):
    # Jacobi: A = D + L + U; x(k+1) = D^-1(b+(L+U)x(k))
    # Jacobi element-based: x_i(k+1) = 1/A_ii * (b_i - SUM(j != i)(a_ij * x_j(k))); i = 1..N
    
    N = len(A)
    xPrev = initializeMat(1, N, 1)
    xNext = initializeMat(1, N, 1)
    norma = 1
    normy = []
    iter = 0
    start = time.perf_counter()
    iterLim = 10000

    while norma > 10e-9:
        iter += 1
        if iter == iterLim:
            print("iterLim")
            break
        for i in range(N):
            suma = 0
            for j in range(N):
                if j != i:
                    suma += A[i][j] * xPrev[j][0]
            xNext[i][0] = 1/A[i][i] * (b[i][0] - suma)
        xPrev = xNext
        # res = (A*x - b)
        # norm = sqrt(SUM(res[i] ^ 2))
        norma = norm(subMat(multiplyMat(A, xNext), b))
        normy.append(norma)

    end = time.perf_counter()
    t[0] = end - start
    plt.figure()
    plt.plot(normy)
    plt.title("Wykres zbieznosci: Jacobi (N={0})".format(N))
    plt.xlabel("iteracja")
    plt.ylabel("norma z residuum")
    # plt.show(block = False)
    plt.savefig("JacobiRes{0}{1}".format(N, A[0][0]))
    print("Jacobi:"
        + "\n\titeracji: {0}".format(iter)
        + "\n\tnorma z residuum: {0}".format(norma)
        + "\n\tczas wykonywania: {0}s".format(t[0])
        + "\n\tN: {0}".format(N))

def gaus_seidel(A, b, t):
    # Gaus-Seidel
    # decomposition: A = L(*) + U (lower triangular + strictly upper triangular)
    # element-wise formula:
    # x(k+1)_i = 1/a_ii * (b_i - SUM(j = 1 -> i-1)(a_ijx(k+1)_j) - SUM(j=i+1 -> N)(a_ijx(k)_j)); i = 1..N
    # L(*)x(k+1) = b - Ux(k)
    # x(k+1) = L(*)^-1 * (b - U*x(k))

    # w metodzie G-S wystarcza jeden wektor

    N = len(A)
    x = initializeMat(1, N, 1)
    norma = 1
    normy = []
    start = time.perf_counter()
    iterLim = 10000
    iter = 0

    while norma > 10e-9:
        iter += 1
        if iter == iterLim:
            print("iterLim")
            break
        for i in range(N):
            suma = 0
            for j in range(N):
                if j != i:
                    suma += A[i][j] * x[j][0]
            x[i][0] = 1 / A[i][i] * (b[i][0] - suma)
        norma = norm(subMat(multiplyMat(A, x), b))
        normy.append(norma)

    end = time.perf_counter()
    t[0] = end - start
    plt.figure()
    plt.plot(normy)
    plt.title("Wykres zbieznosci: Gauss-Seidel (N={0})".format(N))
    plt.xlabel("iteracja")
    plt.ylabel("norma z residuum")
    # plt.show(block = False)
    plt.savefig("GaussRes{0}{1}".format(N, A[0][0]))
    print("\nGauss-Seidel:"
        + "\n\titeracji: {0}".format(iter)
        + "\n\tnorma z residuum: {0}".format(norma)
        + "\n\tczas wykonywania: {0}s".format(t[0])
        + "\n\tN: {0}".format(N))

    # printMat(subMat(np.linalg.solve(A, b), x))

def lu(A, b, t):
    # LU
    # Ax = b -> LUx = b
    # wektor pomocniczy y = Ux
    # uklad rownan Ly = b (podstawianie wprzod)
    # uklad rownan Ux = y (podstawianie wstecz)

    # A -> LU
    N = len(A)
    U = copy.deepcopy(A)
    L = identity(N)
    start = time.perf_counter()
    for k in range(N - 1):
        for j in range(k+1,N):
            L[j][k] = U[j][k]/U[k][k]
            for i in range(k, N):
                U[j][i] = U[j][i] - L[j][k]*U[k][i]

    # Ly = b -> y (forward substitution)
    y = initializeMat(1, N, 0)
    for i in range(N):
        suma = b[i][0]
        for j in range(i): # i - 1
            suma -= L[i][j] * y[j][0]
        y[i][0] = suma / L[i][i]

    # Ux = y -> x (backward substitution)
    x = initializeMat(1, N, 1)
    for i in range(N-1, -1, -1):
        suma = y[i][0]
        for j in range(i + 1, N):
            suma -= U[i][j] * x[j][0]
        x[i][0] = suma / U[i][i]

    end = time.perf_counter()
    t[0] = end - start
    print("\nLU"
        + "\n\tnorma z residuum: {0}".format(norm(subMat(multiplyMat(A, x), b)))
        + "\n\tczas wykonywania: {0}s".format(t[0])
        + "\n\tN: {0}".format(len(A)))