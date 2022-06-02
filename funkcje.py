from math import sin, sqrt


def initializeMat(cols, rows, val):
    ret = []
    for i in range(rows):
        ret.append([])
        for _ in range(cols):
            ret[i].append(val)
    return ret

def printMat(m):
    s = ""
    for row in m:
        s += " "
        for el in row:
            s += str(round(el, 2)) + "\t"
        s = s.removesuffix("\t")
        s += "\n"
    s = s.removesuffix("\n")
    s += "]"
    # s.removeprefix(" ")
    # s = "[" + s
    s.replace(" ", "[", 1)
    print(s)
        
def multiplyMat(A, B):
    # c_ij = SUM(k=1 -> N)(a_ik*b_kj)
    C = initializeMat(len(B[0]), len(A), 0)
    for i in range(len(A)):
        for j in range(len(B[0])):
            suma = 0
            for k in range(len(B)):
                suma += A[i][k] * B[k][j]
            C[i][j] = suma
    return C

def subMat(A, B):
    C = initializeMat(len(A[0]), len(A), 0)
    for i in range(len(A)):
        for j in range(len(A[0])):
            C[i][j] = A[i][j] - B[i][j]
    return C

def identity(N):
    I = []
    for i in range(N):
        I.append([])
        for j in range(N):
            I[i].append( 1 if i == j else 0 )
    return I

def norm(vec):
    suma = 0
    for i in range(len(vec)):
        suma += pow(vec[i][0], 2)
    return sqrt(suma)
    
def initializeA(N, a1, a2, a3):
    A = initializeMat(N, N, 0)
    for row in range(N):
        for col in range(N):
            if row == col:
                A[row][col] = a1
            elif abs(row - col) == 1:
                A[row][col] = a2
            elif abs(row - col) == 2:
                A[row][col] = a3
            # else:
            #      A[row][col] = 0
    return A

def initializeB(N):
    b = initializeMat(1, N, 0)
    for i in range(N):
        b[i][0] = sin(i * 5)
    return b