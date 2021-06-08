import numpy as np
import math

test = np.array([[2, 1, 0, 0],
                 [1, 2, 1, 0],
                 [0, 1, 2, 1],
                 [0, 0, 1, 2]])

test2 = np.array([[5/math.sqrt(5), 4/math.sqrt(5), 1/math.sqrt(5), 0],
                  [0, 3/math.sqrt(5), 2/math.sqrt(5), 0],
                  [0, 1, 2, 1],
                  [0, 0, 1, 2]])

def qrDecompose(matriz):
    tamanho = len(matriz)
    A = matriz
    givensParameters = []
    
    for i in range(tamanho-1):
        givensParameters.append(getGivensParameters(A[i][i], A[i+1][i]))
        A = applyGivensRotation(A, givensParameters[i][0], givensParameters[i][1], i)
        print(A)

def getGivensParameters(alfa, beta):
    if abs(alfa) > abs(beta):
        tau = -beta/alfa
        c = 1/math.sqrt(1+tau**2)
        s = c*tau
    else:
        tau = -alfa/beta
        s = 1/math.sqrt(1+tau**2)
        c = s*tau
    return c, s

def applyGivensRotation(matriz, c, s, k):
    qa = np.identity(len(matriz))
    qa[k][k] = c
    qa[k][k+1] = -s
    qa[k+1][k] = s
    qa[k+1][k+1] = c
    return qa.dot(matriz)

def qrstep(matriz, k):
    """Calcula a matriz da rotação de Givens para alterar a coluna k"""

    alfa = matriz[k][k]
    beta = matriz[k+1][k]

    if abs(alfa) > abs(beta):
        tau = -beta/alfa
        c = 1/math.sqrt(1+tau**2)
        s = c*tau
    else:
        tau = -alfa/beta
        s = 1/math.sqrt(1+tau**2)
        c = s*tau
    
    q = np.identity(len(matriz))
    q[k][k] = c
    q[k][k+1] = -s
    q[k+1][k] = s
    q[k+1][k+1] = c

    return q

qrDecompose(test)