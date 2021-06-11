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

def getGivensMatrix(matriz, linha):
    """Esta função recebe uma matriz quadrada e o número do elemento da diagonal
    principal (matriz[linha][linha]). A função retorna uma matriz de rotação de
    Givens que zera o valor do elemento abaixo do elemento da diagonal principal
    (matriz[linha][linha+1])."""
    alfa = matriz[linha][linha]
    beta = matriz[linha+1][linha]
    if abs(alfa) > abs(beta):
        tau = -beta/alfa
        c = 1/math.sqrt(1+tau**2)
        s = c*tau
    else:
        tau = -alfa/beta
        s = 1/math.sqrt(1+tau**2)
        c = s*tau
    givensRotationMatrixQ = np.identity(len(matriz))
    givensRotationMatrixQ[linha][linha] = c
    givensRotationMatrixQ[linha][linha+1] = -s
    givensRotationMatrixQ[linha+1][linha] = s
    givensRotationMatrixQ[linha+1][linha+1] = c
    return givensRotationMatrixQ

def getQRDecomposition(matriz):
    '''Esta função recebe uma matriz tridiagonal simétrica e retorna uma tupla
    contendo as matrizes Q e R da decomposição QR.'''
    tamanho = len(matriz)
    Qt = np.identity(tamanho) # Esta variável contém a matriz Q transposta
    for i in range(tamanho-1):
        Qi = getGivensMatrix(matriz, i)
        matriz = Qi @ matriz
        Qt = Qi @ Qt
    return Qt.transpose(), matriz

def getMik(matriz):
    '''Calcula o múltiplo da identidade µk que será utilizado no algoritmo QR
    com deslocamento espectral.'''
    n = len(matriz)
    dk = (matriz[n-2][n-2] - matriz[n-1][n-1])/2
    if dk >= 0:
        mik = matriz[n-1][n-1] + dk - math.sqrt(dk**2 + matriz[n-1][n-2]**2)
    else:
        mik = matriz[n-1][n-1] + dk + math.sqrt(dk**2 + matriz[n-1][n-2]**2)
    return mik

Q, R = getQRDecomposition(test)
print(Q)
print(R)
print(Q @ R)