# EP1 de Numérico 2021

Esta pasta contém o código em Python do algorítmo QR com deslocamento
espectral para o cálculo dos auto-valores e auto-vetores de uma matriz
tridiagonal simétrica n × n.

Nesta implementação, foi priorizada a simplicidade do código em detrimento da
eficiência e velocidade de execução.

# Implementação

## Decomposição QR
### **getGivensMatrix(matriz, linha)**
Esta função recebe uma matriz quadrada e o número do elemento da diagonal
principal (matriz[linha][linha]). A função retorna uma matriz de rotação de
Givens que zera o valor do elemento abaixo do elemento da diagonal principal
(matriz[linha][linha+1]).

### **getQRDecomposition(matriz)**
Esta função recebe uma matriz tridiagonal simétrica e retorna uma tupla
contendo as matrizes Q e R da decomposição QR.

## Algoritmo QR com deslocamento espectral
### **getMik(matriz)**
Calcula o múltiplo da identidade µk que será utilizado no algoritmo QR com
deslocamento espectral.