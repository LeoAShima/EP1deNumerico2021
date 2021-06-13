# EP1 de Numérico 2021

Esta pasta contém o código em C do algorítmo QR com deslocamento
espectral para o cálculo dos auto-valores e auto-vetores de uma matriz
tridiagonal simétrica n × n.

# Implementação

## Decomposição QR
### **criarMatriz(int tamanho)**
Esta função aloca dinamicamente uma matriz quadrada *tamanho* x *tamanho* de
*doubles*.
| Nome | Descrição |
| ---- | --------- |
| tamanho | *int*<br>Tamanho da matriz n x n a ser alocada|
#### Retorna
Um ponteiro de ponteiro de doubles. Pode ser utilizada como uma matriz.
Deve ser destruído com *destruirMatriz* depois de usado.

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