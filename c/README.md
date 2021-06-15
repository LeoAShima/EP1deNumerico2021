# EP1 de Numérico 2021

Esta pasta contém o código em C do algorítmo QR com deslocamento
espectral para o cálculo dos auto-valores e auto-vetores de uma matriz
tridiagonal simétrica n × n.

# Implementação

## Matrizes
### **criarMatriz(int tamanho)**
Esta função aloca dinamicamente uma matriz quadrada *tamanho* x *tamanho* de
*doubles*.
| Nome | Descrição |
| ---- | --------- |
| tamanho | *int*<br>Tamanho da matriz n x n a ser alocada|
#### Retorna
Um ponteiro de ponteiro de doubles. Pode ser utilizada como uma matriz.
Deve ser destruído com *destruirMatriz* depois de usado.

### **destruirMatriz(double \*\*matriz, int tamanho)**
Esta função libera a memória da matriz.
| Nome | Descrição |
| ---- | --------- |
| matriz | *double\*\**<br>Ponteiro da matriz alocada com criarMatriz |
| tamanho | *int*<br>Tamanho da matriz a ser desalocada |

### **lerMatriz(double \*\*matriz, int tamanho)**
Esta função lê os valores a serem inseridos na matriz passada como argumento.
São lidos *tamanho* x *tamanho* valores *double*.
| Nome | Descrição |
| ---- | --------- |
| matriz | *double\*\**<br>Ponteiro da matriz a ser escrita |
| tamanho | *int*<br>Tamanho da matriz a ser escrita |

### **printMatriz(double \*\*matriz, int tamanho, int precisao)**
Esta função imprime os valores armazenados na matriz.
| Nome | Descricao |
| ---- | --------- |
| matriz | *double\*\**<br>Ponteiro da matriz a ser impressa |
| tamanho | *int*<br>Tamanho da matriz a ser impressa |
| precisao | *int*<br>Quantidade de digitos após a vírgula a serem exibidos |


## Decomposição QR
### **getQRDecomposition(double \*\*matriz, int tamanho, double \*c, double \*s)**
Esta função sobrescreve a matriz tridiagonal simétrica recebida com a matriz R
resultante de sua decomposição QR. São armazenados nos vetores *c* e *s* os
parâmetros das rotações de Givens utilizados para se obter Q.
| Nome | Descrição |
| ---- | --------- |
| matriz | *double\*\**<br>Ponteiro para matriz a ser decomposta. Deve ser tridiagonal simétrica. Ao final da execução, conterá a matriz R. |
| tamanho | *int*<br>Tamanho da matriz a ser decomposta. |
| c | *double\**<br>Vetor que armazenará os parâmetros c das rotações de Givens da decomposição. |
| s | *double\**<br>Vetor que armazenará os parâmetros s das rotações de Givens da decomposição. |

### **getCSParameters(double alpha, double beta, double \*c, double \*s)**
Esta função recebe o valor da diagonal principal *alpha* e o valor da
subdiagonal inferior *beta*, e calcula os parâmetros *c* e *s* da rotação de
Givens que zera o valor de *beta*.
| Nome | Descrição |
| ---- | --------- |
| alpha | *double*<br>Valor da diagonal principal. |
| beta | *double*<br>Valor da subdiagonal inferior. |
| c | *double\**<br>Ponteiro da variável que conterá o parâmetro c calculado. |
| s | *double\**<br>Ponteiro da variável que conterá o parâmetro s calculado. |

### **applyGivens(double \*\*matriz, int tamanho, int i, double c, double s)**
Esta função utiliza os parâmetros c e s para aplicar a rotação de Givens na
matriz recebida, alterando as linhas i e i+1.
| Nome | Descrição |
| ---- | --------- |
| matriz | *double\*\**<br>Ponteiro da matriz a ser alterada. |
| tamanho | *int*<br>Tamanho da matriz a ser alterada. |
| i | *int*<br>Linha da matriz que será alterada, juntamente com i+1. |
| c | *double*<br>Parâmetro *c* da rotação de Givens a ser aplicada. |
| s | *double*<br>Parâmetro *s* da rotação de Givens a ser aplicada. |


## Algoritmo QR
### **getMiK(double \*\*matriz, int tamanho)**
Esta função calcula o valor μ<sub>k</sub> a ser subtraído da diagonal principal
no algoritmo QR para acelerar a convergência.
| Nome | Descrição |
| ---- | --------- |
| matriz | *double\*\**<br>Ponteiro da matriz. |
| tamanho | *int*<br>Tamanho da matriz a ser considerada.<br>Se o valor for menor que o tamanho total da matriz, serão utilizados para o cálculo apenas o subconjunto de elementos dentro da matriz de dimensão *tamanho* x *tamanho*. |