#include <stdio.h>
#include <stdlib.h>

float** criarMatriz(int tamanho)
{
    float **matriz = malloc(tamanho * sizeof(*matriz));
    for (int i = 0; i < tamanho; i++) {
        matriz[i] = malloc(tamanho * sizeof(*matriz[i]));
    }
    return matriz;
}

float** destruirMatriz(float **matriz, int tamanho)
{
    for (int i = 0; i < tamanho; i++) {
        free(matriz[i]);
    }
    free(matriz);
}

void lerMatriz(float **matriz, int tamanho)
{
    for (int i = 0; i < tamanho; i++) {
        for (int j = 0; j < tamanho; j++) {
            scanf("%f", &matriz[i][j]);
        }
    }
}

void printMatriz(float **matriz, int tamanho)
{
    for (int i = 0; i < tamanho; i++) {
        for (int j = 0; j < tamanho; j++) {
            printf("%.*f\t", 3, matriz[i][j]);
        }
        printf("\n");
    }
}

int main()
{
    int tamanho;
    scanf("%d", &tamanho);
    float **matriz = criarMatriz(tamanho);
    lerMatriz(matriz, tamanho);
    printMatriz(matriz, tamanho);

    destruirMatriz(matriz, tamanho);
    return 0;
}