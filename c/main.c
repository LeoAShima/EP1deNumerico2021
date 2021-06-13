#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double** criarMatriz(int tamanho)
{
    double **matriz = malloc(tamanho * sizeof(*matriz));
    for (int i = 0; i < tamanho; i++) {
        matriz[i] = malloc(tamanho * sizeof(*matriz[i]));
    }
    return matriz;
}

double** destruirMatriz(double **matriz, int tamanho)
{
    for (int i = 0; i < tamanho; i++) {
        free(matriz[i]);
    }
    free(matriz);
}

void lerMatriz(double **matriz, int tamanho)
{
    for (int i = 0; i < tamanho; i++) {
        for (int j = 0; j < tamanho; j++) {
            scanf("%lf", &matriz[i][j]);
        }
    }
}

void printMatriz(double **matriz, int tamanho)
{
    for (int i = 0; i < tamanho; i++) {
        for (int j = 0; j < tamanho; j++) {
            printf("%.*f\t", 3, matriz[i][j]);
        }
        printf("\n");
    }
}

void getCSParameters(double alpha, double beta, double *c, double *s)
{
    if (fabs(alpha) > fabs(beta)) {
        float tau = -beta/alpha;
        *c = 1/sqrt(1+pow(tau, 2));
        *s = *c*tau;
    } else {
        float tau = -alpha/beta;
        *s = 1/sqrt(1+pow(tau, 2));
        *c = *s*tau;
    }
}

void getQRDecomposition(double **matriz, int tamanho, double *c, double *s)
{
    for (int i = 0; i < tamanho-1; i++) {
        getCSParameters(matriz[i][i], matriz[i+1][i], &c[i], &s[i]);
        printf("i=%d c=%f d=%f\n", i, c[i], s[i]);
    }
}

int main()
{
    int tamanho;
    scanf("%d", &tamanho);
    double **matriz = criarMatriz(tamanho);
    lerMatriz(matriz, tamanho);
    //printMatriz(matriz, tamanho);
    double *c = malloc((tamanho-1) * sizeof(double));
    double *s = malloc((tamanho-1) * sizeof(double));
    getQRDecomposition(matriz, tamanho, c, s);

    free(c);
    free(s);
    destruirMatriz(matriz, tamanho);
    return 0;
}