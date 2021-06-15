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

void printMatriz(double **matriz, int tamanho, int precisao)
{
    for (int i = 0; i < tamanho; i++) {
        for (int j = 0; j < tamanho; j++) {
            printf("%.*f\t", precisao, matriz[i][j]);
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

void applyGivens(double **matriz, int tamanho, int i, double c, double s)
{
    for (int j = 0; j < tamanho; j++) {
        double matriz1j = matriz[i][j];
        double matriz2j = matriz[i+1][j];
        matriz[i][j] = c*matriz1j-s*matriz2j;
        matriz[i+1][j] = s*matriz1j+c*matriz2j;
    }
}

void getQRDecomposition(double **matriz, int tamanho, double *c, double *s)
{
    for (int i = 0; i < tamanho-1; i++) {
        getCSParameters(matriz[i][i], matriz[i+1][i], &c[i], &s[i]);
        applyGivens(matriz, tamanho, i, c[i], s[i]);
        //printf("i=%d c=%f d=%f\n", i, c[i], s[i]);
    }
}

double getMiK(double **matriz, int tamanho)
{
    double dk = (matriz[tamanho-2][tamanho-2] - matriz[tamanho-1][tamanho-1])/2;
    if (dk >= 0) {
        return matriz[tamanho-1][tamanho-1] + dk - sqrt(pow(dk,2) + pow(matriz[tamanho-1][tamanho-2],2));
    } else {
        return matriz[tamanho-1][tamanho-1] + dk + sqrt(pow(dk,2) + pow(matriz[tamanho-1][tamanho-2],2));
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
    printMatriz(matriz, tamanho, 3);

    free(c);
    free(s);
    destruirMatriz(matriz, tamanho);
    return 0;
}