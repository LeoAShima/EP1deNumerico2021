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

void setIdentidade(double **matriz, int tamanho)
{
    for (int i = 0; i < tamanho; i++) {
        for (int j = 0; j < tamanho; j++) {
            if (i == j) {
                matriz[i][j] = 1;
            } else {
                matriz[i][j] = 0;
            }
        }
    }
}

void setMatrizTesteA(double **matriz, int tamanho)
{
    for (int i = 0; i < tamanho; i++) {
        for (int j = 0; j < tamanho; j++) {
            if (j == i+1 || i == j+1) {
                matriz[i][j] = -1;
            } else if (i == j) {
                matriz[i][j] = 2;
            } else {
                matriz[i][j] = 0;
            }
        }
    }
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

//multiplicacao de Q*Matriz
void applyGivens(double **matriz, int tamanho, int i, double c, double s)
{
    for (int j = 0; j < tamanho; j++) {
        double matriz1j = matriz[i][j];
        double matriz2j = matriz[i+1][j];
        matriz[i][j] = c*matriz1j-s*matriz2j;
        matriz[i+1][j] = s*matriz1j+c*matriz2j;
    }
}

void applyInverseGivens(double **matriz, int tamanho, int j, double c, double s)
{
    for (int i = 0; i < tamanho; i++) {
        double matrizi1 = matriz[i][j];
        double matrizi2 = matriz[i][j+1];
        matriz[i][j] = c*matrizi1-s*matrizi2;
        matriz[i][j+1] = s*matrizi1+c*matrizi2;
    }
}

void mulRightQ(double **matriz, int tamanho, double *c, double *s)
{
    for (int j = 0; j < tamanho-1; j++) {
        applyInverseGivens(matriz, tamanho, j, c[j], s[j]);
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

void addDiagonal(double **matriz, int tamanho, double MiK)
{
    for (int i = 0; i < tamanho; i++) {
        matriz[i][i] += MiK;
    }
}

int metodoQR(double **matriz, double **V, int tamanho, double precisao, int deslocamento)
{
    double *c = malloc((tamanho-1) * sizeof(double));
    double *s = malloc((tamanho-1) * sizeof(double));
 
    int k = 0;
    double MiK = 0;

    for (int m = tamanho-1; m > 0; m--) {
        do {
            if (k > 0 && deslocamento != 0) {
                MiK = getMiK(matriz, m+1);
            }
            addDiagonal(matriz, m+1, -MiK);
            getQRDecomposition(matriz, m+1, c, s);
            mulRightQ(matriz, m+1, c, s);
            addDiagonal(matriz, m+1, MiK);
            mulRightQ(V, m+1, c, s);
            k++;
        } while (fabs(matriz[m][m-1]) > precisao && k < 9999);
    }

    free(c);
    free(s);
    return k;
}
int main()
{
    int tamanho;
    scanf("%d", &tamanho);
    double **matriz = criarMatriz(tamanho);
    double **V = criarMatriz(tamanho);
    setIdentidade(V, tamanho);
    lerMatriz(matriz, tamanho);
    //printMatriz(matriz, tamanho);
    double *c = malloc((tamanho-1) * sizeof(double));
    double *s = malloc((tamanho-1) * sizeof(double));
    //getQRDecomposition(matriz, tamanho, c, s);
    //printMatriz(matriz, tamanho, 3);

    int k = 0;
    double MiK = 0;

    for (int m = tamanho-1; m > 0; m--) {
        do {
            if (k > 0) {
                MiK = getMiK(matriz, m+1);
            }
            addDiagonal(matriz, m+1, -MiK);
            getQRDecomposition(matriz, m+1, c, s);
            mulRightQ(matriz, m+1, c, s);
            addDiagonal(matriz, m+1, MiK);
            mulRightQ(V, m+1, c, s);
            k++;
        } while (fabs(matriz[m][m-1]) > 0.001 && k < 9999);
        printf("Convergiu m=%d k=%d\n",m,k);
        printMatriz(matriz, tamanho, 3);
    }

    printf("k = %d\n", k);
    printMatriz(matriz, tamanho, 3);
    printMatriz(V, tamanho, 3);

    free(c);
    free(s);
    destruirMatriz(matriz, tamanho);
    destruirMatriz(V, tamanho);
    return 0;
}