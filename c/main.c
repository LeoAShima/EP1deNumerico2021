#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double** criarMatriz(int tamanho)
{
    double **matriz = malloc(tamanho * sizeof(*matriz));
    for (int i = 0; i < tamanho; i++) {
        matriz[i] = malloc(tamanho * sizeof(*matriz[i]));
    }
    for (int i = 0; i < tamanho; i++) {
        for (int j = 0; j < tamanho; j++) {
            matriz[i][j] = 0;
        }
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
            printf("% .*f\t", precisao, matriz[i][j]);
        }
        printf("\n");
    }
}

void printAutovalores(double **matriz, int tamanho, int precisao)
{
    for (int i = 0; i < tamanho; i++) {
        printf("%.*f", precisao, matriz[i][i]);
        if (i < tamanho-1) printf(" / ");
    }
    printf("\n");
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
            mulRightQ(V, tamanho, c, s);
            k++;
        } while (fabs(matriz[m][m-1]) > precisao && k < 9999);
    }

    free(c);
    free(s);
    return k;
}

void tarefaA()
{
    printf("Digite o tamanho da matriz a ser testada: ");
    int tamanho;
    scanf("%d", &tamanho);
    double **matriz = criarMatriz(tamanho);
    double **V = criarMatriz(tamanho);
    setIdentidade(V, tamanho);
    setMatrizTesteA(matriz, tamanho);

    int qtdIteracoes = metodoQR(matriz, V, tamanho, 0.000001, 0);
    printf("\nMetodo QR sem deslocamento espectral\n");
    printf("Convergiu em %d iteracoes\n", qtdIteracoes);
    printf("Autovalores:\n");
    printAutovalores(matriz, tamanho, 3);
    printf("Autovetores:\n");
    printMatriz(V, tamanho, 3);
    
    setIdentidade(V, tamanho);
    setMatrizTesteA(matriz, tamanho);
    qtdIteracoes = metodoQR(matriz, V, tamanho, 0.000001, 1);
    printf("\nMetodo QR com deslocamento espectral\n");
    printf("Convergiu em %d iteracoes\n", qtdIteracoes);
    printf("Autovalores:\n");
    printAutovalores(matriz, tamanho, 3);
    printf("Autovetores:\n");
    printMatriz(V, tamanho, 3);
    
    destruirMatriz(matriz, tamanho);
    destruirMatriz(V, tamanho);
}

void escreverMatrizA(double **A, int tamanho, double massa, double *k)
{
    int i = 0;
    for (i = 0; i < tamanho-1; i++) {
        A[i][i] = (k[i] + k[i+1])/massa;
        A[i][i+1] = (-1*k[i+1])/massa;
        A[i+1][i] = (-1*k[i+1])/massa;
    }
    A[i][i] = (k[i] + k[i+1])/massa;
}

// Realiza a multiplicacao Q^T * X
void QTX(double **Q, int tamanho, double *x, double *y)
{
    //printMatriz(Q, tamanho, 3);
    for (int i = 0; i < tamanho; i++) {
        y[i] = 0;
        for (int j = 0; j < tamanho; j++) {
            y[i] += Q[j][i] * x[j];
            //printf("Qji=%f, xi=%f\n",Q[j][i],x[i]);
            //printf("i=%d, xi=%f\n",i,x[i]);
        }
        //printf("i=%d, xi=%f\n",i,x[i]);
        printf("i=%d, yi=%f\n",i,y[i]);
    }
}

void tarefaB()
{
    int numMassas = 5;
    double massa = 2;
    double *x = malloc((numMassas) * sizeof(double));
    double *y = malloc((numMassas) * sizeof(double));

    //TODO: Ler valores de massa
    x[0] = -2;
    x[1] = -3;
    x[2] = -1;
    x[3] = -3;
    x[4] = -1;

    double *k = malloc((numMassas+1) * sizeof(double));
    for (int i = 0; i < numMassas+1; i++) {
        k[i] = 40 + 2*i;
    }

    double **A = criarMatriz(numMassas);
    double **Q = criarMatriz(numMassas);
    escreverMatrizA(A, numMassas, massa, k);
    setIdentidade(Q, numMassas);
    metodoQR(A, Q, numMassas, 0.000001, 1); //Apos chamada, A se torna LAMBDA
    //printMatriz(Q, numMassas, 3);
    QTX(Q, numMassas, x, y); //obtem y0 a partir de x0, Q contem autovetores

    //debug start
    for (int i = 0; i < numMassas; i++) { //Para cada x
            double xres = 0; //Contem valor de x_i no instante t;
            for (int l = 0; l < numMassas; l++) {
                xres += Q[i][l] * y[l];
            }
            printf("x_%d(0) = %.*f\n",i,3,xres);
    }
    //debug end 


    for (int i = 0; i < numMassas; i++) { //Para cada x
        printf("x_%d = [", i);
        for (int j = 0; j < 401; j++) { //Imprimir 400 valores de x(t)
            double xres = 0; //Contem valor de x_i no instante t;
            for (int l = 0; l < numMassas; l++) {
                xres += Q[i][l] * y[l]*cos(sqrt(A[l][l])*j*0.025);
            }
            printf("%.*f ",3,xres);
        }
        printf("]\n");
    }


    destruirMatriz(A, numMassas);
    destruirMatriz(Q, numMassas);
    free(x);
    free(y);
}

int main()
{
    //tarefaA();
    tarefaB();
    return 0;
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