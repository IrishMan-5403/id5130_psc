#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define TYPE float
#define N 10
#define SMALLVALUE 0.001

void initializeMatrix(TYPE mat[][N])
{
#pragma acc parallel loop collapse(2)
    for (int ii = 0; ii < N; ++ii)
    {
        for (int jj = 0; jj < N && jj < ii; ++jj)
        {
            mat[ii][jj] = (ii + jj) / (float)N / N;
            mat[jj][ii] = (ii + jj) / (float)N / N;
        }
    }

#pragma acc parallel loop
    for (int ii = 0; ii < N; ++ii)
        mat[ii][ii] = 1.0;
}

void printMatrix(TYPE a[][N])
{
#pragma acc parallel loop
    for (int ii = 0; ii < N; ++ii)
    {
        for (int jj = 0; jj < N; ++jj)
            printf("%.2f ", a[ii][jj]);
        printf("\n");
    }
}

void choleskyDecomposition(TYPE a[][N])
{
#pragma acc parallel loop collapse(2)
    for (int ii = 0; ii < N; ++ii)
    {
        for (int jj = 0; jj < ii; ++jj)
        {
            for (int kk = 0; kk < jj; ++kk)
            {
                a[ii][jj] += -a[ii][kk] * a[jj][kk];
            }
            a[ii][jj] /= (a[jj][jj] > SMALLVALUE ? a[jj][jj] : 1);
        }
#pragma acc parallel loop
        for (int kk = 0; kk < ii; ++kk)
        {
            a[ii][ii] += -a[ii][kk] * a[ii][kk];
        }
        a[ii][ii] = sqrt(a[ii][ii]);
    }
}

int main()
{
    TYPE matrix[N][N];
    clock_t timer;
    timer = clock();

#pragma acc data copy(matrix)
    {
        initializeMatrix(matrix);
        choleskyDecomposition(matrix);
#pragma acc wait
        printMatrix(matrix);
    }

    timer = clock() - timer;
    printf("The operation took %ld clicks (%f seconds)\n", timer, ((float)timer) / CLOCKS_PER_SEC);

    return 0;
}