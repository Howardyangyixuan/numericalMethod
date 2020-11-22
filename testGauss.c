#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myMatrix.h"

int main()
{
    Matrix a;
    int m = 3;
    int n = 3;
    Matrix_Init(&a, m, n);
    Vector x;
    VectorInit(&x, m);
    Vector b;
    VectorInit(&b, m);
    a.elem[1][1] = 2;
    a.elem[1][2] = 1;
    a.elem[1][3] = 2;
    a.elem[2][1] = 5;
    a.elem[2][2] = -1;
    a.elem[2][3] = 1;
    a.elem[3][1] = 1;
    a.elem[3][2] = -3;
    a.elem[3][3] = -4;
    b.elem[1] = 5;
    b.elem[2] = 8;
    b.elem[3] = -4;
    PrintMatrix(a);
    GaussElimination(a, b, &x);
    // PrintMatrix(a);
    // PrintVector(x);
    return 0;
}
