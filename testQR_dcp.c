#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myMatrix.h"

int main()
{
    Matrix a;
    int m = 10;
    int n = 10;
    Matrix_Init(&a, m, n);
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            if (i != j)
            {
                a.elem[i][j] = sin(0.5 * i + 0.2 * j);
            }
            else
            {
                a.elem[i][j] = 1.52 * cos(i + 1.2 * j);
            }
        }
    }
    PrintMatrix(a);
    Hess(&a);
    Matrix q;
    Matrix_Init(&q, m, n);
    QR_dcp(&a, &q);
    PrintMatrix(a);
    PrintMatrix(q);
    return 0;
}
