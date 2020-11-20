#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myMatrix.h"

int main()
{
    // Matrix_diag a;
    // int n = 5;
    // int s = 2;
    // int r = 2;
    // Matrix_diagInit(&a, n, s, r);
    // //r上三角赋值
    // for (int i = 1; i <= r; ++i)
    // {
    //     for (int j = r - i + 2; j <= n; ++j)
    //     {
    //         a.elem[i][j] = 1;
    //     }
    // }
    // //s下三角赋值
    // for (int i = 1; i <= s; ++i)
    // {
    //     for (int j = 1; j <= n - i; ++j)
    //     {
    //         a.elem[i + r + 1][j] = 3;
    //     }
    // }
    // //1主对角线赋值
    // for (int i = r + 1; i <= r + 1; ++i)
    // {
    //     for (int j = 1; j <= n; ++j)
    //     {
    //         a.elem[i][j] = 2;
    //     }
    // }
    Matrix a, d;
    int m = 2;
    int n = 3;
    Matrix_Init(&a, m, n);
    Matrix_Init(&d, m, n);
    int cnt = 0;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            a.elem[i + 1][j + 1] = cnt;
            cnt++;
        }
    }
    PrintMatrix(a);
    MatrixCopy(a, &d);
    printf("D:\n");
    PrintMatrix(d);
    Matrix b;
    Matrix_Init(&b, m, n);
    cnt = 2;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            b.elem[i + 1][j + 1] = cnt;
            cnt++;
        }
    }
    PrintMatrix(b);
    Matrix c;
    Matrix_Init(&c, m, n);
    MmM(a, b, &c);
    PrintMatrix(c);
    M_M(a, b, &c);
    PrintMatrix(c);
    Vector p;
    VectorInit(&p, m);
    Vector pp;
    VectorInit(&pp, m);
    for (int i = 1; i < m; i++)
    {
        p.elem[i] = 0.5;
    }
    PrintVector(p);
    V_M(p, a, &pp);
    PrintVector(pp);
    Vector q;
    VectorInit(&q, n);
}
