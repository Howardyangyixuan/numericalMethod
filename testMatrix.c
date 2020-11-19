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
    // PrintMatrix_diag(a);
    // Matrix_diag b;
    // Matrix_diagCopy(a, &b);
    // PrintMatrix_diag(b);
    // b.elem[2][1] = 2;
    // PrintMatrix_diag(b);
    // Trans(&b, 2);
    // PrintMatrix_diag(b);
    Matrix a;
    int m = 10;
    int n = 10;
    Matrix_Init(&a, m, n);
    //    int cnt = 0;
    //    for (int i = 0; i < m; i++)
    //    {
    //        for (int j = 0; j < n; j++)
    //        {
    //            a.elem[i + 1][j + 1] = cnt;
    //            cnt++;
    //        }
    //    }
    //    PrintMatrix(a);
    //    Matrix b;
    //    Matrix_Init(&b, m, n);
    //    cnt = 2;
    //    for (int i = 0; i < m; i++)
    //    {
    //        for (int j = 0; j < n; j++)
    //        {
    //            b.elem[i + 1][j + 1] = cnt;
    //            cnt++;
    //        }
    //    }
    //    PrintMatrix(b);
    //    Matrix c;
    //    Matrix_Init(&c, m, n);
    //    MmM(a, b, &c);
    //    PrintMatrix(c);
    //    M_M(a, b, &c);
    //    PrintMatrix(c);
    //    Vector p;
    //    VectorInit(&p, m);
    //    Vector pp;
    //    VectorInit(&pp, m);
    //    for (int i = 1; i < m; i++)
    //    {
    //        p.elem[i] = 0.5;
    //    }
    //    PrintVector(p);
    //    V_M(p, a, &pp);
    //    PrintVector(pp);
    //    Vector q;
    //    VectorInit(&q, n);
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
    PrintMatrix(a);
    return 0;
}
