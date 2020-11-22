#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myMatrix.h"
int main()
{
    //题目常量
    const int n = 10;
    int i, j;
    double times;

    Matrix A;
    Matrix_Init(&A, n, n);
    complex lambda[11];
    //用第一个复数的部分表示表中下一个存储空间
    lambda[0].real = 1;

    //初始赋值
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            if (i != j)
            {
                A.elem[i][j] = sin(0.5 * i + 0.2 * j);
            }
            else
            {
                A.elem[i][j] = 1.52 * cos(i + 1.2 * j);
            }
        }
    }
    printf("A的全部特征值\n");
    QR2Tran(A, lambda);
    // //A的备份，用于求特征向量
    // Matrix A1;
    // Matrix_Init(&A1, n, n);
    // MatrixCopy(A, &A1);
    // //A的备份，用于双步QR求特征值
    // Matrix A2;
    // Matrix_Init(&A2, n, n);
    // MatrixCopy(A, &A2);
}