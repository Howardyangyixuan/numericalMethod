#include "myMatrix.h"
int cmp(const void *a, const void *b)
{
    return ((*(double *)a - *(double *)b > 0) ? 1 : -1);
}
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
    //A的备份，用于求特征向量
    Matrix A1;
    Matrix_Init(&A1, n, n);
    MatrixCopy(A, &A1);
    //A的备份，用于双步QR求特征值
    Matrix A2;
    Matrix_Init(&A2, n, n);
    MatrixCopy(A, &A2);
    //1.求拟上三角化：
    Hess(&A2);
    printf("1.A拟上三角化后所得矩阵A(n-1):\n");
    PrintMatrix(A2);
    printf("\n");
    //2.直接QR分解
    printf("2.A(n-1)在第一次QR分解后结果:\n");
    //A(n-1)的备份，用于QR方法
    // Matrix A22;
    // Matrix_Init(&A22, n, n);
    // MatrixCopy(A, &A22);
    Matrix Q;
    Matrix_Init(&Q, n, n);
    Matrix B;
    Matrix_Init(&B, n, n);
    QR_dcp(&A2, &Q);
    printf("A(n-1)在第一次QR分解后所得矩阵R:\n");
    PrintMatrix(A2);
    printf("A(n-1)在第一次QR分解后所得矩阵Q:\n");
    PrintMatrix(Q);
    //3.矩阵A的全部特征值
    printf("A的全部特征值(格式：实部+i*虚部):\n");
    QR2Tran(A, lambda);
    printf("\n");
    //4.矩阵A的全部实特征值的特征向量
    printf("A的全部特征值\n");
    double real[6];
    for (i = 1, j = 0; i <= 10; ++i)
    {
        if (lambda[i].imaginary == 0)
            real[j++] = lambda[i].real;
    }
    //排列
    qsort(real, 6, sizeof(double), cmp);
    for (i = 1; i <= 5; ++i)
    {
        printf("特征值%lf的一个特征向量\n", real[i]);
        IPMTran_eigenvalue(A1, real[i]);
    }
    printf("\n");
    return 0;
}