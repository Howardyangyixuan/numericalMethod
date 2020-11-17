#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ElemType double
#define ERROR -2
#define MISMATCH -3  //(维数)不匹配
#define UNORDERED -4 //不能用顺序法求解
#define STRANGE 5    //矩阵奇异
#define E 1.0e-2     //相对误差限

typedef struct
{
    int dimension; //矩阵维数
    int s;         //上半带宽
    int r;         //下半带宽
    double **elem; //矩阵，二维数组
} Matrix_diag;     //n维s+1+r对角方阵

typedef struct
{
    int m;
    int n;
    double **elem;
} Matrix;

//输入&M，n：初始化m*n维矩阵:
void Matrix_Init(Matrix *a, int m, int n)
{
    if (m < 2 | n < 2)
    {
        printf("矩阵错误：m或n小于零");
        exit(MISMATCH);
    }
    a->m = m;
    a->n = n;
    a->elem = (double **)malloc(sizeof(double *) * (m + 1));
    for (int i = 1; i <= m; ++i)
    {
        (a->elem)[i] = (double *)malloc(sizeof(double) * (n + 1));
        //for mac non-zero init
        for (int j = 1; j <= n; ++j)
        {
            a->elem[i][j] = 0;
        }
    }
    return;
}

//输出矩阵A:
void PrintMatrix(Matrix A)
{
    for (int i = 1; i <= A.m; ++i)
    {
        for (int j = 1; j <= A.n; ++j)
        {
            printf("%lf ", A.elem[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    return;
}

//输入&M，n：初始化n维s+1+r对角方阵M(用(1~s+1+r)*(1~n)):
void Matrix_diagInit(Matrix_diag *a, int n, int s, int r)
{
    int m = s + 1 + r;
    if (m > n)
    {
        printf("对角方阵错误：n+1+r大于n");
        exit(MISMATCH);
    }
    a->dimension = n;
    a->s = s;
    a->r = r;
    a->elem = (double **)malloc(sizeof(double *) * (m + 1));
    for (int i = 1; i <= m; ++i)
    {
        (a->elem)[i] = (double *)malloc(sizeof(double) * (n + 1));
        //for mac non-zero init
        for (int j = 1; j <= n; ++j)
        {
            a->elem[i][j] = 0;
        }
    }
    return;
}

//输出带状方阵A:
void PrintMatrix_diag(Matrix_diag A)
{
    for (int i = 1; i <= A.s + 1 + A.r; ++i)
    {
        for (int j = 1; j <= A.dimension; ++j)
        {
            printf("%lf ", A.elem[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    return;
}

//把A值复制给A1（A为带状阵，A1不需要初始化）
void Matrix_diagCopy(Matrix_diag A, Matrix_diag *a1)
{
    Matrix_diagInit(a1, A.dimension, A.s, A.r);
    for (int i = 1; i <= A.s + A.r + 1; ++i)
    {
        for (int j = 1; j <= A.dimension; ++j)
        {
            a1->elem[i][j] = A.elem[i][j];
        }
    }
    return;
}

//矩阵平移(A->A-pI)
void Trans(Matrix_diag *a, double p)
{
    for (int i = 1; i <= a->dimension; ++i)
    {
        a->elem[a->s + 1][i] -= p;
    }
    return;
}

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
    int m = 3;
    int n = 2;
    Matrix_Init(&a, m, n);
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

    return 0;
}