#include "myMatrix.h"
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

//幂法(求矩阵按模最大的特征值):
double PM_eigenvalue(Matrix_diag A)
{
    //本题常量
    double b = 0.16;
    double c = -0.064;
    double a[502];
    int i, j;
    for (j = 1; j <= 501; ++j)
    {
        a[j] = (1.64 - 0.024 * j) * sin(0.2 * j) - 0.64 * exp(0.1 / j);
    }
    const int n = A.dimension;
    const int r = A.r;
    const int s = A.s;
    const int m = s + 1 + r;
    int count, k;
    double norm, e, lambda_old, lambda_new = 0;
    Vector U;
    Vector *u = &U;
    VectorInit(&U, n);
    Vector V;
    Vector *v = &V;
    VectorInit(&V, n);
    for (i = 1; i <= n; ++i)
    {
        U.elem[i] = 1;
    }
    e = 1;
    count = 0;
    while (e > E)
    {
        norm = sqrt(DotProduct(U, U));
        v = Vector_Num(U, 1 / norm, &V);
        //等价于Vector_Num(U, 1 / norm, &V);
        U.elem[1] = a[1] * V.elem[1] + b * V.elem[2] + c * V.elem[3];
        U.elem[2] = b + V.elem[1] + a[2] * V.elem[2] + b * V.elem[3] + c * V.elem[4];
        for (i = 3; i <= 499; ++i)
        {
            U.elem[i] = a[i] * V.elem[i] + b * (V.elem[i - 1] + V.elem[i + 1]) + c * (V.elem[i - 2] + V.elem[i + 2]);
        }
        U.elem[500] = b * V.elem[501] + a[500] * V.elem[500] + b * V.elem[499] + c * V.elem[498];
        U.elem[501] = a[501] * V.elem[501] + b * V.elem[500] + c * V.elem[499];
        lambda_old = lambda_new;
        lambda_new = DotProduct(V, U);
        e = fabs(lambda_new - lambda_old) / fabs(lambda_new); //求相对误差
        ++count;
        printf("epoch : %d, error : %.12lf\n", count, e);
    }
    free(V.elem);
    free(U.elem);
    return lambda_new;
}

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
//矩阵 * 向量 A * x = b
void M_V(Matrix A, Vector x, Vector *b)
{
    if (A.n != x.dimension || x.dimension != b->dimension)
    {
        printf("维数不匹配");
        exit(MISMATCH);
    }
    for (int i = 1; i <= A.m; i++)
    {
        double sum = 0;
        for (int j = 1; j <= A.n; j++)
        {
            sum += A.elem[i][j] * x.elem[j];
        }
        b->elem[i] = sum;
    }
}
//向量 * 矩阵 x * A = b
void V_M(Vector x, Matrix A, Vector *b)
{
    if (A.m != x.dimension || A.n != b->dimension)
    {
        printf("维数不匹配");
        exit(MISMATCH);
    }
    for (int i = 1; i <= A.n; i++)
    {
        double sum = 0;
        for (int j = 1; j <= A.m; j++)
        {
            sum += A.elem[i][j] * x.elem[j];
        }
        b->elem[i] = sum;
    }
}
//矩阵 * 矩阵 A * B = C
void M_M(Matrix A, Matrix B, Matrix *c)
{
    if (A.n != B.m || A.m != c->m || B.n != c->n)
    {
        printf("维数不匹配");
        exit(MISMATCH);
    }
    for (int i = 1; i <= A.m; i++)
    {
        for (int j = 1; j <= A.n; j++)
        {
            double sum = 0;
            for (int k = 1; k <= A.n; k++)
            {
                sum += A.elem[i][k] * B.elem[k][j];
            }
            c->elem[i][j] = sum;
        }
    }
}
//矩阵 - 矩阵 A - B = C
void MmM(Matrix A, Matrix B, Matrix *c)
{
    if (A.m != B.m || B.n != B.n || c->m != A.m || c->n != A.n)
    {
        printf("维数不匹配");
        exit(MISMATCH);
    }
    for (int i = 1; i <= A.m; i++)
    {
        for (int j = 1; j <= A.n; j++)
        {
            c->elem[i][j] = A.elem[i][j] - B.elem[i][j];
        }
    }
}