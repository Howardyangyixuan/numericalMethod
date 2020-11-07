#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ElemType double
#define ERROR -2
#define MISMATCH -3  //(维数)不匹配
#define UNORDERED -4 //不能用顺序法求解
#define STRANGE 5    //矩阵奇异
#define E 1.0e-12    //相对误差限

typedef struct
{
    int dimension;  //维数
    ElemType *elem; //向量,一维数组
} Vector;

typedef struct
{
    int dimension; //矩阵维数
    int s;         //上半带宽
    int r;         //下半带宽
    double **elem; //矩阵，二维数组
} Matrix_diag;     //n维s+1+r对角方阵

//输入&V，n；初始化n维向量V（用1~n):
void VectorInit(Vector *v, int n)
{
    v->dimension = n;

    v->elem = (ElemType *)malloc(sizeof(ElemType) * (n + 1));
    return;
}

void VectorCopy(Vector V, Vector *x)
{
    if (V.dimension != x->dimension)
        exit(MISMATCH);
    for (int i = 1; i <= V.dimension; ++i)
    {
        (x->elem)[i] = V.elem[i];
    }
    return;
}
void PrintVector(Vector V)
{
    for (int i = 1; i <= V.dimension; ++i)
    {
        printf("%lf ", V.elem[i]);
    }
    printf("\n");
    return;
}

//V = 向量U * 数k(V已初始化)
Vector *Vector_Num(Vector U, double k, Vector *v)
{
    if (U.dimension == v->dimension)
    {
        for (int i = 1; i <= U.dimension; ++i)
        {
            v->elem[i] = U.elem[i] * k;
        }
    }
    else
    {
        exit(MISMATCH);
    }
    return v;
}

//求向量内积
double DotProduct(Vector X, Vector Y)
{
    if (X.dimension == Y.dimension)
    {
        double p = 0;
        for (int i = 1; i <= X.dimension; ++i)
        {
            p += X.elem[i] * Y.elem[i];
        }
        return p;
    }
    else
    {
        exit(MISMATCH);
    }
}

//输入&M，n：初始化n维s+1+r对角方阵M(用(1~s+1+r)*(1~n)):
void Matrix_diagInit(Matrix_diag *a, int n, int s, int r)
{
    int m = s + 1 + r;
    a->dimension = n;
    a->s = s;
    a->r = r;
    a->elem = (double **)malloc(sizeof(double *) * (m + 1));
    for (int i = 1; i <= m; ++i)
    {
        (a->elem)[i] = (double *)malloc(sizeof(double) * (n + 1));
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

double max(double x, double y)
{
    return x > y ? x : y;
}

double min(double x, double y)
{
    return x < y ? x : y;
}

double max3(double x, double y, double z)
{
    double s;
    s = max(x, y);
    return max(s, z);
}
void DiagLU(Matrix_diag *a)
{
    const int n = a->dimension;
    const int r = a->r;
    const int s = a->s;
    const int m = s + 1 + r;
    double det, temp;
    int i, j, k, t;
    for (k = 1; k <= n; ++k)
    {
        for (j = k; j <= min(k + s, n); ++j)
        {
            for (temp = 0, t = max3(1, k - r, j - s); t <= k - 1; ++t)
            {
                temp += a->elem[k - t + s + 1][t] * a->elem[t - j + s + 1][j];
            }
            a->elem[k - j + s + 1][j] -= temp;
        }
        if (k == n)
            break;
        for (i = k + 1; i <= min(k + r, n); ++i)
        {
            for (temp = 0, t = max3(1, i - r, k - s); t <= k - 1; ++t)
            {
                temp += a->elem[i - t + s + 1][t] * a->elem[t - k + s + 1][k];
            }
            a->elem[i - k + s + 1][k] = (a->elem[i - k + s + 1][k] - temp) / a->elem[s + 1][k];
        }
    }
    return;
}

//带状方程组AX=B三角分解之后解方程(X需初始化)
void solveDiagLU(Matrix_diag *a, Vector B, Vector *x)
{
    const int n = a->dimension;
    const int r = a->r;
    const int s = a->s;
    const int m = s + 1 + r;
    double det, temp;
    int i, j, k, t;
    //解方程
    for (i = 2; i <= n; ++i)
    {
        for (t = max(1, i - r), temp = 0; t <= i - 1; ++t)
        {
            temp += a->elem[i - t + s + 1][t] * B.elem[t];
        }
        B.elem[i] -= temp;
    }
    x->elem[n] = B.elem[n] / a->elem[s + 1][n];
    for (i = n - 1; i >= 1; --i)
    {
        for (t = i + 1, temp = 0; t <= min(i + s, n); ++t)
        {
            temp += a->elem[i - t + s + 1][t] * x->elem[t];
        }
        x->elem[i] = (B.elem[i] - temp) / a->elem[s + 1][i];
    }
    return;
}

double IPM_eigenvalue(Matrix_diag A)
{
    const int n = A.dimension;
    const int r = A.r;
    const int s = A.s;
    const int m = s + 1 + r;
    int i, count;
    double look, norm, e, lambda_old, lambda_new = 1;
    Vector U;
    Vector *u = &U;
    VectorInit(&U, n);
    Vector V;
    Vector *v = &V;
    VectorInit(&V, n);
    for (i = 1; i <= n; ++i)
    {
        U.elem[i] = 1;
        V.elem[i] = 1;
    }
    // for (i = 200; i <= n; ++i)
    // {
    //     U.elem[i] = 0;
    // }
    e = 1;
    count = 0;
    Matrix_diag A1;
    Matrix_diagCopy(A, &A1);
    DiagLU(&A1);
    for (int j = 0; j < 2; j++)
    {
        printf("--------------------------------------------------------------------------\n");
        printf("U: ");
        PrintVector(U);
        printf("V: ");
        PrintVector(V);
        norm = sqrt(DotProduct(U, U));
        v = Vector_Num(U, 1 / norm, &V);
        printf("norm %.12e\n", norm);
        printf("after U: ");
        PrintVector(U);
        printf("after V: ");
        PrintVector(V);
        Vector tmp;
        VectorInit(&tmp, 3);
        VectorCopy(V, &tmp);
        solveDiagLU(&A1, tmp, &U);
        printf("after U: ");
        PrintVector(U);
        printf("after V: ");
        PrintVector(V);
        lambda_old = lambda_new;
        lambda_new = DotProduct(V, U);
        e = fabs((lambda_new - lambda_old) / lambda_new); //求相对误差
        ++count;
        printf("epoch : %d, error : %.12e\n", count, e);
        look = 1 / lambda_new;
    }
    free(V.elem);
    free(U.elem);
    free(A1.elem[1]);
    free(A1.elem[2]);
    free(A1.elem[3]);
    free(A1.elem);
    return 1 / lambda_new;
}
int main()
{
    // const int n = 501;
    // const int s = 2;
    // const int r = 2;
    // double b = 0.16;
    // double c = -0.064;
    // double a[502];
    // int i, j;
    // for (j = 1; j <= n; ++j)
    // {
    //     a[j] = (1.64 - 0.024 * j) * sin(0.2 * j) - 0.64 * exp(0.1 / j);
    // }
    // double lambda[10];
    // Matrix_diag A;
    // Matrix_diagInit(&A, n, s, r);
    // i = 1;
    // for (j = 1; j <= n; ++j)
    // {
    //     A.elem[i][j] = c;
    // }
    // i = 2;
    // for (j = 1; j <= n; ++j)
    // {
    //     A.elem[i][j] = b;
    // }
    // i = 3;
    // for (j = 1; j <= n; ++j)
    // {
    //     A.elem[i][j] = a[j];
    // }
    // i = 4;
    // for (j = 1; j <= n; ++j)
    // {
    //     A.elem[i][j] = b;
    // }
    // i = 5;
    // for (j = 1; j <= n; ++j)
    // {
    //     A.elem[i][j] = c;
    // }
    // Matrix_diag a;
    // int n = 5;
    // int s = 2;
    // int r = 2;
    // Matrix_diagInit(&a, n, s, r);
    // //r上三角赋值
    // for (int i = 1; i <= n; ++i)
    // {
    //     for (int j = 1; j <= n; ++j)
    //     {
    //         a.elem[i][j] = 1;
    //     }
    // }
    // //s下三角赋值
    // for (int i = 1; i <= s; ++i)
    // {
    //     for (int j = 1; j <= n; ++j)
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
    Matrix_diag a;
    int n = 3;
    int s = 1;
    int r = 1;
    Matrix_diagInit(&a, n, s, r);
    //r上三角赋值
    for (int i = 1; i <= r; ++i)
    {
        for (int j = r - i + 2; j <= n; ++j)
        {
            a.elem[i][j] = 1;
        }
    }
    //s下三角赋值
    for (int i = 1; i <= s; ++i)
    {
        for (int j = 1; j <= n - i; ++j)
        {
            a.elem[i + r + 1][j] = 2;
        }
    }
    //1主对角线赋值
    for (int i = r + 1; i <= r + 1; ++i)
    {
        for (int j = 1; j <= n; ++j)
        {
            a.elem[i][j] = 3;
        }
    }
    // PrintMatrix_diag(a);
    // DiagLU(&a);
    PrintMatrix_diag(a);
    double tmp = IPM_eigenvalue(a);
    // double tmp = IPM_eigenvalue(A);
    printf("IPM result: %.12e\n", tmp);
    // printf("PMTrans result: %.12e\n", IPMTran_eigenvalue(a, tmp));
    // PrintMatrix_diag(A1);
    return 0;
}