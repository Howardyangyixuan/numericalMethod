#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ElemType double
#define ERROR -2
#define MISMATCH -3  //(维数)不匹配
#define UNORDERED -4 //不能用顺序法求解
#define STRANGE -5   //矩阵奇异
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
        a->elem[(a->s) + 1][i] -= p;
    }
    return;
}

double max(double x, double y)
{
    // return x >= y ? x : y;
    if (x >= y)
        return x;
    else
        return y;
}

double min(double x, double y)
{
    // return x <= y ? x : y;
    if (x <= y)
        return x;
    else
        return y;
}

double max3(double x, double y, double z)
{
    double s;
    s = max(x, y);
    return max(s, z);
}

//三角分解法求带状方阵的行列式
double Diag_LU(Matrix_diag A)
{
    //原矩阵复制到A
    //Matrix_diag A;
    //A = Matrix_diagCopy(M);
    //定义变量
    const int n = A.dimension;
    const int r = A.r;
    const int s = A.s;
    double det, temp;
    int i, j, k, t;
    for (k = 1; k <= n; ++k)
    {
        for (j = k; j <= min(k + s, n); ++j)
        {
            temp = 0;
            for (t = max3(1, k - r, j - s); t <= k - 1; ++t)
            {
                temp += A.elem[k - t + s + 1][t] * A.elem[t - j + s + 1][j];
            }
            A.elem[k - j + s + 1][j] -= temp;
        }
        if (k == n)
            break;
        for (i = k + 1; i <= min(k + r, n); ++i)
        {
            temp = 0;
            for (t = max3(1, i - r, k - s); t <= k - 1; ++t)
            {
                temp += A.elem[i - t + s + 1][t] * A.elem[t - k + s + 1][k];
            }
            A.elem[i - k + s + 1][k] = (A.elem[i - k + s + 1][k] - temp) / A.elem[s + 1][k];
        }
    }
    for (i = 1, det = 1; i <= n; ++i)
    {
        det *= A.elem[s + 1][i];
    }
    //释放所有指针
    // for(i=0;i<=m;++i){
    //     free(A.elem[i]);
    //     A.elem[i] = NULL;
    // }
    // free(A.elem);
    // A.elem = NULL;
    return det;
}

void DiagLU(Matrix_diag *a)
{
    const int n = a->dimension;
    const int r = a->r;
    const int s = a->s;
    int i, j, k, t;
    double temp;
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
    double temp;
    int i, t;
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
    int count;
    double norm, e, lambda_old, lambda_new = 0;
    Vector U;
    VectorInit(&U, n);
    Vector V;
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
        Vector_Num(U, 1 / norm, &V);
        U.elem[1] = a[1] * V.elem[1] + b * V.elem[2] + c * V.elem[3];
        U.elem[2] = b * V.elem[1] + a[2] * V.elem[2] + b * V.elem[3] + c * V.elem[4];
        for (i = 3; i <= 499; ++i)
        {
            U.elem[i] = a[i] * V.elem[i] + b * (V.elem[i - 1] + V.elem[i + 1]) + c * (V.elem[i - 2] + V.elem[i + 2]);
        }
        U.elem[500] = b * V.elem[501] + a[500] * V.elem[500] + b * V.elem[499] + c * V.elem[498];
        U.elem[501] = a[501] * V.elem[501] + b * V.elem[500] + c * V.elem[499];
        lambda_old = lambda_new;
        lambda_new = DotProduct(V, U);
        e = fabs((lambda_new - lambda_old) / lambda_new); //求相对误差
        ++count;
    }
    free(V.elem);
    free(U.elem);
    return lambda_new;
}

double PMTran_eigenvalue(Matrix_diag A, double p)
{
    double b = 0.16;
    double c = -0.064;
    double a[502];
    int i, j;
    for (j = 1; j <= 501; ++j)
    {
        a[j] = (1.64 - 0.024 * j) * sin(0.2 * j) - 0.64 * exp(0.1 / j) - p;
    }
    //定义变量
    const int n = A.dimension;
    int count;
    double norm, e, lambda_old, lambda_new = 0;
    Vector U;
    VectorInit(&U, n);
    Vector V;
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
        Vector_Num(U, 1 / norm, &V);
        U.elem[1] = a[1] * V.elem[1] + b * V.elem[2] + c * V.elem[3];
        U.elem[2] = b * V.elem[1] + a[2] * V.elem[2] + b * V.elem[3] + c * V.elem[4];
        for (i = 3; i <= 499; ++i)
        {
            U.elem[i] = a[i] * V.elem[i] + b * (V.elem[i - 1] + V.elem[i + 1]) + c * (V.elem[i - 2] + V.elem[i + 2]);
        }
        U.elem[500] = b * V.elem[501] + a[500] * V.elem[500] + b * V.elem[499] + c * V.elem[498];
        U.elem[501] = a[501] * V.elem[501] + b * V.elem[500] + c * V.elem[499];
        lambda_old = lambda_new;
        lambda_new = DotProduct(V, U);
        e = fabs((lambda_new - lambda_old) / (lambda_new + p)); //求相对误差
        ++count;
    }
    free(V.elem);
    free(U.elem);
    return lambda_new + p;
}

//反幂法
double IPM_eigenvalue(Matrix_diag A)
{
    int const n = A.dimension;
    const int r = A.r;
    const int s = A.s;
    int i, count;
    double look, norm, e, lambda_old, lambda_new = 1;
    Vector U;
    VectorInit(&U, n);
    Vector V;
    VectorInit(&V, n);
    for (i = 1; i <= n; ++i)
    {
        U.elem[i] = 1;
    }
    e = 1;
    count = 0;
    Matrix_diag A1;
    Matrix_diagCopy(A, &A1);
    DiagLU(&A1);
    do
    {
        norm = sqrt(DotProduct(U, U));
        Vector_Num(U, 1 / norm, &V);
        Vector tmp;
        VectorInit(&tmp, n);
        VectorCopy(V, &tmp);
        solveDiagLU(&A1, tmp, &U);
        lambda_old = lambda_new;
        lambda_new = DotProduct(V, U);
        e = fabs((lambda_new - lambda_old) / lambda_old); //求相对误差
        ++count;
        look = 1 / lambda_new;
    } while (e > E);
    free(V.elem);
    free(U.elem);
    free(A1.elem[1]);
    free(A1.elem[2]);
    free(A1.elem[3]);
    free(A1.elem);
    return 1 / lambda_new;
}

double IPMTran_eigenvalue(Matrix_diag A, double p)
{
    double l;
    Trans(&A, p);
    l = IPM_eigenvalue(A);
    Trans(&A, (-1) * p);
    return l + p;
}

int main()
{
    const int n = 501;
    const int s = 2;
    const int r = 2;
    double b = 0.16;
    double c = -0.064;
    double a[502];
    int i, j;
    for (j = 1; j <= n; ++j)
    {
        a[j] = (1.64 - 0.024 * j) * sin(0.2 * j) - 0.64 * exp(0.1 / j);
    }
    double lambda[10];
    Matrix_diag A;
    Matrix_diagInit(&A, n, s, r);
    double det;
    i = 1;
    for (j = 1; j <= n; ++j)
    {
        A.elem[i][j] = c;
    }
    i = 2;
    for (j = 1; j <= n; ++j)
    {
        A.elem[i][j] = b;
    }
    i = 3;
    for (j = 1; j <= n; ++j)
    {
        A.elem[i][j] = a[j];
    }
    i = 4;
    for (j = 1; j <= n; ++j)
    {
        A.elem[i][j] = b;
    }
    i = 5;
    for (j = 1; j <= n; ++j)
    {
        A.elem[i][j] = c;
    }
    Matrix_diag A1;
    Matrix_diagCopy(A, &A1);

    lambda[0] = IPM_eigenvalue(A);
    lambda[1] = PM_eigenvalue(A);
    lambda[2] = PMTran_eigenvalue(A, lambda[1]);
    lambda[3] = (lambda[1] < lambda[2]) ? lambda[1] : lambda[2];
    lambda[4] = (lambda[1] > lambda[2]) ? lambda[1] : lambda[2];
    printf("2020.11.7 杨义轩ZY2006330 数值分析上机实习一\n");
    printf("lambda_1 = %.12e\n", lambda[3]);
    printf("lambda_501 = %.12e\n", lambda[4]);
    printf("lambda_s = %.12e\n", lambda[0]);

    double mu[40];
    for (i = 1; i <= 39; ++i)
    {
        mu[i] = lambda[3] + i / 40.0 * (lambda[4] - lambda[3]);
        mu[i] = IPMTran_eigenvalue(A, mu[i]);
        printf("Lambda_i_%d = %.12e\n", i, mu[i]);
    }
    printf("cond(A)2 = %.12e\n", fabs(lambda[1] / lambda[0]));

    det = Diag_LU(A1);
    printf("det(A) = %.12e\n", det);
    return 0;
}
