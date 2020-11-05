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

//三角分解法求带状方阵的行列式
double DiagLU_det(Matrix_diag A)
{
    //原矩阵复制到A
    //Matrix_diag A;
    //A = Matrix_diagCopy(M);
    //定义变量
    const int n = A.dimension;
    const int r = A.r;
    const int s = A.s;
    const int m = s + 1 + r;
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
    PrintMatrix_diag(A);
    for (i = 1, det = 1; i <= n; ++i)
    {
        det *= A.elem[i][i];
        printf("epoch: %d, det: %.12lf\n", i, det);
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
int main()
{
    Matrix_diag a;
    int n = 5;
    int s = 2;
    int r = 2;
    Matrix_diagInit(&a, n, s, r);
    //r上三角赋值
    for (int i = 1; i <= n; ++i)
    {
        for (int j = 1; j <= n; ++j)
        {
            a.elem[i][j] = 1;
        }
    }
    //s下三角赋值
    for (int i = 1; i <= s; ++i)
    {
        for (int j = 1; j <= n; ++j)
        {
            a.elem[i + r + 1][j] = 3;
        }
    }
    //1主对角线赋值
    for (int i = r + 1; i <= r + 1; ++i)
    {
        for (int j = 1; j <= n; ++j)
        {
            a.elem[i][j] = 2;
        }
    }
    PrintMatrix_diag(a);
    // DiagLU(&a);
    // PrintMatrix_diag(a);
    // Vector b;
    // VectorInit(&b, 5);
    // b.elem[1] = 1;
    // b.elem[2] = 1;
    // b.elem[3] = 2;
    // b.elem[4] = 3;
    // b.elem[5] = 3;
    // Vector ans;
    // VectorInit(&ans, 5);
    // solveDiagLU(&a, b, &ans);
    // PrintVector(ans);
    double tmp = DiagLU_det(a);
    printf("det: %.12e", tmp);
    return 0;
}