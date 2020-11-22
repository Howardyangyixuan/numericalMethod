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
            // printf("%.12e& ", A.elem[i][j]);
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
            sum += A.elem[j][i] * x.elem[j];
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
//向量 * 向量 a * bT = C
Matrix *V_V(Vector a, Vector b, Matrix *c)
{
    if (a.dimension != b.dimension || c->m != a.dimension || c->n != a.dimension)
    {
        printf("维数不匹配");
        exit(MISMATCH);
    }
    for (int i = 1; i <= c->m; i++)
    {
        for (int j = 1; j <= c->n; j++)
        {
            c->elem[i][j] = a.elem[i] * b.elem[j];
        }
    }
    return c;
}
//矩阵拷贝
void MatrixCopy(Matrix A, Matrix *a)
{
    if (a->m != A.m || a->n != A.n)
    {
        printf("维数不匹配");
        exit(MISMATCH);
    }
    for (int i = 1; i <= A.m; i++)
    {
        for (int j = 1; j <= A.n; j++)
        {
            a->elem[i][j] = A.elem[i][j];
        }
    }
}
//拟上三角化
void Hess(Matrix *a)
{
    if (a->m != a->n)
    {
        printf("a不是方阵");
        exit(MISMATCH);
    }
    const int n = a->m;
    int i, j, r;
    double d, temp, c, h;
    Vector U;
    VectorInit(&U, n);

    Vector V;
    VectorInit(&V, n);
    Vector P;
    VectorInit(&P, n);
    Vector Q;
    VectorInit(&Q, n);

    Vector Vtemp;
    VectorInit(&Vtemp, n);
    Matrix Mtemp;
    Matrix_Init(&Mtemp, n, n);
    for (r = 1; r <= n - 2; ++r)
    {
        //求c,h
        for (d = 0, i = r + 2; i <= n; ++i)
        {
            d += (a->elem[i][r]) * (a->elem[i][r]);
        }
        if (d == 0)
        {
            //不需要拟对角化
            continue;
        }
        c = sqrt(d + (a->elem[r + 1][r]) * (a->elem[r + 1][r]));
        if (a->elem[r + 1][r] > 0)
        {
            c = (-1) * c;
        }
        h = c * (c - a->elem[r + 1][r]);
        //求U
        for (i = 1; i <= n; ++i)
        {
            if (i <= r)
            {
                U.elem[i] = 0;
            }
            else if (i == r + 1)
            {
                U.elem[i] = a->elem[i][r] - c;
            }
            else
            {
                U.elem[i] = a->elem[i][r];
            }
        }
        //计算A(r+1)
        Vector_Num(U, 1 / h, &V); //V=U/h
        V_M(V, *a, &P);           //P=UA
        M_V(*a, V, &Q);           //Q=AU
        double tmp = DotProduct(P, V);
        Vector_Num(U, tmp, &Vtemp);
        VmV(Q, Vtemp, &Vtemp); //Vtemp=Q-V(PU)
        V_V(Vtemp, U, &Mtemp);
        MmM(*a, Mtemp, a); //A=(A-VP)-(Vtemp*V)
        V_V(U, P, &Mtemp);
        MmM(*a, Mtemp, a); //A=(A-VP)-(Vtemp*V)
    }
    free(U.elem);
    free(V.elem);
    free(P.elem);
    free(Q.elem);
    free(Vtemp.elem);
    for (i = 1; i <= n; ++i)
    {
        free(Mtemp.elem[i]);
    }
    free(Mtemp.elem);
    return;
}
//QR分解
void QR_dcp(Matrix *a, Matrix *q)
{
    int i, j, r, n;
    if (a->m != a->n || q->m != q->n)
    {
        printf("a 或 q不是方阵");
        exit(MISMATCH);
    }
    if (a->m == q->m)
    {
        n = a->m;
    }
    else
    {
        exit(MISMATCH);
    }
    double d, c, h;
    Vector U;
    VectorInit(&U, n);
    Vector V;
    VectorInit(&V, n);
    Vector P;
    VectorInit(&P, n);
    Matrix M;
    Matrix_Init(&M, n, n);

    //Q=I
    for (i = 1; i <= n; ++i)
    {
        for (j = 1; j <= n; ++j)
        {
            q->elem[i][j] = (i == j) ? 1 : 0;
        }
    }
    // PrintMatrix(*q);
    for (r = 1; r <= n - 1; ++r)
    {
        //求c,h
        for (d = 0, i = r + 1; i <= n; ++i)
        {
            d += (a->elem[i][r]) * (a->elem[i][r]);
        }
        if (d == 0)
        {
            continue; //无需QR
        }
        c = sqrt(d + (a->elem[r][r]) * (a->elem[r][r]));
        if (a->elem[r][r] > 0)
        {
            c = (-1) * c; //令c的符合与A[r][r]相反
        }
        h = c * (c - a->elem[r][r]);
        //求U
        for (i = 1; i <= n; ++i)
        {
            if (i < r)
            {
                U.elem[i] = 0;
            }
            else if (i == r)
            {
                U.elem[i] = a->elem[i][r] - c;
            }
            else
            {
                U.elem[i] = a->elem[i][r];
            }
        }
        //计算Q，R：
        Vector_Num(U, 1 / h, &V);
        M_V(*q, U, &P);
        V_V(P, V, &M);
        MmM(*q, M, q);
        V_M(U, *a, &P);
        V_V(V, P, &M);
        MmM(*a, M, a);
    }
    free(U.elem);
    free(V.elem);
    free(P.elem);
    for (i = 1; i < n; ++i)
    {
        free(M.elem[i]);
    }
    free(M.elem);
}
//QR方法，A最终变为对角线一二阶子块的分块上三角阵
void QR(Matrix *a, Matrix *q)
{
    int i, j, r, n, count, flag;
    if (a->m != a->n || q->m != q->n)
    {
        printf("a 或 q不是方阵");
        exit(MISMATCH);
    }
    if (a->m == q->m)
    {
        n = a->m;
    }
    else
    {
        exit(MISMATCH);
    }
    Matrix B;
    Matrix_Init(&B, n, n);
    flag = n, count = 0;
    while (flag > 1)
    {
        QR_dcp(a, q);
        M_M(*a, *q, &B);
        if (flag >= 2 && fabs(B.elem[flag][flag - 1]) < E)
        {
            flag--;
        }
        else if (flag >= 3 && fabs(B.elem[flag - 1][flag - 2]) < E && (B.elem[flag - 1][flag - 1] - B.elem[flag][flag]) * (B.elem[flag - 1][flag - 1] - B.elem[flag][flag]) + 4 * B.elem[flag][flag - 1] * B.elem[flag - 1][flag] < 0)
        {
            flag -= 2;
        }
        ++count;
        MatrixCopy(B, a);
    }
    for (i = 1; i <= n; ++i)
    {
        free(B.elem[i]);
    }
    free(B.elem);
    return;
}
//双步位移QR，对M进行分解
void QR2Tran_QR(Matrix *b, Matrix *a)
{
    int n, i, r;
    double d, c, temp, h;
    if (a->m != a->n || b->m != b->n)
    {
        printf("a 或 b不是方阵");
        exit(MISMATCH);
    }
    if (a->m == b->m)
    {
        n = a->m;
    }
    else
    {
        exit(MISMATCH);
    }
    Vector U;
    VectorInit(&U, n);

    Vector V;
    VectorInit(&V, n);
    Vector P;
    VectorInit(&P, n);
    Vector Q;
    VectorInit(&Q, n);

    Vector Vtemp;
    VectorInit(&Vtemp, n);
    Matrix Mtemp;
    Matrix_Init(&Mtemp, n, n);
    for (r = 1; r <= n - 1; ++r)
    {
        //求c,h
        for (d = 0, i = r + 1; i <= n; ++i)
        {
            d += (b->elem[i][r]) * (b->elem[i][r]);
        }
        if (d == 0)
        {
            continue;
        }
        //求c
        c = sqrt(d + (b->elem[r][r]) * (b->elem[r][r]));
        if (b->elem[r][r] > 0)
        {
            c = (-1) * c;
        }
        h = c * (c - b->elem[r][r]);
        //求U
        for (i = 1; i <= n; ++i)
        {
            if (i <= r - 1)
            {
                U.elem[i] = 0;
            }
            else if (i == r)
            {
                U.elem[i] = b->elem[i][r] - c;
            }
            else
            {
                U.elem[i] = b->elem[i][r];
            }
        }
        //计算B(r+1)
        // Vector_Num(U, 1 / h, &V);
        // M_V(*b, V, &P);
        // V_V(U, P, &Mtemp);
        // MmM(*b, Mtemp, b);
        // V_M(V, *a, &P);
        // M_V(*a, V, &Q);
        // double tmp = DotProduct(P, V);
        // Vector_Num(U, tmp, &Vtemp);
        // VmV(Q, Vtemp, &Vtemp);
        // V_V(Vtemp, U, &Mtemp);
        // MmM(*a, Mtemp, a);
        // V_V(U, P, &Mtemp);
        // MmM(*a, Mtemp, a);
        Vector_Num(U, 1 / h, &V);
        M_V(*b, V, &P);
        V_V(U, P, &Mtemp);
        MmM(*b, Mtemp, b);
        V_M(U, *a, &P);
        M_V(*a, U, &Q);
        V_V(V, P, &Mtemp);
        MmM(*a, Mtemp, a);
        double tmp = DotProduct(P, V);
        Vector_Num(V, tmp, &Vtemp);
        VmV(Q, Vtemp, &Vtemp);
        V_V(Vtemp, V, &Mtemp);
        MmM(*a, Mtemp, a);
    }
    free(U.elem);
    free(V.elem);
    free(P.elem);
    free(Q.elem);
    free(Vtemp.elem);
    for (i = 1; i <= n; ++i)
    {
        free(Mtemp.elem[i]);
    }
    free(Mtemp.elem);
    return;
}