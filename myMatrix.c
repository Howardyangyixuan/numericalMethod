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
            printf("%.12e ", A.elem[i][j]);
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

void Trans(Matrix *a, double p)
{
    if (a->m != a->n)
    {
        printf("a不是方阵");
        exit(MISMATCH);
    }
    for (int i = 1; i <= a->n; ++i)
    {
        a->elem[i][i] -= p;
    }
    return;
}

//高斯消去Ax = b, 求解x
void GaussElimination(Matrix A, Vector b, Vector *x)
{
    if (A.m != A.n)
    {
        printf("a不是方阵");
        exit(MISMATCH);
    }
    int n = A.m;
    int i, j, k, maxline;
    double maxUnit, m;
    //消元过程
    for (k = 1; k <= n; ++k)
    {
        maxUnit = fabs(A.elem[k][k]);
        maxline = k;
        for (j = k + 1; j <= n; ++j)
        {
            if (fabs(A.elem[j][k]) > maxUnit)
            {
                maxUnit = fabs(A.elem[j][k]);
                maxline = j;
            }
        }
        //交换两行
        double tmp;
        for (j = 1; j <= n; ++j)
        {
            tmp = A.elem[k][j];
            A.elem[k][j] = A.elem[maxline][j];
            A.elem[maxline][j] = tmp;
        }
        tmp = b.elem[k];
        b.elem[k] = b.elem[maxline];
        b.elem[maxline] = tmp;
        //如果A[k][k]太小，做分母会引入计算误差，因列主元过小时,直接跳过进入回带过程
        if (fabs(A.elem[k][k]) < E)
        {
            continue;
        }
        //开始消元
        for (i = k + 1; i <= n; ++i)
        {
            m = A.elem[i][k] / A.elem[k][k];
            for (j = k + 1; j <= n; ++j)
            {
                A.elem[i][j] -= m * A.elem[k][j];
            }
            b.elem[i] -= m * b.elem[k];
        }
    }
    //回带过程
    x->elem[n] = b.elem[n] / A.elem[n][n];
    for (k = n - 1; k >= 1; --k)
    {
        x->elem[k] = b.elem[k];
        for (j = k + 1; j <= n; ++j)
        {

            x->elem[k] -= A.elem[k][j] * x->elem[j];
        }
        x->elem[k] /= A.elem[k][k];
    }
}

//幂法(求矩阵按模最大的特征值):
double PM_eigenvalue(Matrix A)
{
    if (A.m != A.n)
    {
        printf("a不是方阵");
        exit(MISMATCH);
    }
    int const n = A.m;
    int i, count;
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
        double tmp = DotProduct(U, U);
        norm = sqrt(tmp);
        Vector_Num(U, 1 / norm, &V);
        M_V(A, V, &U);
        lambda_old = lambda_new;
        lambda_new = DotProduct(V, U);
        e = fabs((lambda_new - lambda_old) / lambda_new);
        ++count;
    }
    //首项归一化
    // Vector_Num(V, 1 / V.elem[1], &V);
    PrintVector(V);
    free(V.elem);
    free(U.elem);
    return lambda_new;
}

//反幂法
double IPM_eigenvalue(Matrix A)
{
    if (A.m != A.n)
    {
        printf("a不是方阵");
        exit(MISMATCH);
    }
    int const n = A.m;
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
    Vector tmp;
    VectorInit(&tmp, n);
    Matrix m;
    Matrix_Init(&m, n, n);
    do
    {
        norm = sqrt(DotProduct(U, U));
        Vector_Num(U, 1 / norm, &V);
        VectorCopy(V, &tmp);
        MatrixCopy(A, &m);
        GaussElimination(m, tmp, &U);
        lambda_old = lambda_new;
        lambda_new = DotProduct(V, U);
        e = fabs((lambda_new - lambda_old) / lambda_old); //求相对误差
        ++count;
    } while (e > E);
    PrintVector(V);
    free(V.elem);
    free(U.elem);
    return 1 / lambda_new;
}
//幂法平移(求矩阵距离p最远的特征值，并打印特征向量)
double PMTran_eigenvalue(Matrix A, double p)
{
    double l;
    Trans(&A, p);
    l = PM_eigenvalue(A);
    Trans(&A, (-1) * p);
    return l + p;
}
//反幂法平移(求矩阵距离p最近的特征值，并打印特征向量)
double IPMTran_eigenvalue(Matrix A, double p)
{
    double l;
    Trans(&A, p);
    l = IPM_eigenvalue(A);
    Trans(&A, (-1) * p);
    return l + p;
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
            printf("%.12e ", A.elem[i][j]);
            // printf("%.12e& ", A.elem[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    return;
}
//方阵 * 向量 A * x = b
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
        double tmp = DotProduct(P, U);
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
//双步位移QR方法中，解二次方程，并存储两个特征值
void QR2Tran_solve(double x11, double x12, double x21, double x22, complex *lambda)
{
    double det, tr, delta;
    tr = x11 + x22;
    det = x11 * x22 - x12 * x21;
    delta = tr * tr - 4 * det;
    if (delta >= 0)
    {
        printf("%.12e\n", (sqrt(delta) + tr) / 2);
        lambda[(int)lambda[0].real].real = (sqrt(delta) + tr) / 2;
        lambda[(int)lambda[0].real++].imaginary = 0;
        printf("%.12e\n", (-sqrt(delta) + tr) / 2);
        lambda[(int)lambda[0].real].real = (-sqrt(delta) + tr) / 2;
        lambda[(int)lambda[0].real++].imaginary = 0;
    }
    else
    {
        printf("%.12e + i * %.12e\n", tr / 2, sqrt(-delta) / 2);
        lambda[(int)lambda[0].real].real = tr / 2;
        lambda[(int)lambda[0].real++].imaginary = sqrt(-delta) / 2;
        printf("%.12e + i * %.12e\n", tr / 2, -sqrt(-delta) / 2);
        lambda[(int)lambda[0].real].real = tr / 2;
        lambda[(int)lambda[0].real++].imaginary = -sqrt(-delta) / 2;
    }
    return;
}
//双步位移QR方法中，求全部特征值，自带拟上三角化
void QR2Tran(Matrix A, complex *lambda)
{
    if (A.m != A.n)
    {
        printf("a不是方阵");
        exit(MISMATCH);
    }
    //最大迭代次数
    Hess(&A);
    // PrintMatrix(A);
    int count = 1;
    int n = A.m;
    int loc = lambda[0].real;

    double tr, det;
    Matrix A1;
    Matrix_Init(&A1, n, n);
    Matrix M;
    Matrix_Init(&M, n, n);

    while (1)
    {
        //降一阶
        if (n == 1 || fabs(A.elem[n][n - 1]) < E)
        {
            printf("%.12e\n", A.elem[n][n]);
            lambda[(int)lambda[0].real].real = A.elem[n][n];
            lambda[(int)lambda[0].real++].imaginary = 0; //储存特征值
            free(A.elem[n]);
            free(A1.elem[n]);
            free(M.elem[n]);
            --n;
            M.m = A1.m = A.m = n;
            M.n = A1.n = A.n = n;
            if (n == 0)
            {
                free(A.elem);
                free(A1.elem);
                free(M.elem);
                break;
            }
            continue;
        }
        //降二阶
        if (n == 2 || fabs(A.elem[n - 1][n - 2]) < E)
        {
            //存储特征值
            QR2Tran_solve(A.elem[n - 1][n - 1], A.elem[n - 1][n], A.elem[n][n - 1], A.elem[n][n], lambda);
            free(A.elem[n]);
            free(A.elem[n - 1]);
            free(A1.elem[n]);
            free(A1.elem[n - 1]);
            free(M.elem[n]);
            free(M.elem[n - 1]);
            n -= 2;
            //矩阵降维
            M.m = A1.m = A.m = n;
            M.n = A1.n = A.n = n;
            if (n == 0)
            {
                free(A.elem);
                free(A1.elem);
                free(M.elem);
                break;
            }
            continue;
        }
        //QR分解
        tr = A.elem[n - 1][n - 1] + A.elem[n][n];
        det = A.elem[n - 1][n - 1] * A.elem[n][n] - A.elem[n - 1][n] * A.elem[n][n - 1];
        MatrixCopy(A, &A1);
        Trans(&A1, tr);
        M_M(A, A1, &M);
        Trans(&M, -det);
        QR2Tran_QR(&M, &A);
        ++count;
        if (count >= L)
        {
            printf("迭代次数过多，结束\n");
            exit(TOOMANYTIMES);
        }
    }
}
//Newton法，解非线性方程
void NewtonMethod(double x, double y, double TU[2])
{
    double t, u, v, w;
    Vector F;
    Matrix f;
    Vector X;
    Vector Y;
    VectorInit(&F, 4);
    Matrix_Init(&f, 4, 4);
    VectorInit(&X, 4);
    VectorInit(&Y, 4);
    //初始向量X，可变
    X.elem[1] = 0.5;
    X.elem[2] = 1;
    X.elem[3] = y - 1;
    X.elem[4] = x + 2;
    //初始值
    t = X.elem[1];
    u = X.elem[2];
    v = X.elem[3];
    w = X.elem[4];
    //-F
    F.elem[1] = (0.5 * cos(t) + u + v + w - x - 2.67) * (-1);
    F.elem[2] = (t + 0.5 * sin(u) + v + w - y - 1.07) * (-1);
    F.elem[3] = (0.5 * t + u + cos(v) + w - x - 3.74) * (-1);
    F.elem[4] = (t + 0.5 * u + v + sin(w) - y - 0.79) * (-1);
    //f
    f.elem[1][1] = -0.5 * sin(t);
    f.elem[1][2] = 1;
    f.elem[1][3] = 1;
    f.elem[1][4] = 1;
    f.elem[2][1] = 1;
    f.elem[2][2] = 0.5 * cos(u);
    f.elem[2][3] = 1;
    f.elem[2][4] = 1;
    f.elem[3][1] = 0.5;
    f.elem[3][2] = 1;
    f.elem[3][3] = -sin(v);
    f.elem[3][4] = 1;
    f.elem[4][1] = 1;
    f.elem[4][2] = 0.5;
    f.elem[4][3] = 1;
    f.elem[4][4] = cos(w);
    //
    double nx, ny;
    double e = 1;
    int count = 0;
    while (e > E)
    {
        t = X.elem[1];
        u = X.elem[2];
        v = X.elem[3];
        w = X.elem[4];
        F.elem[1] = (0.5 * cos(t) + u + v + w - x - 2.67) * (-1);
        F.elem[2] = (t + 0.5 * sin(u) + v + w - y - 1.07) * (-1);
        F.elem[3] = (0.5 * t + u + cos(v) + w - x - 3.74) * (-1);
        F.elem[4] = (t + 0.5 * u + v + sin(w) - y - 0.79) * (-1);
        f.elem[1][1] = -0.5 * sin(t);
        f.elem[2][2] = 0.5 * cos(u);
        f.elem[3][3] = -sin(v);
        f.elem[4][4] = cos(w);
        GaussElimination(f, F, &Y);
        ny = norm_inf(Y);
        printf("ny:%lf\n", ny);
        //x+y
        Vector_Num(Y, -1, &Y);
        VmV(X, Y, &X);
        nx = norm_inf(X);
        e = ny / nx;
        if (++count > L)
        {
            printf("迭代次数过多，结束\n");
            exit(TOOMANYTIMES);
        }
    }
    //(t,u)
    TU[0] = X.elem[1];
    printf("TU[0]:%lf", TU[0]);
    TU[1] = X.elem[2];
    printf("TU[1]:%lf", TU[1]);
    int i;
    free(X.elem);
    free(Y.elem);
    free(F.elem);
    for (i = 1; i <= 4; ++i)
    {
        free(f.elem[i]);
    }
    free(f.elem);
    return;
}
//矩阵转置,B = AT
void Matrix_T(Matrix A, Matrix *b)
{
    int m, n, i, j, k;
    if (b->m == A.n && b->n == A.m)
    {
        m = A.m;
        n = A.n;
    }
    else
    {
        exit(MISMATCH);
    }
    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            b->elem[j][i] = A.elem[i][j];
        }
    }
}
//向量Y=方阵A*向量X (m,n)*(n,1)=(m,1)
void MtV(Matrix A, Vector X, Vector *y)
{
    int m, n, i, j;
    if (X.dimension == A.n && y->dimension == A.m)
    {
        m = A.m;
        n = A.n;
    }
    else
    {
        exit(MISMATCH);
    }
    for (i = 0; i < m; ++i)
    {
        for (y->elem[i + 1] = 0, j = 0; j < n; ++j)
        {
            y->elem[i + 1] += A.elem[i][j] * X.elem[j + 1];
        }
    }
    return;
}
//解矩阵方程AXB (m,m)*(m,n)=(m,n)
void Gauss_matrix(Matrix A, Matrix B, Matrix *x)
{
    int m, n, i, j;
    if (x->n == B.n && x->m == A.m)
    {
        m = x->m;
        n = x->n;
    }
    else
    {
        exit(MISMATCH);
    }
    Vector VX;
    Vector VB;
    VectorInit(&VX, m);
    VectorInit(&VB, m);
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            VB.elem[i + 1] = B.elem[i][j];
        }
        GaussElimination(A, VB, &VX);
        for (i = 0; i < m; ++i)
        {
            x->elem[i][j] = VX.elem[i + 1];
        }
    }
    free(VB.elem);
    free(VX.elem);
    return;
}
//插值(t,u)->z
double Interpotation(double t, double u, double T[6], double U[6], double Z[6][6])
{
    int jt, ju, i, j;
    double z;
    double l__jt, l_jt_, ljt__;
    double s__ju, s_ju_, sju__;
    //求出t,u落在的插值中节点(jt,ju)
    switch ((int)(t / 0.2 + 0.5))
    {
    case 0:
    case 1:
    {
        jt = 1;
        break;
    }
    case 2:
    {
        jt = 2;
        break;
    }
    case 3:
    {
        jt = 3;
        break;
    }
    case 4:
    case 5:
    {
        jt = 4;
        break;
    }
    default:
        exit(ERROR);
        break;
    }
    switch ((int)(u / 0.4 + 0.5))
    {
    case 0:
    case 1:
    {
        ju = 1;
        break;
    }
    case 2:
    {
        ju = 2;
        break;
    }
    case 3:
    {
        ju = 3;
        break;
    }
    case 4:
    case 5:
    {
        ju = 4;
        break;
    }
    default:
        exit(ERROR);
        break;
    }
    printf("Lagrange\n");
    //求Lagrange插值基函数
    l__jt = (t - T[jt]) * (t - T[jt + 1]) / (T[jt - 1] - T[jt]) / (T[jt - 1] - T[jt + 1]);
    l_jt_ = (t - T[jt - 1]) * (t - T[jt + 1]) / (T[jt] - T[jt - 1]) / (T[jt] - T[jt + 1]);
    ljt__ = (t - T[jt - 1]) * (t - T[jt]) / (T[jt + 1] - T[jt - 1]) / (T[jt + 1] - T[jt]);
    s__ju = (u - U[ju]) * (u - U[ju + 1]) / (U[ju - 1] - U[ju]) / (U[ju - 1] - U[ju + 1]);
    s_ju_ = (u - U[ju - 1]) * (u - U[ju + 1]) / (U[ju] - U[ju - 1]) / (U[ju] - U[ju + 1]);
    sju__ = (u - U[ju - 1]) * (u - U[ju]) / (U[ju + 1] - U[ju - 1]) / (U[ju + 1] - U[ju]);
    //求z
    z = l__jt * (Z[jt - 1][ju - 1] * (s__ju) + Z[jt - 1][ju] * (s_ju_) + Z[jt - 1][ju + 1] * (sju__)) + l_jt_ * (Z[jt][ju - 1] * (s__ju) + Z[jt][ju] * (s_ju_) + Z[jt][ju + 1] * (sju__)) +
        ljt__ * (Z[jt + 1][ju - 1] * (s__ju) + Z[jt + 1][ju] * (s_ju_) + Z[jt + 1][ju + 1] * (sju__));
    return z;
}
//z=f(x,y)
double func(double x, double y, double T[6], double U[6], double Z[6][6])
{
    double TU[2];
    double z;
    //(x,y)->(t,u)
    NewtonMethod(x, y, TU);
    printf("Inter");
    //(t,u)->z
    z = Interpotation(TU[0], TU[1], T, U, Z);
    return z;
}
//求解满足精度要求的拟合次数，计算并输出正交基下的拟合系数矩阵，返回p*ij
int fitting(double f[11][21], double p_star[9][6])
{
    int i, j, t1, t2;
    double x[11];
    double y[21];
    double P[11][21];
    Matrix F;
    Matrix_Init(&F, 11, 21);
    for (i = 0; i <= 10; ++i)
    {
        x[i] = 0.08 * i;
    }
    for (j = 0; j <= 20; ++j)
    {
        y[j] = 0.5 + 0.05 * j;
    }
    for (i = 0; i <= 10; ++i)
    {
        for (j = 0; j <= 20; ++j)
        {
            F.elem[i][j] = f[i][j];
        }
    }
    double sigma;
    int k = 2; //k-1为拟合次数
    while (1)
    {
        //求正交基函数
        //求矩阵BG：
        Matrix BT;
        Matrix_Init(&BT, k, 11);
        Matrix GT;
        Matrix_Init(&GT, k, 21);
        for (i = 0; i <= k - 1; ++i)
        {
            for (j = 0; j <= 10; ++j)
            {
                BT.elem[i][j] = pow(x[j], i);
            }
        }
        for (i = 0; i <= k - 1; ++i)
        {
            for (j = 0; j <= 20; ++j)
            {
                GT.elem[i][j] = pow(y[j], i);
            }
        }
        //解方程并求出系数矩阵C
        Matrix B;
        Matrix_Init(&B, 11, k);
        Matrix G;
        Matrix_Init(&G, 21, k);
        Matrix U;
        Matrix_Init(&U, k, 21);
        Matrix BTB;
        Matrix_Init(&BTB, k, k);
        Matrix GTG;
        Matrix_Init(&GTG, k, k);
        Matrix A;
        Matrix_Init(&A, k, 21);
        Matrix DT;
        Matrix_Init(&DT, k, 21);
        Matrix D;
        Matrix_Init(&D, 21, k);
        Matrix C;
        Matrix_Init(&C, k, k);
        //计算
        Matrix_T(BT, &B);
        Matrix_T(GT, &G);
        M_M(BT, B, &BTB);
        M_M(GT, G, &GTG);
        M_M(BT, F, &U);
        Gauss_matrix(BTB, U, &A);
        Gauss_matrix(GTG, GT, &GT);
        Matrix_T(DT, &D);
        M_M(A, D, &C);
        //free
        for (i = 0; i < BT.m; ++i)
            free(BT.elem[i]);
        free(BT.elem);
        for (i = 0; i < GT.m; ++i)
            free(GT.elem[i]);
        free(GT.elem);
        for (i = 0; i < B.m; ++i)
            free(B.elem[i]);
        free(B.elem);
        for (i = 0; i < G.m; ++i)
            free(G.elem[i]);
        free(G.elem);
        for (i = 0; i < DT.m; ++i)
            free(DT.elem[i]);
        free(BT.elem);
        for (i = 0; i < D.m; ++i)
            free(D.elem[i]);
        free(D.elem);
        for (i = 0; i < U.m; ++i)
            free(U.elem[i]);
        free(U.elem);
        for (i = 0; i < A.m; ++i)
            free(A.elem[i]);
        free(A.elem);
        for (i = 0; i < BTB.m; ++i)
            free(BTB.elem[i]);
        free(BTB.elem);
        for (i = 0; i < GTG.m; ++i)
            free(GTG.elem[i]);
        free(GTG.elem);
        //计算p[i][j]
        for (i = 0; i <= 10; ++i)
        {
            for (j = 0; j <= 20; ++j)
            {
                for (P[i][j] = 0, t1 = 0; t1 <= k - 1; ++t1)
                {
                    for (t2 = 0; t2 <= k - 1; ++t2)
                    {
                        P[i][j] += C.elem[t1 + 1][t2 + 1] * pow(x[i], t1) * pow(y[j], t2);
                    }
                }
            }
        }
        //判断误差
        sigma = 0;
        for (i = 0; i <= 10; ++i)
        {
            for (j = 0; j <= 20; ++j)
            {
                sigma += (F.elem[i][j] - P[i][j]) * (F.elem[i][j] - P[i][j]);
            }
        }
        if (sigma < 1.0e-7)
        {
            printf("\n最终的k和o:\n");
            printf("k = %d\t\to = %.12e\n", k - 1, sigma);
            printf("\n系数矩阵((cij):\n");
            PrintMatrix(C);
            //计算p*[i][j]
            for (i = 0; i <= 8; ++i)
            {
                for (j = 1; j <= 5; ++j)
                {
                    for ((p_star)[i][j] = 0, t1 = 0; t1 <= k - 1; ++t1)
                    {
                        for (t2 = 0; t2 <= k - 1; ++t2)
                        {
                            (p_star)[i][j] += C.elem[t1 + 1][t2 + 1] * pow(0.1 * i, t1) * pow(0.5 + 0.2 * j, t2);
                        }
                    }
                }
            }
            for (i = 1; i <= C.m; ++i)
            {
                free(C.elem[i]);
            }
            free(C.elem);
            break;
        }
        else
        {
            //输出选择过程中的k和o
            printf("k = %d\t\to = %.12e\n", k - 1, sigma);
            ++k;
            for (i = 1; i <= C.m; ++i)
            {
                free(C.elem[i]);
            }
            free(C.elem);
        }
    }
    //把C换算成系数crs
    for (i = 0; i < F.m; ++i)
    {
        free(F.elem[i]);
    }
    free(F.elem);
    return k - 1;
}
//方阵C=矩阵A*矩阵B (m,n)*(n,m)=(m,m)
void mmM(Matrix A, Matrix B, Matrix *c)
{
    int m, n, i, j, k;
    if (B.m == A.n && A.m == B.n && c->m == A.m && c->n == A.n)
    {
        m = A.m;
        n = A.n;
    }
    else
    {
        exit(MISMATCH);
    }
    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            for (c->elem[i + 1][j + 1] = 0, k = 0; k < n; ++k)
            {
                c->elem[i + 1][j + 1] += A.elem[i][k] * B.elem[k][j];
            }
        }
    }
    return;
}
//矩阵C=矩阵A*矩阵B (m,k)*(k,n)=(m,n)
void mmm(Matrix A, Matrix B, Matrix *c)
{
    int m, n, k, i, j, t;
    if (B.m == A.n && c->m == A.m && c->n == B.n)
    {
        m = A.m;
        n = B.n;
        k = A.n;
    }
    else
    {
        exit(MISMATCH);
    }
    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            for (c->elem[i][j] = 0, t = 0; t < k; ++t)
            {
                c->elem[i][j] += A.elem[i][t] * B.elem[t][j];
            }
        }
    }
    return;
}