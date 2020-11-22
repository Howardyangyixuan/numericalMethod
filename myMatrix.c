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
        printf("%lf\n", (sqrt(delta) + tr) / 2);
        lambda[(int)lambda[0].real].real = (sqrt(delta) + tr) / 2;
        lambda[(int)lambda[0].real++].imaginary = 0;
        printf("%lf\n", (-sqrt(delta) + tr) / 2);
        lambda[(int)lambda[0].real].real = (-sqrt(delta) + tr) / 2;
        lambda[(int)lambda[0].real++].imaginary = 0;
    }
    else
    {
        printf("%lf + i * %lf\n", tr / 2, sqrt(-delta) / 2);
        lambda[(int)lambda[0].real].real = tr / 2;
        lambda[(int)lambda[0].real++].imaginary = sqrt(-delta) / 2;
        printf("%lf + i * %lf\n", tr / 2, -sqrt(-delta) / 2);
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
    const int L = 10e3;
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
            printf("%lf\n", A.elem[n][n]);
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
            return;
        }
    }
}