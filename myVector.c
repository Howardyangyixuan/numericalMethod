#include "myVector.h"
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
//输入&V，n；初始化n维向量V（用1~n):
void VectorInit(Vector *v, int n)
{
    v->dimension = n;
    v->elem = (ElemType *)malloc(sizeof(ElemType) * (n + 1));
    return;
}
//复制向量
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

//打印向量
void PrintVector(Vector V)
{
    for (int i = 1; i <= V.dimension; ++i)
    {
        printf("%lf ", V.elem[i]);
    }
    printf("\n");
    return;
}

//向量数乘，V=向量U*数k（要求V已经初始化），传入V的指针，改变V后，返回V的指针
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

//向量内积
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

//减法X - Y = V
void VmV(Vector X, Vector Y, Vector *v)
{
    if (X.dimension != Y.dimension)
    {
        printf("维数不匹配");
        exit(MISMATCH);
    }
    for (int i = 1; i <= X.dimension; i++)
    {
        v->elem[i] = X.elem[i] - Y.elem[i];
    }
}
//向量的无穷范数
double norm_inf(Vector X)
{
    int n = X.dimension;
    int i;
    double M = fabs(X.elem[1]);
    for (i = 2.; i <= n; ++i)
    {
        M = max(M, fabs(X.elem[i]));
    }
    return M;
}