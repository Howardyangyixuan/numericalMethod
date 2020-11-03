#include "myVector.h"
void VectorInit(Vector *v, int n)
{
    v->dimension = n;
    v->elem = (ElemType *)malloc(sizeof(ElemType) * (n + 1));
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