#include "myVector.h"
int main()
{
    Vector a;
    VectorInit(&a, 1);
    a.elem[1] = 1;
    PrintVector(a);
    return 0;
}