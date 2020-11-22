#include "myMatrix.h"
int main()
{
    Matrix a;
    int n = 501;
    double c = -0.064;
    double b = 0.16;
    Matrix_Init(&a, n, n);
    for (int i = 1; i <= n; ++i)
    {
        a.elem[i][i] = (1.64 - 0.024 * i) * sin(0.2 * i) - 0.64 * exp(0.1 / i);
        if (i + 1 <= n && i >= 1)
            a.elem[i + 1][i] = b;
        if (i - 1 >= 1 && i <= n)
            a.elem[i - 1][i] = b;
        if (i + 2 <= n && i >= 1)
            a.elem[i + 2][i] = c;
        if (i - 2 >= 1 && i <= n)
            a.elem[i - 2][i] = c;
    }
    double tmp = PM_eigenvalue(a);
    printf("PM result: %.12e\n", tmp);
    printf("PMTrans result: %.12e\n", PMTran_eigenvalue(a, tmp));
    printf("IPMTrans result: %.12e\n", IPMTran_eigenvalue(a, tmp));
    return 0;
}