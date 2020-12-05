#include "myMatrix.h"
int main()
{
    int i, j;
    //定义T[i] U[j] Z[ti][uj]:
    double T[6];
    double U[6];
    double Z[6][6];
    for (i = 0; i < 6; ++i)
    {
        T[i] = 0.2 * i;
    }
    for (j = 0; j < 6; ++j)
    {
        U[j] = 0.4 * j;
    }
    Z[0][0] = -0.5;
    Z[0][1] = -0.34;
    Z[0][2] = 0.14;
    Z[0][3] = 0.94;
    Z[0][4] = 2.06;
    Z[0][5] = 3.5;
    Z[1][0] = -0.42;
    Z[1][1] = -0.5;
    Z[1][2] = -0.26;
    Z[1][3] = 0.3;
    Z[1][4] = 1.18;
    Z[1][5] = 2.38;
    Z[2][0] = -0.18;
    Z[2][1] = -0.5;
    Z[2][2] = -0.5;
    Z[2][3] = -0.18;
    Z[2][4] = 0.46;
    Z[2][5] = 1.42;
    Z[3][0] = 0.22;
    Z[3][1] = -0.34;
    Z[3][2] = -0.58;
    Z[3][3] = -0.5;
    Z[3][4] = -0.1;
    Z[3][5] = 0.62;
    Z[4][0] = 0.78;
    Z[4][1] = -0.02;
    Z[4][2] = -0.5;
    Z[4][3] = -0.66;
    Z[4][4] = -0.5;
    Z[4][5] = -0.02;
    Z[5][0] = 1.5;
    Z[5][1] = 0.46;
    Z[5][2] = -0.26;
    Z[5][3] = -0.66;
    Z[5][4] = -0.74;
    Z[5][5] = -0.5;
    //计算f(xi,yj)并输出数表(xi，xj，f(xi,yj))：
    double f[11][21];
    for (i = 0; i <= 10; ++i)
    {
        for (j = 0; j <= 20; ++j)
        {
            f[i][j] = func(0.08 * i, 0.5 + 0.05 * j, T, U, Z);
        }
    }
    printf("\n数表：(xi，xj，f(xi,yj))\n");
    for (i = 0; i <= 10; ++i)
    {
        for (j = 0; j <= 20; ++j)
        {
            printf("%lf, %lf, %.12e\t", 0.08 * i, 0.05 * j, f[i][j]);
        }
        printf("\n");
    }
    //计算f*(xi,yj)
    double f_star[9][6];
    for (i = 1; i <= 8; ++i)
    {
        for (j = 1; j <= 5; ++j)
        {
            f_star[i][j] = func(0.1 * i, 0.5 + 0.2 * j, T, U, Z);
        }
    }
    //计算并输出ij,计算p(x,y),输出数表
    double p_star[9][6];
    printf("\n选择过程中k和o：\n");
    fitting(f, p_star);
    printf("\n数表:(x*j,y*j,f*(xi,yj),p*(i,j))\n");
    for (i = 1; i <= 8; ++i)
    {
        for (j = 1; j <= 5; ++j)
        {
            printf("%lf, %lf,%.12e,%.12e\t", 0.1 * i, 0.5 + 0.2 * j, f_star[i][j], p_star[i][j]);
        }
    }
    return 0;
}