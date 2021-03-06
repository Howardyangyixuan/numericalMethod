#ifndef _MY_MATRIX_H_
#define _MY_MATRIX_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myVector.h"

typedef struct
{
    double real;
    double imaginary;
} complex;
typedef struct
{
    int dimension; //矩阵维数
    int s;         //上半带宽
    int r;         //下半带宽
    double **elem; //矩阵，二维数组
} Matrix_diag;     //n维s+1+r对角方阵

typedef struct
{
    int m;
    int n;
    double **elem;
} Matrix;

//输入&M，n：初始化n维s+1+r对角方阵M(用(1~s+1+r)*(1~n)):
void Matrix_diagInit(Matrix_diag *a, int n, int s, int r);
//输出带状方阵A:
void PrintMatrix_diag(Matrix_diag A);
//把A值复制给A1（A为带状阵，A1不需要初始化）
void Matrix_diagCopy(Matrix_diag A, Matrix_diag *a1);
//矩阵平移(A->A-pI)
void Trans(Matrix *a, double p);
//高斯消去Ax = b, 求解x
void GaussElimination(Matrix A, Vector b, Vector *x);
//幂法(求矩阵按模最大的特征值)
double PM_eigenvalue(Matrix A);
//幂法平移(求矩阵距离p最远的特征值，并打印特征向量)
double PMTran_eigenvalue(Matrix A, double p);
//反幂法(求矩阵按模最小的特征值，并打印特征向量)
double IPM_eigenvalue(Matrix A);
//反幂法平移(求矩阵距离p最近的特征值，并打印特征向量)
double IPMTran_eigenvalue(Matrix A, double p);
//输入&M，n：初始化m*n维矩阵
void Matrix_Init(Matrix *a, int m, int n);
//输出矩阵A:
void PrintMatrix(Matrix A);
//矩阵 * 向量 A * x = b
void M_V(Matrix A, Vector x, Vector *b);
//向量 * 矩阵 x * A = b
void V_M(Vector v, Matrix A, Vector *b);
//矩阵 * 矩阵 A * B = C
void M_M(Matrix A, Matrix B, Matrix *c);
//向量 * 向量 a * b = C
Matrix *V_V(Vector a, Vector b, Matrix *c);
//矩阵 - 矩阵 A - B = C
void MmM(Matrix A, Matrix B, Matrix *c);
//拟上三角化
void Hess(Matrix *a);
//QR分解
void QR_dcp(Matrix *a, Matrix *q);
//QR方法，A最终变为对角线一二阶子块的分块上三角阵
void QR(Matrix *a, Matrix *q);
//矩阵拷贝
void MatrixCopy(Matrix A, Matrix *a);
//双步位移QR，对M进行分解，并迭代AK
void QR2Tran_QR(Matrix *b, Matrix *a);
//双步位移QR方法中，解二次方程，并存储两个特征值
void QR2Tran_solve(double x11, double x12, double x21, double x22, complex *lambda);
//双步位移QR方法中，求全部特征值，自带拟上三角化
void QR2Tran(Matrix A, complex *lambda);
//Newton法，解非线性方程
void NewtonMethod(double x, double y, double TU[2]);
//矩阵转置,B = AT
void Matrix_T(Matrix A, Matrix *b);
//向量Y=方阵A*向量X (m,n)*(n,1)=(m,1)
void MtV(Matrix A, Vector X, Vector *y);
//方阵C=矩阵A*矩阵B (m,n)*(n,m)=(m,m)
void mmM(Matrix A, Matrix B, Matrix *c);
//矩阵C=矩阵A*矩阵B (m,k)*(k,n)=(m,n)
void mmm(Matrix A, Matrix B, Matrix *c);
//解矩阵方程AXB (m,m)*(m,n)=(m,n)
void Gauss_matrix(Matrix A, Matrix B, Matrix *x);
//插值(t,u)->z
double Interpotation(double t, double u, double T[6], double U[6], double Z[6][6]);
//z=f(x,y)
double func(double x, double y, double T[6], double U[6], double Z[6][6]);
//求解满足精度要求的拟合次数，计算并输出正交基下的拟合系数矩阵，返回p*ij
int fitting(double f[11][21], double p_star[9][6]);
#endif