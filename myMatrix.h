#ifndef _MY_MATRIX_H_
#define _MY_MATRIX_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myVector.h"

typedef struct
{
    int dimension; //矩阵维数
    int s;         //上半带宽
    int r;         //下半带宽
    double **elem; //矩阵，二维数组
} Matrix_diag;     //n维s+1+r对角方阵

//输入&M，n：初始化n维s+1+r对角方阵M(用(1~s+1+r)*(1~n)):
void Matrix_diagInit(Matrix_diag *a, int n, int s, int r);
//输出带状方阵A:
void PrintMatrix_diag(Matrix_diag A);
//把A值复制给A1（A为带状阵，A1不需要初始化）
void Matrix_diagCopy(Matrix_diag A, Matrix_diag *a1);
//矩阵平移(A->A-pI)
void Trans(Matrix_diag *a, double p);
//幂法(求矩阵按模最大的特征值)
double PM_eigenvalue(Matrix_diag A);
#endif