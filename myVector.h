#ifndef _MYVECTOR_H_
#define _MYVECTOR_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ElemType double
#define ERROR -2
#define MISMATCH -3  //(维数)不匹配
#define UNORDERED -4 //不能用顺序法求解
#define STRANGE 5    //矩阵奇异
#define E 1.0e-12    //相对误差限
typedef struct
{
    int dimension;  //维数
    ElemType *elem; //向量,一维数组
} Vector;

typedef struct
{
    int dimension; //矩阵维数
    int s;         //上半带宽
    int r;         //下半带宽
    double **elem; //矩阵，二维数组
} Matrix_diag;     //n维s+1+r对角方阵
void VectorInit(Vector *v, int n);
void VectorCopy(Vector V, Vector *x);
void PrintVector(Vector V);
#endif