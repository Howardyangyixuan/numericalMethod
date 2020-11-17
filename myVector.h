#ifndef _MY_VECTOR_H_
#define _MY_VECTOR_H_
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

void VectorInit(Vector *v, int n);
void VectorCopy(Vector V, Vector *x);
void PrintVector(Vector V);
Vector *Vector_Num(Vector U, double k, Vector *v);
double DotProduct(Vector X, Vector Y);
void VmV(Vector X, Vector Y, Vector *v); //X - Y = V
#endif