#include "curve_param.h"
#include "string.h"

#define TRUE 1
#define FALSE 0

void 
get_curve_param(Params *para)
{
  //init parameters
  //prime
  para->p[4] = 0x00000000;
  para->p[3] = 0xFFFFFFFD;
  para->p[2] = 0xFFFFFFFF;
  para->p[1] = 0xFFFFFFFF;
  para->p[0] = 0xFFFFFFFF;
    
  memset(para->omega, 0, NUMWORDS * NN_DIGIT_LEN);
  para->omega[0] = 0x00000001;
  para->omega[3] = 0x00000002;	  
  //cure that will be used
  //a
  para->E.a[4] = 0x00000000;
  para->E.a[3] = 0xFFFFFFFD;
  para->E.a[2] = 0xFFFFFFFF;
  para->E.a[1] = 0xFFFFFFFF;
  para->E.a[0] = 0xFFFFFFFC;
   
  para->E.a_minus3 = TRUE;
  para->E.a_zero = FALSE;
   
  //b
  para->E.b[4] = 0x00000000;
  para->E.b[3] = 0xE87579C1;
  para->E.b[2] = 0x1079F43D;
  para->E.b[1] = 0xD824993C;
  para->E.b[0] = 0x2CEE5ED3;
          
  //base point
  para->G.x[4] =  0x00000000;
  para->G.x[3] =  0x161FF752;
  para->G.x[2] =  0x8B899B2D;
  para->G.x[1] =  0x0C28607C;
  para->G.x[0] =  0xA52C5B86;
   
  para->G.y[4] =  0x00000000;
  para->G.y[3] =  0xCF5AC839;
  para->G.y[2] =  0x5BAFEB13;
  para->G.y[1] =  0xC02DA292;
  para->G.y[0] =  0xDDED7A83;
          	
  //prime divide the number of points
  para->r[4] = 0x00000000;
  para->r[3] = 0xFFFFFFFE;
  para->r[2] = 0x00000000;
  para->r[1] = 0x75A30D1B;
  para->r[0] = 0x9038A115; 

}

NN_UINT 
omega_mul(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *omega, NN_UINT digits)
{
  NN_Assign(a, b, digits);
  a[digits + 3] += NN_AddDigitMult(&a[3], &a[3], omega[3], b, digits);
  return (digits + 4);
}

