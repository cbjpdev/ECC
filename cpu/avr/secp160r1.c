#include "curve_param.h"
#include "string.h"

#define TRUE 1
#define FALSE 0

void 
get_curve_param(Params *para)
{
  //init parameters
  //prime
  para->p[5] = 0x00000000;
  para->p[4] = 0xFFFFFFFF;
  para->p[3] = 0xFFFFFFFF;
  para->p[2] = 0xFFFFFFFF;
  para->p[1] = 0xFFFFFFFF;
  para->p[0] = 0x7FFFFFFF;
  memset(para->omega, 0, NUMWORDS);        
  para->omega[0] = 0x80000001;
    
  //cure that will be used
  //a
  para->E.a[5] = 0x00000000;
  para->E.a[4] = 0xFFFFFFFF;
  para->E.a[3] = 0xFFFFFFFF;
  para->E.a[2] = 0xFFFFFFFF;
  para->E.a[1] = 0xFFFFFFFF;
  para->E.a[0] = 0x7FFFFFFC;
    	
  para->E.a_minus3 = TRUE;
  para->E.a_zero = FALSE;
   
  //b
  para->E.b[5] = 0x00000000;
  para->E.b[4] = 0x1C97BEFC;
  para->E.b[3] = 0x54BD7A8B;
  para->E.b[2] = 0x65ACF89F;
  para->E.b[1] = 0x81D4D4AD;
  para->E.b[0] = 0xC565FA45;
        
  //base point
  para->G.x[5] = 0x00000000;
  para->G.x[4] = 0x4A96B568;
  para->G.x[3] = 0x8EF57328;
  para->G.x[2] = 0x46646989;
  para->G.x[1] = 0x68C38BB9;
  para->G.x[0] = 0x13CBFC82;
  para->G.y[5] = 0x00000000;
  para->G.y[4] = 0x23A62855;
  para->G.y[3] = 0x3168947D;
  para->G.y[2] = 0x59DCC912;
  para->G.y[1] = 0x04235137;
  para->G.y[0] = 0x7AC5FB32;
        	
  //prime divide the number of points
  para->r[5] = 0x00000001;
  para->r[4] = 0x00000000;
  para->r[3] = 0x00000000;
  para->r[2] = 0x0001F4C8;
  para->r[1] = 0xF927AED3;
  para->r[0] = 0xCA752257;

}

NN_UINT 
omega_mul(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *omega, NN_UINT digits)
{
  //NN_Assign(a, b, digits);
  a[digits] += NN_AddDigitMult(&a[0], &a[0], omega[0], b, digits);
  return (digits + 1);
}

