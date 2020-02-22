#ifndef __ECC_H__
#define __ECC_H__

#include "nn.h"

//The size of sliding window, must be power of 2 (change this if you
//want to use other window size, for example: 2 or 4)
#define W_BITS 4

//basic mask used to generate mask array (you need to change this if
//you want to change window size)
//For example: if W_BITS is 2, BASIC_MASK must be 0x03;
//             if W_BITS is 4, BASIC_MASK must be 0x0f
//			   if W_BITS is 8, BASIC_MASK must be 0xff
#define BASIC_MASK 0x0f

//number of windows in one digit, NUM_MASKS = NN_DIGIT_BITS/W_BITS
#define NUM_MASKS (NN_DIGIT_BITS/W_BITS)

//number of points for precomputed points, NN_POINTS = 2^W_BITS - 1
#define NUM_POINTS ((1 << W_BITS) - 1)

//voltage (in mV) and frequency (in MHz) of the imote2 processor
#define CORE_VOLT 950
#define CORE_FREQ 104

/*
 * The data structure define the elliptic curve.
 */
typedef struct Curve
{
    // curve's coefficients
    NN_DIGIT a[NUMWORDS];
    NN_DIGIT b[NUMWORDS];

    //whether a is -3
    char a_minus3;

    //whether a is zero
    char a_zero;

} Curve;

typedef struct Point
{
    // point's coordinates
    NN_DIGIT x[NUMWORDS];
    NN_DIGIT y[NUMWORDS];
} Point;

typedef struct ZCoordinate
{
  NN_DIGIT z[NUMWORDS];
} ZCoordinate;


/*
 * All the parameters needed for elliptic curve operation
 */
typedef struct Params
{
    // prime modulus
    NN_DIGIT p[NUMWORDS];

    // Omega, p = 2^m -omega
    NN_DIGIT omega[NUMWORDS];

    // curve over which ECC will be performed
    Curve E;

    // base point, a point on E of order r
    Point G;

    // a positive, prime integer dividing the number of points on E
    NN_DIGIT r[NUMWORDS];

    // a positive prime integer, s.t. k = #E/r
//    NN_DIGIT k[NUMWORDS];
} Params;



/*
 * Init the parameters and base point array for sliding window method
 * the first function to call
 */
void ecc_init();
    
/*
 * Provide order of curve for the modules which need to know
 */
void ecc_get_order(NN_DIGIT * order);
    
/*
 * Point addition, P0 = P1 + P2
 */
void ecc_add(Point * P0, Point * P1, Point * P2);

/*
 * Projective point addition, P0 = P1 + P2
 */
void ecc_add_proj(Point * P0, NN_DIGIT *Z0, Point * P1, NN_DIGIT * Z1, Point * P2, NN_DIGIT * Z2);

/*
 * Projective point doubleing, P0 = 2*P1
 */
void ecc_dbl_proj(Point * P0, NN_DIGIT *Z0, Point * P1, NN_DIGIT * Z1);
    
/*
 * Scalar point multiplication P0 = n * P1
 */
void ecc_mul(Point * P0, Point * P1, NN_DIGIT * n);
    
/*
 * Precompute the points for sliding window method
 */ 
void ecc_win_precompute(Point * baseP, Point * pointArray);

/*    
 * Scalr point multiplication using slide window method
 * P0 = n * Point, this Point may not be the base point of curve
 * pointArray is constructed by call win_precompute(Point, pointArray)
 */
void ecc_win_mul(Point * P0, NN_DIGIT * n, Point * pointArray);

/*    
 * Scalr point multiplication using slide window method, P0 = n * basePoint of curve
 */
void ecc_win_mul_base(Point * P0, NN_DIGIT * n);

//----------------------------------------------------------------Edits----------------------------------//

//int p_iszero(Point * P0);
int p_iszero(Point * P0);

int ecc_point2octet(uint8_t *octet, NN_UINT octet_len, Point *P, int compress);

void ecc_win_precompute_Z(Point * baseP, Point * pointArray, ZCoordinate *ZList);

void ecc_win_mul_Z(Point *P0, NN_DIGIT *n, Point *pointArray, ZCoordinate *ZList);

void c_m_dbl_projective(Point * P0, NN_DIGIT *Z0, uint8_t m);

int ecc_octet2point(Point *P, uint8_t *octet, int octet_len);

int ecc_check_point(Point *P);

void ecc_gen_public_key(Point *PublicKey, NN_DIGIT *PrivateKey);

void ecc_gen_private_key(NN_DIGIT *PrivateKey);
//----------------------------------------------------------------Edits----------------------------------//

/*
 * Get base point
 */
Point * ecc_get_base_p();

/* 
 * Get parameters
 */

Params * ecc_get_param();


#endif /* __ECC_H__ */
