#ifndef __NN_H__
#define __NN_H__

#include <stdint.h>

#if defined (SECP128R1) || defined (SECP128R2)
#define KEY_BIT_LEN 128
#else 
#if defined (SECP160K1) || defined (SECP160R1) || defined (SECP160R2) 
#define KEY_BIT_LEN 160
#else 
#if defined (SECP192K1) || defined (SECP192R1)
#define KEY_BIT_LEN 192
#else
#define KEY_BIT_LEN 128
#endif /* 192 */
#endif /* 160 */
#endif /* 128 */
/*---------------------------------------------------------------------------*/

#ifdef EIGHT_BIT_PROCESSOR

/* Type definitions */
typedef uint8_t NN_DIGIT;
typedef uint16_t NN_DOUBLE_DIGIT;

/* Types for length */
typedef uint8_t NN_UINT;
typedef uint16_t NN_UINT2;

/* Length of digit in bits */
#define NN_DIGIT_BITS 8

/* Length of digit in bytes */
#define NN_DIGIT_LEN (NN_DIGIT_BITS/8)

/* Maximum value of digit */
#define MAX_NN_DIGIT 0xff

/* Number of digits in key
 * used by optimized mod multiplication (ModMultOpt) and optimized mod square (ModSqrOpt)
 *
 */
#define KEYDIGITS (KEY_BIT_LEN/NN_DIGIT_BITS) //20

/* Maximum length in digits */
#define MAX_NN_DIGITS (KEYDIGITS+1)

/* buffer size
 *should be large enough to hold order of base point
 */
#define NUMWORDS MAX_NN_DIGITS

#define MOD_SQR_MASK1 0x8000
#define MOD_SQR_MASK2 0x0100

#endif /* EIGHT_BIT_PROCESSOR */

/*---------------------------------------------------------------------------*/
#ifdef SIXTEEN_BIT_PROCESSOR

/* Type definitions */
typedef uint16_t NN_DIGIT;
typedef uint32_t NN_DOUBLE_DIGIT;

/* Types for length */
typedef uint8_t NN_UINT;
typedef uint16_t NN_UINT2;

/* Length of digit in bits */
#define NN_DIGIT_BITS 16

/* Length of digit in bytes */
#define NN_DIGIT_LEN (NN_DIGIT_BITS/8)

/* Maximum value of digit */
#define MAX_NN_DIGIT 0xffff

/* Number of digits in key
 * used by optimized mod multiplication (ModMultOpt) and optimized mod square (ModSqrOpt)
 *
 */
#define KEYDIGITS (KEY_BIT_LEN/NN_DIGIT_BITS) //10

/* Maximum length in digits */
#define MAX_NN_DIGITS (KEYDIGITS+1)

/* buffer size
 *should be large enough to hold order of base point
 */
#define NUMWORDS MAX_NN_DIGITS

#define MOD_SQR_MASK1 0x80000000
#define MOD_SQR_MASK2 0x00010000

#endif /* SIXTEEN_BIT_PROCESSOR */
/*---------------------------------------------------------------------------*/

#ifdef THIRTYTWO_BIT_PROCESSOR

/* Type definitions */
typedef uint32_t NN_DIGIT;
typedef uint64_t NN_DOUBLE_DIGIT;

/* Types for length */
typedef uint8_t NN_UINT;
typedef uint16_t NN_UINT2;

/* Length of digit in bits */
#define NN_DIGIT_BITS 32

/* Length of digit in bytes */
#define NN_DIGIT_LEN (NN_DIGIT_BITS/8)

/* Maximum value of digit */
#define MAX_NN_DIGIT 0xffffffff

/* Number of digits in key
 * used by optimized mod multiplication (ModMultOpt) and optimized mod square (ModSqrOpt)
 *
 */
#define KEYDIGITS (KEY_BIT_LEN/NN_DIGIT_BITS) //5

/* Maximum length in digits */
#define MAX_NN_DIGITS (KEYDIGITS+1)

/* buffer size
 *should be large enough to hold order of base point
 */
#define NUMWORDS MAX_NN_DIGITS

/* the mask for ModSqrOpt */
#define MOD_SQR_MASK1 0x8000000000000000ll
#define MOD_SQR_MASK2 0x0000000100000000ll

#endif /* THIRTYTWO_BIT_PROCESSOR */


/*---------------------------------------------------------------------------*/
/* Conversion functions */

  /*
    CONVERSIONS
   NN_Decode (a, digits, b, len)   Decodes character string b into a.
   NN_Encode (a, len, b, digits)   Encodes a into character string b.

   ASSIGNMENTS
   NN_Assign (a, b, digits)        Assigns a = b.
   NN_ASSIGN_DIGIT (a, b, digits)  Assigns a = b, where b is a digit.
   NN_AssignZero (a, digits)    Assigns a = 0.
   NN_Assign2Exp (a, b, digits)    Assigns a = 2^b.
     
   ARITHMETIC OPERATIONS
   NN_Add (a, b, c, digits)        Computes a = b + c.
   NN_Sub (a, b, c, digits)        Computes a = b - c.
   NN_Mult (a, b, c, digits)       Computes a = b * c.
   NN_LShift (a, b, c, digits)     Computes a = b * 2^c.
   NN_RShift (a, b, c, digits)     Computes a = b / 2^c.
   NN_Div (a, b, c, cDigits, d, dDigits)  Computes a = c div d and b = c mod d.

   NUMBER THEORY
   NN_Mod (a, b, bDigits, c, cDigits)  Computes a = b mod c.
   NN_ModMult (a, b, c, d, digits) Computes a = b * c mod d.
   NN_ModExp (a, b, c, cDigits, d, dDigits)  Computes a = b^c mod d.
   NN_ModInv (a, b, c, digits)     Computes a = 1/b mod c.
   NN_Gcd (a, b, c, digits)        Computes a = gcd (b, c).

   OTHER OPERATIONS
   NN_EVEN (a, digits)             Returns 1 iff a is even.
   NN_Cmp (a, b, digits)           Returns sign of a - b.
   NN_EQUAL (a, digits)            Returns 1 iff a = b.
   NN_Zero (a, digits)             Returns 1 iff a = 0.
   NN_Digits (a, digits)           Returns significant length of a in digits.
   NN_Bits (a, digits)             Returns significant length of a in bits.
  */


/* Sets a = b / c, where a and c are digits.
 *
 * Lengths: b[2].
 * Assumes b[1] < c and HIGH_HALF (c) > 0. For efficiency, c should be
 * normalized.
 */
void NN_DigitDiv(NN_DIGIT *a, NN_DIGIT b[2], NN_DIGIT c);

void NN_Decode(NN_DIGIT *a, NN_UINT digits, unsigned char *b, NN_UINT len);

void NN_Encode(unsigned char *a, NN_UINT len, NN_DIGIT *b, NN_UINT digits);

void NN_Assign(NN_DIGIT *a, NN_DIGIT *b, NN_UINT digits);

void NN_AssignZero(NN_DIGIT *a, NN_UINT digits);

void NN_Assign2Exp(NN_DIGIT *a, NN_UINT2 b, NN_UINT digits);

NN_DIGIT NN_Add(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *c, NN_UINT digits);

NN_DIGIT NN_Sub(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *c, NN_UINT digits);

void NN_Mult(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *c, NN_UINT digits);

NN_DIGIT NN_LShift(NN_DIGIT *a, NN_DIGIT *b, NN_UINT c, NN_UINT digits);

NN_DIGIT NN_RShift(NN_DIGIT *a, NN_DIGIT *b, NN_UINT c, NN_UINT digits);

void NN_Div(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *c, NN_UINT cDigits, NN_DIGIT *d, NN_UINT dDigits);

void NN_Mod(NN_DIGIT *a, NN_DIGIT *b, NN_UINT bDigits, NN_DIGIT *c, NN_UINT cDigits);

void NN_ModAdd(NN_DIGIT * a, NN_DIGIT * b, NN_DIGIT * c, NN_DIGIT * d, NN_UINT digits);

void NN_ModMult(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *c, NN_DIGIT *d, NN_UINT digits);

void NN_ModExp(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *c, NN_UINT cDigits, NN_DIGIT *d, NN_UINT dDigits);

void NN_ModInv(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *c, NN_UINT digits);

void NN_Gcd(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *c, NN_UINT digits);

int NN_Cmp(NN_DIGIT *a, NN_DIGIT *b, NN_UINT digits);

int NN_Equal(NN_DIGIT *a, NN_DIGIT *b, NN_UINT digits);

int NN_Zero(NN_DIGIT *a, NN_UINT digits);

unsigned int NN_Bits(NN_DIGIT *a, NN_UINT digits);

unsigned int NN_Digits(NN_DIGIT *a, NN_UINT digits);

NN_DIGIT NN_AddDigitMult(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT c, NN_DIGIT *d, NN_UINT digits);

NN_DIGIT NN_SubDigitMult(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT c, NN_DIGIT *d, NN_UINT digits);

unsigned int NN_DigitBits(NN_DIGIT a);

void NN_ModMultOpt(NN_DIGIT * a, NN_DIGIT * b, NN_DIGIT * c, NN_DIGIT * d, NN_DIGIT * omega, NN_UINT digits);

void NN_ModSqrOpt(NN_DIGIT * a, NN_DIGIT * b, NN_DIGIT * d, NN_DIGIT * omega, NN_UINT digits);

void NN_ModSub(NN_DIGIT * a, NN_DIGIT * b, NN_DIGIT * c, NN_DIGIT * d, NN_UINT digits);

void NN_ModSmall(NN_DIGIT * b, NN_DIGIT * c, NN_UINT digits);

void NN_AssignDigit(NN_DIGIT * a, NN_DIGIT b, NN_UINT digits);

NN_UINT omega_mul(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *omega, NN_UINT digits);

//Edits

int NN_ModSqrRootOpt(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *c, NN_UINT digits, NN_DIGIT *omega);

void Lucas_Sequence(NN_DIGIT *V0, NN_DIGIT *Q0, NN_DIGIT *P, NN_DIGIT *Q, NN_DIGIT *k, NN_DIGIT *p, NN_DIGIT *omega);

void NN_ModSqr(NN_DIGIT * a, NN_DIGIT * b, NN_DIGIT * d, NN_UINT digits);

NN_DIGIT b_testbit(NN_DIGIT * a, int16_t i);


#endif /* __NN_H__ */
