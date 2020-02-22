#include "nn.h"
#include <string.h>

#define TRUE 1
#define FALSE 0

#define MODINVOPT

#define MAX(a,b) ((a) < (b) ? (b) : (a))
#define DIGIT_MSB(x) (NN_DIGIT)(((x) >> (NN_DIGIT_BITS - 1)) & 1)
#define DIGIT_2MSB(x) (NN_DIGIT)(((x) >> (NN_DIGIT_BITS - 2)) & 3)

#define NN_ASSIGN_DIGIT(a, b, digits) {NN_AssignZero (a, digits); a[0] = b;}
#define NN_EQUAL(a, b, digits) (! NN_Cmp (a, b, digits))
#define NN_EVEN(a, digits) (((digits) == 0) || ! (a[0] & 1))

/*---------------------------------------------------------------------------*/
/* 
 * Return b * c, where b and c are NN_DIGITs. 
 * The result is a NN_DOUBLE_DIGIT
 */

#define NN_DigitMult(b, c) (NN_DOUBLE_DIGIT)(b) * (c)

/*---------------------------------------------------------------------------*/
void 
NN_DigitDiv(NN_DIGIT *a, NN_DIGIT b[2], NN_DIGIT c)
{
  NN_DOUBLE_DIGIT t;
  t = (((NN_DOUBLE_DIGIT)b[1]) << NN_DIGIT_BITS) ^ ((NN_DOUBLE_DIGIT)b[0]);
  *a = t/c;
}
/*---------------------------------------------------------------------------*/

/* Decodes character string b into a, where character string is ordered
 * from most to least significant.
 *
 * Lengths: a[digits], b[len].
 * Assumes b[i] = 0 for i < len - digits * NN_DIGIT_LEN. (Otherwise most
 * significant bytes are truncated.)
 */
void 
NN_Decode(NN_DIGIT *a, NN_UINT digits, unsigned char *b, NN_UINT len)
{
  NN_DIGIT t;
  int j;
  unsigned int i, u;
  
  for(i = 0, j = len - 1; i < digits && j >= 0; i++) {
    t = 0;
    for(u = 0; j >= 0 && u < NN_DIGIT_BITS; j--, u += 8) {
      t |= ((NN_DIGIT)b[j]) << u;
    }
    a[i] = t;
  }
  
  for(; i < digits; i++) {
    a[i] = 0;
  }
}
/*---------------------------------------------------------------------------*/

/* Encodes b into character string a, where character string is ordered
 * from most to least significant.
 *
 * Lengths: a[len], b[digits].
 * Assumes NN_Bits (b, digits) <= 8 * len. (Otherwise most significant
 * digits are truncated.)
 */
void 
NN_Encode(unsigned char *a, NN_UINT len, NN_DIGIT *b, NN_UINT digits)
{
  NN_DIGIT t;
  int j;
  unsigned int i, u;

  for(i = 0, j = len - 1; i < digits && j >= 0; i++) {
    t = b[i];
    for(u = 0; j >= 0 && u < NN_DIGIT_BITS; j--, u += 8) {
      a[j] = (unsigned char)(t >> u);
    }
  }

  for(; j >= 0; j--) {
    a[j] = 0;
  }
}

/*---------------------------------------------------------------------------*/
/* Assigns a = b.
 * Lengths: a[digits], b[digits].
 */
void 
NN_Assign(NN_DIGIT *a, NN_DIGIT *b, NN_UINT digits)
{
  memcpy(a, b, digits*NN_DIGIT_LEN);
}

/*---------------------------------------------------------------------------*/
/* Assigns a = 0.
 * Lengths: a[digits].
 */
void 
NN_AssignZero(NN_DIGIT *a, NN_UINT digits)
{
  uint8_t i;

  for(i = 0; i < digits; i++) {
    a[i] = 0;
  }
}

/*---------------------------------------------------------------------------*/
/* Assigns a = 2^b.
 * Lengths: a[digits].
 *  Requires b < digits * NN_DIGIT_BITS.
 */
void 
NN_Assign2Exp(NN_DIGIT *a, NN_UINT2 b, NN_UINT digits)
{
  NN_AssignZero(a, digits);

  if(b >= digits * NN_DIGIT_BITS) {
    return;
  }

  a[b / NN_DIGIT_BITS] = (NN_DIGIT)1 << (b % NN_DIGIT_BITS);
}

/*---------------------------------------------------------------------------*/

/* Computes a = b + c. Returns carry.
 * a, b ,c can be same
 * Lengths: a[digits], b[digits], c[digits].
 */
NN_DIGIT 
NN_Add(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *c, NN_UINT digits)
{
  NN_DIGIT carry, ai;
  NN_UINT i;

  carry = 0;

  for(i = 0; i < digits; i++) {
    if((ai = b[i] + carry) < carry) {
      ai = c[i];
    } else if ((ai += c[i]) < c[i]) {
      carry = 1;
    } else {
      carry = 0;
    }
    a[i] = ai;
  }

  return carry;
}

/*---------------------------------------------------------------------------*/
/* Computes a = b - c. Returns borrow.
 * a, b, c can be same
 * Lengths: a[digits], b[digits], c[digits].
 */
NN_DIGIT 
NN_Sub(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *c, NN_UINT digits)
{
  NN_DIGIT ai, borrow;
  NN_UINT i;

  borrow = 0;

  for(i = 0; i < digits; i++) {
    if((ai = b[i] - borrow) > (MAX_NN_DIGIT - borrow)) {
      ai = MAX_NN_DIGIT - c[i];
    } else if((ai -= c[i]) > (MAX_NN_DIGIT - c[i])) {
        borrow = 1;
    } else {
        borrow = 0;
    }
    a[i] = ai;
  }

  return borrow;
}
/*---------------------------------------------------------------------------*/
/* Computes a = b * c.
 * a, b, c can be same
 * Lengths: a[2*digits], b[digits], c[digits].
 * Assumes digits < MAX_NN_DIGITS.
 */
void 
NN_Mult(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *c, NN_UINT digits)
{
  NN_DIGIT t[2 * MAX_NN_DIGITS];
  unsigned int b_digits, c_digits, i;

  NN_AssignZero (t, 2 * digits);
  
  b_digits = NN_Digits (b, digits);
  c_digits = NN_Digits (c, digits);

  for (i = 0; i < b_digits; i++)
    t[i + c_digits] += NN_AddDigitMult (&t[i], &t[i], b[i], c, c_digits);
  
  NN_Assign (a, t, 2 * digits);

}
/*---------------------------------------------------------------------------*/

void 
NN_Sqr(NN_DIGIT *a, NN_DIGIT *b, NN_UINT digits)
{
  NN_DIGIT t[2 * MAX_NN_DIGITS];
  NN_UINT b_digits, i;

  NN_AssignZero (t, 2 * digits);
  
  b_digits = NN_Digits (b, digits);
    
  for (i = 0; i < b_digits; i++) {
    t[i + b_digits] += NN_AddDigitMult (&t[i], &t[i], b[i], b, b_digits);
  }
  
  NN_Assign (a, t, 2 * digits);
}

/*---------------------------------------------------------------------------*/

/* Computes a = b * 2^c (i.e., shifts left c bits), returning carry.
 * a, b can be same
 * Lengths: a[digits], b[digits].
 * Requires c < NN_DIGIT_BITS.
 */
NN_DIGIT 
NN_LShift(NN_DIGIT *a, NN_DIGIT *b, NN_UINT c, NN_UINT digits)
{
  NN_DIGIT bi, carry;
  NN_UINT i, t;
  
  if(c >= NN_DIGIT_BITS) {
    return (0);
  }
  
  t = NN_DIGIT_BITS - c;

  carry = 0;

  for(i = 0; i < digits; i++) {
    bi = b[i];
    a[i] = (bi << c) | carry;
    carry = c ? (bi >> t) : 0;
  }
  
  return carry;
}
/*---------------------------------------------------------------------------*/
/* Computes a = b div 2^c (i.e., shifts right c bits), returning carry.
 * a, b can be same
 * Lengths: a[digits], b[digits].
 * Requires: c < NN_DIGIT_BITS.
 */
NN_DIGIT 
NN_RShift(NN_DIGIT *a, NN_DIGIT *b, NN_UINT c, NN_UINT digits)
{
  NN_DIGIT bi, carry;
  int i;
  NN_UINT t;
  
  if(c >= NN_DIGIT_BITS) {
    return (0);
  }
  
  t = NN_DIGIT_BITS - c;

  carry = 0;

  for(i = digits - 1; i >= 0; i--) {
    bi = b[i];
    a[i] = (bi >> c) | carry;
    carry = c ? (bi << t) : 0;
  }
  
  return carry;
}

/*---------------------------------------------------------------------------*/
/* Computes a = c div d and b = c mod d.
 * a, c, d can be same
 * b, c, d can be same
 * Lengths: a[cDigits], b[dDigits], c[cDigits], d[dDigits].
 * Assumes d > 0, cDigits < 2 * MAX_NN_DIGITS,
 *         dDigits < MAX_NN_DIGITS.
 */
void 
NN_Div(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *c, NN_UINT c_digits, NN_DIGIT *d, NN_UINT d_digits)
{
  NN_DIGIT ai, cc[2 * MAX_NN_DIGITS+1], dd[MAX_NN_DIGITS], t;

  int i;
  int dd_digits, shift;
  
  dd_digits = NN_Digits (d, d_digits);
  if(dd_digits == 0) {
    return;
  }
  
  /* Normalize operands. */
  shift = NN_DIGIT_BITS - NN_DigitBits (d[dd_digits-1]);
  NN_AssignZero (cc, dd_digits);
  cc[c_digits] = NN_LShift (cc, c, shift, c_digits);
  NN_LShift (dd, d, shift, dd_digits);
  t = dd[dd_digits - 1];

  if(a != NULL) {
    NN_AssignZero (a, c_digits);
  }

  for(i = c_digits - dd_digits; i >= 0; i--) {
    /* Underestimate quotient digit and subtract. */
    if (t == MAX_NN_DIGIT) {
      ai = cc[i + dd_digits];
    } else {
      NN_DigitDiv(&ai, &cc[i + dd_digits-1], t + 1);
    }
    cc[i + dd_digits] -= NN_SubDigitMult(&cc[i], &cc[i], ai, dd, dd_digits);

    /* Correct estimate. */
    while(cc[i + dd_digits] || (NN_Cmp (&cc[i], dd, dd_digits) >= 0)) {
      ai++;
      cc[i+dd_digits] -= NN_Sub (&cc[i], &cc[i], dd, dd_digits);
    }
    if(a != NULL) {
      a[i] = ai;
    }
  }    
    /* Restore result. */
    NN_AssignZero (b, d_digits);
    NN_RShift (b, cc, shift, dd_digits);
  //}
}
/*---------------------------------------------------------------------------*/

/* Computes a = b mod c.
 * Lengths: a[cDigits], b[bDigits], c[cDigits].
 * Assumes c > 0, bDigits < 2 * MAX_NN_DIGITS, cDigits < MAX_NN_DIGITS.
 */

void 
NN_Mod(NN_DIGIT *a, NN_DIGIT *b, NN_UINT b_digits, NN_DIGIT *c, NN_UINT c_digits)
{  
  NN_Div(NULL, a, b, b_digits, c, c_digits);
}

/*---------------------------------------------------------------------------*/
/* Computes a = b * c mod d.
 * a, b, c can be same
 * Lengths: a[digits], b[digits], c[digits], d[digits].
 * Assumes d > 0, digits < MAX_NN_DIGITS.
 */
void 
NN_ModMult(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *c, NN_DIGIT *d, NN_UINT digits)
{
  NN_DIGIT t[2 * MAX_NN_DIGITS];
    
  //memset(t, 0, 2*MAX_NN_DIGITS*NN_DIGIT_LEN);
  t[2 * MAX_NN_DIGITS-1] = 0;
  t[2 * MAX_NN_DIGITS-2] = 0;
  NN_Mult(t, b, c, digits);
  NN_Mod(a, t, 2 * digits, d, digits);
}

/*---------------------------------------------------------------------------*/
/* Computes a = b^c mod d.
 * Lengths: a[dDigits], b[dDigits], c[cDigits], d[dDigits].
 * Assumes d > 0, cDigits > 0, dDigits < MAX_NN_DIGITS.
 */
void 
NN_ModExp(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *c, NN_UINT c_digits, NN_DIGIT *d, NN_UINT d_digits)
{
  NN_DIGIT b_power[3][MAX_NN_DIGITS], ci, t[MAX_NN_DIGITS];
  int i;
  uint8_t ci_bits, j, s;

  /* Store b, b^2 mod d, and b^3 mod d. */

  NN_Assign(b_power[0], b, d_digits);
  NN_ModMult(b_power[1], b_power[0], b, d, d_digits);
  NN_ModMult(b_power[2], b_power[1], b, d, d_digits);
  
  NN_ASSIGN_DIGIT(t, 1, d_digits);

  c_digits = NN_Digits(c, c_digits);

  for(i = c_digits - 1; i >= 0; i--) {
    ci = c[i];
    ci_bits = NN_DIGIT_BITS;
      
    /* Scan past leading zero bits of most significant digit. */
    if(i == (int)(c_digits - 1)) {
      while (! DIGIT_2MSB (ci)) {
        ci <<= 2;
        ci_bits -= 2;
      }
    }

    for(j = 0; j < ci_bits; j += 2, ci <<= 2) {
      /* Compute t = t^4 * b^s mod d, where s = two MSB's of ci. */
      NN_ModMult(t, t, t, d, d_digits);
      NN_ModMult(t, t, t, d, d_digits);
      if ((s = DIGIT_2MSB (ci)) != 0) {
        NN_ModMult (t, t, b_power[s-1], d, d_digits);
      }
    }
  }
  
  NN_Assign(a, t, d_digits);
}

/*---------------------------------------------------------------------------*/
/* Compute a = 1/b mod c, assuming inverse exists.
 * a, b, c can be same
 * Lengths: a[digits], b[digits], c[digits].
 * Assumes gcd (b, c) = 1, digits < MAX_NN_DIGITS.
   */
void 
NN_ModInv(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *c, NN_UINT digits)
{
  NN_DIGIT q[MAX_NN_DIGITS], t1[MAX_NN_DIGITS], t3[MAX_NN_DIGITS],
          u1[MAX_NN_DIGITS], u3[MAX_NN_DIGITS], v1[MAX_NN_DIGITS],
          v3[MAX_NN_DIGITS], w[2 * MAX_NN_DIGITS];
  int u1Sign;

  /* Apply extended Euclidean algorithm, modified to avoid negative numbers. */
  NN_ASSIGN_DIGIT(u1, 1, digits);
  NN_AssignZero(v1, digits);
  NN_Assign(u3, b, digits);
  NN_Assign(v3, c, digits);
  u1Sign = 1;

  while (!NN_Zero(v3, digits)) {
    NN_Div (q, t3, u3, digits, v3, digits);
    NN_Mult (w, q, v1, digits);
    NN_Add (t1, u1, w, digits);
    NN_Assign (u1, v1, digits);
    NN_Assign (v1, t1, digits);
    NN_Assign (u3, v3, digits);
    NN_Assign (v3, t3, digits);
    u1Sign = -u1Sign;
  }
  
  /* Negate result if sign is negative. */
  if (u1Sign < 0) {
      NN_Sub (a, c, u1, digits);
  } else {
      NN_Assign (a, u1, digits);
  }

}

/*---------------------------------------------------------------------------*/
/*
 * a= b/c mod d
 * algorithm in "From Euclid's GCD to Montgomery Multiplication to the Great Divide"
 * 
 */
void 
NN_ModDivOpt(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *c, NN_DIGIT *d, NN_UINT digits)
{
  NN_DIGIT A[MAX_NN_DIGITS], B[MAX_NN_DIGITS], U[MAX_NN_DIGITS], V[MAX_NN_DIGITS];
  int tmp_even;

  NN_Assign(A, c, digits);
  NN_Assign(B, d, digits);
  NN_Assign(U, b, digits);
  NN_AssignZero(V, digits);
    
  while((tmp_even = NN_Cmp(A, B, digits)) != 0) {
    if(NN_EVEN(A, digits)) {
      NN_RShift(A, A, 1, digits);
      if(NN_EVEN(U, digits)){
        NN_RShift(U, U, 1, digits);
      } else {
        NN_Add(U, U, d, digits);
        NN_RShift(U, U, 1, digits);
      }
    } else if(NN_EVEN(B, digits)) {
      NN_RShift(B, B, 1, digits);
      if(NN_EVEN(V, digits)) {
        NN_RShift(V, V, 1, digits);
      } else {
        NN_Add(V, V, d, digits);
        NN_RShift(V, V, 1, digits);
      }
    } else if(tmp_even > 0) {
      NN_Sub(A, A, B, digits);
      NN_RShift(A, A, 1, digits);
      if(NN_Cmp(U, V, digits) < 0) {
        NN_Add(U, U, d, digits);
      }
      NN_Sub(U, U, V, digits);
      if(NN_EVEN(U, digits)) {
        NN_RShift(U, U, 1, digits);
      } else {
        NN_Add(U, U, d, digits);
        NN_RShift(U, U, 1, digits);
      }
    } else {
      NN_Sub(B, B, A, digits);
      NN_RShift(B, B, 1, digits);
      if(NN_Cmp(V, U, digits) < 0) {
        NN_Add(V, V, d, digits);
      }
      NN_Sub(V, V, U, digits);
      if(NN_EVEN(V, digits)) {
        NN_RShift(V, V, 1, digits);
      } else {
        NN_Add(V, V, d, digits);
        NN_RShift(V, V, 1, digits);
      }
    }
  }

  NN_Assign(a, U, digits);
}

/*---------------------------------------------------------------------------*/
/* Computes a = gcd(b, c).
 * a, b, c can be same
 * Lengths: a[digits], b[digits], c[digits].
 * Assumes b > c, digits < MAX_NN_DIGITS.
 */
void 
NN_Gcd(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *c, NN_UINT digits)
{
  NN_DIGIT t[MAX_NN_DIGITS], u[MAX_NN_DIGITS], v[MAX_NN_DIGITS];

  NN_Assign(u, b, digits);
  NN_Assign(v, c, digits);

  while(!NN_Zero(v, digits)) {
    NN_Mod(t, u, digits, v, digits);
    NN_Assign(u, v, digits);
    NN_Assign(v, t, digits);
  }

  NN_Assign(a, u, digits);

}

/*---------------------------------------------------------------------------*/

/* Returns sign of a - b.
 * Lengths: a[digits], b[digits].
 */
int 
NN_Cmp(NN_DIGIT *a, NN_DIGIT *b, NN_UINT digits)
{
  int i;
  
  for(i = digits - 1; i >= 0; i--) { 
    if(a[i] > b[i]) {
      return 1;
    } else if (a[i] < b[i]) {
      return -1; 
    }
  }

  return 0;
}

/*---------------------------------------------------------------------------*/
/* Returns nonzero iff a is zero.
 *  Lengths: a[digits].
 */
int 
NN_Zero(NN_DIGIT *a, NN_UINT digits)
{
  NN_UINT i;
  
  for(i = 0; i < digits; i++) {
    if(a[i]) {
      return 0;
    }
  }
  return 1;
}
/*---------------------------------------------------------------------------*/
/* 
 * returns 1 iff a = 1
 */
int 
NN_One(NN_DIGIT * a, NN_UINT digits)
{
  uint8_t i;
    
  for(i = 1; i < digits; i++) {
    if(a[i]) {
      return FALSE;
    }
    if(a[0] == 1) {
      return TRUE;
    }
  }
    
  return FALSE;
}
/*---------------------------------------------------------------------------*/
/* Returns the significant length of a in bits.
 * Lengths: a[digits].
 */
unsigned int 
NN_Bits(NN_DIGIT *a, NN_UINT digits)
{
  if((digits = NN_Digits(a, digits)) == 0) {
    return 0;
  }
  return ((digits - 1) * NN_DIGIT_BITS + NN_DigitBits (a[digits-1]));
}

/*---------------------------------------------------------------------------*/

/* Returns the significant length of a in digits.
 * Lengths: a[digits].
 */
unsigned int 
NN_Digits(NN_DIGIT *a, NN_UINT digits)
{
  int i;
  
  for(i = digits - 1; i >= 0; i--) {
    if(a[i]) {
      break;
    }
  }

  return i + 1;
}

/*---------------------------------------------------------------------------*/
/* Computes a = b + c*d, where c is a digit. Returns carry.
 * a, b, c can be same
 * Lengths: a[digits], b[digits], d[digits].
 */
NN_DIGIT 
NN_AddDigitMult(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT c, NN_DIGIT *d, NN_UINT digits)
{
  NN_DIGIT carry;
  unsigned int i;
  NN_DOUBLE_DIGIT t;

  //Should copy b to a
  if(c == 0) {
    return (0);
  }

  carry = 0;

  for(i = 0; i < digits; i++) {
    t = NN_DigitMult (c, d[i]);
    if ((a[i] = b[i] + carry) < carry) {
      carry = 1;
    } else {
      carry = 0;
    }
    if((a[i] += (t & MAX_NN_DIGIT)) < (t & MAX_NN_DIGIT)) {
      carry++;
    }
    carry += (NN_DIGIT)(t >> NN_DIGIT_BITS);
  }

  return carry;
}

/*---------------------------------------------------------------------------*/
/* Computes a = b - c*d, where c is a digit. Returns borrow.
 * a, b, d can be same
 * Lengths: a[digits], b[digits], d[digits].
 */
NN_DIGIT 
NN_SubDigitMult(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT c, NN_DIGIT *d, NN_UINT digits)
{
  NN_DIGIT borrow;
  unsigned int i;

  NN_DOUBLE_DIGIT t;


  if (c == 0) {
    return 0;
  }

  borrow = 0;

  for(i = 0; i < digits; i++) {
    t = NN_DigitMult (c, d[i]);
    if ((a[i] = b[i] - borrow) > (MAX_NN_DIGIT - borrow)) {
      borrow = 1;
    } else {
      borrow = 0;
    }
    if ((a[i] -= (t & MAX_NN_DIGIT)) > (MAX_NN_DIGIT - (t & MAX_NN_DIGIT))) {
      borrow++;
    }
    borrow += (NN_DIGIT)(t >> NN_DIGIT_BITS);
  }
    
  return borrow;
}

/*---------------------------------------------------------------------------*/

/* Returns the significant length of a in bits, where a is a digit.
 */
unsigned int 
NN_DigitBits(NN_DIGIT a)
{
  unsigned int i;
  
  for(i = 0; i < NN_DIGIT_BITS; i++, a >>= 1) {
    if(a == 0) {
      break;
    }
  }  
  return i;
}
/*---------------------------------------------------------------------------*/
int 
NN_Equal(NN_DIGIT *a, NN_DIGIT *b, NN_UINT digits)
{
  return !NN_Cmp(a, b, digits);
}

/*---------------------------------------------------------------------------*/
/*
 * Computes a = b * c mod d, d is generalized mersenne prime, d = 2^KEYBITS - omega
 */
void 
NN_ModMultOpt(NN_DIGIT * a, NN_DIGIT * b, NN_DIGIT * c, NN_DIGIT * d, NN_DIGIT * omega, NN_UINT digits)
{
  NN_DIGIT t1[2*MAX_NN_DIGITS];
  NN_DIGIT t2[2*MAX_NN_DIGITS];
  NN_DIGIT *pt1;
  NN_UINT len_t2, len_t1;

  //memset(t1, 0, 2*MAX_NN_DIGITS*NN_DIGIT_LEN);
  //memset(t2+KEYDIGITS*NN_DIGIT_LEN, 0, (2*MAX_NN_DIGITS-KEYDIGITS)*NN_DIGIT_LEN);
  t1[2*MAX_NN_DIGITS-1]=0;
  t1[2*MAX_NN_DIGITS-2]=0;
  t2[2*MAX_NN_DIGITS-1]=0;
  t2[2*MAX_NN_DIGITS-2]=0;

  NN_Mult(t1, b, c, KEYDIGITS);
     
  pt1 = &(t1[KEYDIGITS]);
  len_t2 = 2 * KEYDIGITS;
     
    //the "Curve-Specific Optimizations" algorithm in "Comparing Elliptic Curve Cryptography and RSA on 8-bit CPUs"
  while(!NN_Zero(pt1, KEYDIGITS)) {
    memset(t2, 0, len_t2*NN_DIGIT_LEN);
    len_t2 -= KEYDIGITS;
    len_t1 = len_t2;
    len_t2 = omega_mul(t2, pt1, omega, len_t2);
    memset(pt1, 0, len_t1*NN_DIGIT_LEN);
    NN_Add(t1, t2, t1, MAX(KEYDIGITS,len_t2)+1);
  }

  while(NN_Cmp(t1, d, digits) > 0) {
    NN_Sub(t1, t1, d, digits);      
  }

  NN_Assign(a, t1, digits);
     
}
/*---------------------------------------------------------------------------*/
void 
NN_ModSqrOpt(NN_DIGIT * a, NN_DIGIT * b, NN_DIGIT * d, NN_DIGIT * omega, NN_UINT digits)
{
  NN_DIGIT t1[2*MAX_NN_DIGITS];
  NN_DIGIT t2[2*MAX_NN_DIGITS];
  NN_DIGIT *pt1;
  NN_UINT len_t1, len_t2;

  t1[2*MAX_NN_DIGITS-1]=0;
  t1[2*MAX_NN_DIGITS-2]=0;
  t2[2*MAX_NN_DIGITS-1]=0;
  t2[2*MAX_NN_DIGITS-2]=0;

  NN_Sqr(t1, b, KEYDIGITS);
     
  pt1 = &(t1[KEYDIGITS]);
  len_t2 = 2*KEYDIGITS;
  //the "Curve-Specific Optimizations" algorithm in "Comparing Elliptic Curve Cryptography and RSA on 8-bit CPUs"
  while(!NN_Zero(pt1, KEYDIGITS)) {
    memset(t2, 0, len_t2*NN_DIGIT_LEN);
    len_t2 -= KEYDIGITS;
    len_t1 = len_t2;
    len_t2 = omega_mul(t2, pt1, omega, len_t2);
    memset(pt1, 0, len_t1*NN_DIGIT_LEN);
    NN_Add(t1, t2, t1, MAX(KEYDIGITS,len_t2)+1);
  }
     
  while(NN_Cmp(t1, d, digits) > 0) {
    NN_Sub(t1, t1, d, digits);
  }
  NN_Assign (a, t1, digits);

}
/*---------------------------------------------------------------------------*/
/*
 * Computes a = (b - c) mod d.
 * Assume b and c are all smaller than d
 * always return positive value
 */
void 
NN_ModSub(NN_DIGIT * a, NN_DIGIT * b, NN_DIGIT * c, NN_DIGIT * d, NN_UINT digits)
{
  NN_DIGIT tmp[MAX_NN_DIGITS];
  NN_DIGIT borrow;
    
  borrow = NN_Sub(tmp, b, c, digits);
  if(borrow) {
    NN_Add(a, tmp, d, digits);
  } else {
    NN_Assign(a, tmp, digits);
  }
    
}
/*---------------------------------------------------------------------------*/
/*
 * Computes a = (b + c) mod d.
 * a, b, c can be same
 * Assumption: b,c is in [0, d)
 */
void 
NN_ModAdd(NN_DIGIT * a, NN_DIGIT * b, NN_DIGIT * c, NN_DIGIT * d, NN_UINT digits)
{
  NN_DIGIT tmp[MAX_NN_DIGITS];
  NN_DIGIT carry;
    
  carry = NN_Add(tmp, b, c, digits);
  if(carry) {
    NN_Sub(a, tmp, d, digits);
  } else if(NN_Cmp(tmp, d, digits) >= 0) {
      NN_Sub(a, tmp, d, digits);
  } else {
    NN_Assign(a, tmp, digits);
  }
    
}
/*---------------------------------------------------------------------------*/
void 
NN_ModSmall(NN_DIGIT * b, NN_DIGIT * c, NN_UINT digits)
{
  while(NN_Cmp(b, c, digits) > 0) {
    NN_Sub(b, b, c, digits);
  }
}
/*---------------------------------------------------------------------------*/
void 
NN_AssignDigit(NN_DIGIT * a, NN_DIGIT b, NN_UINT digits)
{
  NN_AssignZero(a, digits);
  a[0] = b;
}
/*---------------------------------------------------------------------------*/
NN_UINT 
omega_mul(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *omega, NN_UINT digits)
{
  //NN_Assign(a, b, digits);
  a[digits] += NN_AddDigitMult(&a[0], &a[0], omega[0], b, digits);
  return (digits + 1);
}
/*-----------------MYEDIT---------------------------*/

NN_DIGIT b_testbit(NN_DIGIT * a, int16_t i)
  {
      return (*(a + (i / NN_DIGIT_BITS)) & ((NN_DIGIT)1 << (i % NN_DIGIT_BITS)));
  }
  
void Lucas_Sequence(NN_DIGIT *V0, NN_DIGIT *Q0, NN_DIGIT *P, NN_DIGIT *Q, NN_DIGIT *k, NN_DIGIT *p, NN_DIGIT *omega){
    NN_DIGIT v0[NUMWORDS], v1[NUMWORDS], q0[NUMWORDS], q1[NUMWORDS];
    NN_DIGIT tmp[NUMWORDS];
    int r;

    memset(v0, 0, NUMWORDS*NN_DIGIT_LEN);
    v0[0] = 2;
    memcpy(v1, P, NUMWORDS*NN_DIGIT_LEN);
    memset(q0, 0, NUMWORDS*NN_DIGIT_LEN);
    q0[0] = 1;
    memset(q1, 0, NUMWORDS*NN_DIGIT_LEN);
    q1[1] = 1;

    r = NN_Bits(k, NUMWORDS) - 1;

    while(r >= 0){
      NN_ModMultOpt(q0, q0, q1, p, omega, NUMWORDS);
      if(b_testbit(k, r)){
	NN_ModMultOpt(q1, q0, Q, p, omega, NUMWORDS);
	NN_ModMultOpt(tmp, P, q0, p, omega, NUMWORDS);
	NN_ModMultOpt(v0, v0, v1, p, omega, NUMWORDS);
	NN_ModSub(v0, v0, tmp, p, NUMWORDS);
	NN_ModSqrOpt(v1, v1, p, omega, NUMWORDS);
	NN_ModSub(v1, v1, q1, p, NUMWORDS);
	NN_ModSub(v1, v1, q1, p, NUMWORDS);
      }else{
	memcpy(q1, q0, NUMWORDS*NN_DIGIT_LEN);
	NN_ModMultOpt(v1, v1, v0, p, omega, NUMWORDS);
	NN_ModMultOpt(tmp, P, q0, p, omega, NUMWORDS);
	NN_ModSub(v1, v1, tmp, p, NUMWORDS);
	NN_ModSqrOpt(v0, v0, p, omega, NUMWORDS);
	NN_ModSub(v0, v0, q0, p, NUMWORDS);
	NN_ModSub(v0, v0, q0, p, NUMWORDS);
      }
    }
    memcpy(V0, v0, NUMWORDS*NN_DIGIT_LEN);
    memcpy(Q0, q0, NUMWORDS*NN_DIGIT_LEN);

  }
  
  //Computes a = b^2 mod d, The Standard Squaring Algorithm in "High-Speed RSA Implementation"
void NN_ModSqr(NN_DIGIT * a, NN_DIGIT * b, NN_DIGIT * d, NN_UINT digits)
  {
#ifdef BARRETT_REDUCTION

    NN_DIGIT q2[2*MAX_NN_DIGITS+8], x[2*MAX_NN_DIGITS+8], r2[2*MAX_NN_DIGITS+8], m[MAX_NN_DIGITS+4];
    //int i; //for debug

    //b^2
    NN_Sqr(x, b, digits);
    memset(x+2*digits, 0, (2*MAX_NN_DIGITS+8-2*digits)*NN_DIGIT_LEN);

    memcpy(m, d, digits*NN_DIGIT_LEN);
    memset(m+digits, 0, (MAX_NN_DIGITS+4-digits)*NN_DIGIT_LEN);

    //q_2=q_1*mu
    NN_Mult(q2, x+pBarrett->km-1, pBarrett->mu, pBarrett->mu_len);

    //q_3*m
    NN_Mult(r2, q2+pBarrett->km+1, m, pBarrett->mu_len);
    memset(r2+pBarrett->km+1, 0, (2*MAX_NN_DIGITS+8-pBarrett->km-1)*NN_DIGIT_LEN);

    memset(x+pBarrett->km+1, 0, (2*MAX_NN_DIGITS+8-pBarrett->km-1)*NN_DIGIT_LEN);
    if (NN_Cmp(x, r2, pBarrett->km+1) < 0)
      x[pBarrett->km+1] = 1;
    NN_Sub(x, x, r2, pBarrett->km+2);
    while(NN_Cmp(x, m, digits) >= 0)
      NN_Sub(x, x, m, digits);
    memcpy(a, x, digits*NN_DIGIT_LEN);

#else
    NN_DIGIT t[2*MAX_NN_DIGITS];
    
    NN_Sqr (t, b, digits);
    NN_Mod (a, t, 2 * digits, d, digits);
#endif
  }

int NN_ModSqrRootOpt(NN_DIGIT *a, NN_DIGIT *b, NN_DIGIT *c, NN_UINT digits, NN_DIGIT *omega){
    NN_DIGIT m[NUMWORDS];
    NN_DIGIT k[NUMWORDS], residue[NUMWORDS];
    NN_DIGIT gam[NUMWORDS], i[NUMWORDS], twog[NUMWORDS];

    //case 1.
    memset(m, 0, NUMWORDS*NN_DIGIT_LEN);
    m[0] = 4;
    NN_Div(k, residue, c, NUMWORDS, m, NUMWORDS);
    if (residue[0] == 3){
      m[0] = 1;
      NN_Add(k, k, m, NUMWORDS);
      NN_ModExp(a, b, k, NUMWORDS, c, digits);
      return 1;
    }

    //case 2.
    m[0] = 8;
    NN_Div(k, residue, c, NUMWORDS, m, NUMWORDS);
    if(residue[0] == 5){
      NN_LShift(twog, b, 1, NUMWORDS);
      if(NN_Cmp(twog, c, NUMWORDS) >=0)
	NN_Sub(twog, twog, c, NUMWORDS);
      NN_ModExp(gam, twog, k, NUMWORDS, c, digits);
#ifdef CURVE_OPT
      NN_ModSqrOpt(i, gam, c, omega, digits);
      NN_ModMultOpt(i, i, twog, c, omega, digits);
      m[0] = 1;
      NN_ModSub(i, i, m, c, digits);
      NN_ModMultOpt(a, b, gam, c, omega, digits);
      NN_ModMultOpt(a, a, i, c, omega, digits);
#else
      NN_ModSqr(i, gam, c, digits);
      NN_ModMult(i, i, twog, c, digits);
      m[0] = 1;
      NN_ModSub(i, i, m, c, digits);
      NN_ModMult(a, b, gam, c, digits);
      NN_ModMult(a, a, i, c, digits);
#endif
      return 1;
    }

    //case 3.
    m[0] = 8;
    NN_Div(k, residue, c, NUMWORDS, m, NUMWORDS);
    if(residue[0] == 1){
      //Q = b
      
      //0<P<p
      memset(i, 0, NUMWORDS*NN_DIGIT_LEN);
      i[0] = 1;
      m[0] = 1;
      //k = (p+1)/2
      memcpy(k, c, NUMWORDS*NN_DIGIT_LEN);
      NN_Add(k, k, m, NUMWORDS);
      NN_RShift(k, k, 1, NUMWORDS);

      while (NN_Cmp(i, c, NUMWORDS) < 0){
	//compute V
	Lucas_Sequence(gam, twog, i, b, k, c, omega); //V = gam, Q = twog
	//z = V/2 mod p
	if ((gam[0] & 0x1) == 1)
	  NN_Add(gam, gam, c, NUMWORDS);	  
	NN_RShift(gam, gam, 1, NUMWORDS);
	//check whether z^2 mod p = b
	NN_ModSqrOpt(residue, gam, c, omega, NUMWORDS);
	if(NN_Cmp(residue, b, NUMWORDS) == 0){
	  memcpy(a, gam, NUMWORDS*NN_DIGIT_LEN);
	  return 1;
	}
	if(NN_Cmp(twog, m, NUMWORDS) > 0){
	  NN_Add(twog, twog, m, NUMWORDS);
	  if(NN_Cmp(twog, c, NUMWORDS) < 0)
	    return -1;
	}
	NN_Add(i, i, m, NUMWORDS);
      }
    }

    return -1;
  }
