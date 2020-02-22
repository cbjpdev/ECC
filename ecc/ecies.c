#include "ecc.h"
#include "sha1.h"
#include "lib/random.h"
#include <stdio.h>
#include <string.h>

#define MAX_M_LEN 41
#define HMAC_LEN 20

#define TRUE 1
#define FALSE 0

//#ifdef SLIDING_WIN
	static Point baseArray[NUM_POINTS];
//#ifdef PROJECTIVE
 	static ZCoordinate ZList[NUM_POINTS];
/*#endif*/
/*#endif*/

void ecies_init(){
    ecc_init();
  }

void KDF(uint8_t *K, int K_len, uint8_t *Z){
    int len, i;
    uint8_t z[KEYDIGITS*NN_DIGIT_LEN+4];
    SHA1Context ctx;
    uint8_t sha1sum[20];

    memcpy(z, Z, KEYDIGITS*NN_DIGIT_LEN);
    memset(z + KEYDIGITS*NN_DIGIT_LEN, 0, 3);
    //KDF
    //|z|+|ShareInfo|+4 < 2^64, no need to check
    //keydatalen < 20*(2^32-1), no need to check
    len = K_len;
    i = 1;
    while(len > 0){
      z[KEYDIGITS*NN_DIGIT_LEN + 3] = i;
      sha1_reset(&ctx);
      sha1_update(&ctx, z, KEYDIGITS*NN_DIGIT_LEN+4);
      sha1_digest(&ctx, sha1sum);
      if(len >= 20)
      {
				memcpy(K+(i-1)*20, sha1sum, 20);
      }else
      {
				memcpy(K+(i-1)*20, sha1sum, len);
      }
      i++;
      len = len - 20;
    }
  }
  
  //refer to RFC 2104
  void hmac_sha1(uint8_t *text, int text_len, uint8_t *key, int key_len, uint8_t *digest){
    SHA1Context context;
    uint8_t k_ipad[65];    /* inner padding - * key XORd with ipad		    */
    uint8_t k_opad[65];    /* outer padding - * key XORd with opad		    */
    uint8_t tk[20];
    int i;
    /* if key is longer than 64 bytes reset it to key=MD5(key) */
    if (key_len > 64) 
    {  
      SHA1Context      tctx;
      
      sha1_reset(&tctx);
      sha1_update(&tctx, key, key_len);
      sha1_digest(&tctx, tk);
      
      key = tk;
      key_len = 20;
    }
    
    /*
     * the HMAC_SHA1 transform looks like:
     *
     * SHA1(K XOR opad, SHA1(K XOR ipad, text))
     *
     * where K is an n byte key
     * ipad is the byte 0x36 repeated 64 times
     
     * opad is the byte 0x5c repeated 64 times
     * and text is the data being protected
     */
    
    /* start out by storing key in pads */
    memcpy(k_ipad, key, key_len);
		memset(k_ipad + key_len, 0, 65 - key_len);
    memcpy(k_opad, key, key_len);
    memset(k_opad + key_len, 0, 65 - key_len);
    
    /* XOR key with ipad and opad values */
    for (i=0; i<64; i++) {
      k_ipad[i] ^= 0x36;
      k_opad[i] ^= 0x5c;
    }
    /*
     * perform inner SHA1
     */
    sha1_reset(&context);                   /* init context for 1st pass */
    sha1_update(&context, k_ipad, 64);      /* start with inner pad */
    sha1_update(&context, text, text_len); /* then text of datagram */
    sha1_digest(&context, digest);          /* finish up 1st pass */
    /*
     * perform outer SHA1
     */
    sha1_reset(&context);                   /* init context for 2nd pass */
    sha1_update(&context, k_opad, 64);     /* start with outer pad */
    sha1_update(&context, digest, 20);
    sha1_digest(&context, digest);         /* then results of 1st hash */

  }
  
void gen_random(NN_DIGIT *a, uint8_t length)
{
/*
  a[5] = 0x00000000;
  a[4] = 0x7b012db7;
  a[3] = 0x681a3f28;
  a[2] = 0xb9185c8b;
  a[1] = 0x2ac5d528;
  a[0] = 0xdecd52da;
  uint8_t ri;
*/	
/*  int ri;*/
/*  for(ri=0; ri<length; ri++) { */
/*    srand(100);*/
/*    a[ri] = ((uint32_t)rand() << 16)^((uint32_t)rand());*/
/*  }*/
} 

int encrypt(uint8_t *C, int C_len, uint8_t *M, int M_len, Point *PublicKey){
    NN_DIGIT k[NUMWORDS];
    uint8_t z[KEYDIGITS*NN_DIGIT_LEN];
    Point R, P;
    //uint8_t octet_buf[2*KEYDIGITS*NN_DIGIT_LEN];
    int octet_len;
    uint8_t K[MAX_M_LEN + HMAC_LEN];
    int i;

    if(C_len < KEYDIGITS*NN_DIGIT_LEN+1+HMAC_LEN)
      return -1;

    if(C_len < 2*KEYDIGITS*NN_DIGIT_LEN+1+HMAC_LEN)
      return -1;
      
    //1. select key pair
    
		//gen_random(k, NUMWORDS);
/*		    k[5] = 0x0;*/
/*		    k[4] = 0x7b012db7;*/
/*		    k[3] = 0x681a3f28;*/
/*		    k[2] = 0xb9185c8b;*/
/*		    k[1] = 0x2ac5d528;*/
/*		    k[0] = 0xdecd52da;*/
		   
/*		   	k[5] = 0x00000000;*/
/*				k[4] = 0xc36e3e96;*/
/*				k[3] = 0xc26c3d91;*/
/*				k[2] = 0xc7ec7db1;*/
/*				k[1] = 0xd7e47933;*/
/*				k[0] = 0x16020c0d;*/

 		//random
		ecc_gen_private_key(k);
   
    ecc_gen_public_key(&R, k);

    //2. convert R to octet string
#ifdef POINT_COMPRESS
    octet_len = ecc_point2octet(C, C_len, &R, TRUE);
#else  //no point compression
    octet_len = ecc_point2octet(C, C_len, &R, FALSE);
#endif
    //3. derive shared secret z=P.x
#ifdef SLIDING_WIN
#ifdef PROJECTIVE
    ecc_win_precompute_Z(PublicKey, baseArray, ZList);
    ecc_win_mul_Z(&P, k, baseArray, ZList);
#else
    ecc_win_precompute(PublicKey, baseArray);
    ecc_win_mul(&P, k, baseArray);
#endif //PROJECTIVE
#else  //SLIDING_WIN
    ecc_mul(&P, PublicKey, k);
#endif  //SLIDING_WIN

    if (p_iszero(&P))
      return -1;

    //4. convert z to octet string Z
    NN_Encode(z, KEYDIGITS*NN_DIGIT_LEN, P.x, NUMWORDS);

    //5. use KDF to generate K of length enckeylen + mackeylen octets from Z
    //enckeylen = M_len, mackeylen = 20
    KDF(K, M_len+HMAC_LEN, z);

    //6. the left most enckeylen octets of K is EK, right most mackeylen octets is MK

    //7. encrypt EM
    for (i=0; i<M_len; i++){
      C[octet_len+i] = M[i] ^ K[i];
    }

    //8. generate mac D
    hmac_sha1(C + octet_len, M_len, K + M_len, HMAC_LEN, C + octet_len + M_len);

    //9. output C = R||EM||D
    return (octet_len + M_len + HMAC_LEN);    
  }
  
int decrypt(uint8_t *M, int M_len, uint8_t *C, int C_len, NN_DIGIT *d){

    uint8_t z[KEYDIGITS*NN_DIGIT_LEN];
    Point R,P; 
    int octet_len;
    uint8_t K[MAX_M_LEN + HMAC_LEN];
    int i;    
    uint8_t hmac_tmp[HMAC_LEN];

    //1. parse R||EM||D
    
    //2. get the point R
    octet_len = ecc_octet2point(&R, C, C_len);

    //3. make sure R is valid
    if (ecc_check_point(&R) != 1)
      return -1;
 
    //4. use private key to generate shared secret z
#ifdef SLIDING_WIN
#ifdef PROJECTIVE
    ecc_win_precompute_Z(&R, baseArray, ZList);
    ecc_win_mul_Z(&P, d, baseArray, ZList);
#else
    ecc_win_precompute(&R, baseArray);
    ecc_win_mul(&P, d, baseArray);
#endif //PROJECTIVE
#else  //SLIDING_WIN
    ecc_mul(&P, &R, d);
#endif  //SLIDING_WIN

    if (p_iszero(&P))
      return -1;    

    //5. convert z to octet string Z
    NN_Encode(z, KEYDIGITS*NN_DIGIT_LEN, P.x, NUMWORDS);

    //6. use KDF to derive EK and MK
    KDF(K, C_len - octet_len, z);

    //7. check D first
    if (M_len < C_len - HMAC_LEN - octet_len)
      return -1;
      
    M_len = C_len - HMAC_LEN - octet_len;
    hmac_sha1(C + octet_len, M_len, K + M_len, HMAC_LEN, hmac_tmp);

    for (i=0; i<HMAC_LEN; i++){
      if (hmac_tmp[i] != C[octet_len + M_len + i])
				return -2;
    }
    
    //8. decrypt
    for(i=0; i<M_len; i++){
      M[i] = C[octet_len+i] ^ K[i];
    }
    
    return M_len;
  }
