#ifndef _ECIES_H_
#define _ECIES_H_

#define MAX_M_LEN 41
#define HMAC_LEN 20

#include "nn.h"

void ecies_init();
void KDF(uint8_t *K, int K_len, uint8_t *Z);
void hmac_sha1(uint8_t *text, int text_len, uint8_t *key, int key_len, uint8_t *digest);
void gen_random(NN_DIGIT *a, uint8_t length);
int encrypt(uint8_t *C, int C_len, uint8_t *M, int M_len, Point *PublicKey);
int decrypt(uint8_t *M, int M_len, uint8_t *C, int C_len, NN_DIGIT *d);

#endif
