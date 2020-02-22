#ifndef __EDSA_H__
#define __EDSA_H__

#include "nn.h"

/*
 * Initialize the ECDSA, pKey is public key used to verify the signature.
 * we assume that the node has already know the public key when node is deployed 	
 */
int ecdsa_init(Point * pKey);

/*
 * 
 */
void ecdsa_sign(uint8_t *msg, uint16_t len, NN_DIGIT *r, NN_DIGIT *s, NN_DIGIT *d);

uint8_t ecdsa_verify(uint8_t *msg, uint16_t len, NN_DIGIT *r, NN_DIGIT *s, Point *Q);


#endif /* __EDSA_H__ */
