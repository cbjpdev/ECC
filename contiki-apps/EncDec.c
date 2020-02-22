#include "ecc.h"
#include "ecies.h"
#include "ecdsa.h"
#include "contiki.h"
#include "curve_param.h"
#include "lib/random.h"
#include "net/rime.h"
#include "dev/button-sensor.h"
#include "dev/leds.h"

#include "messages.h"

#include <stdio.h> /* For printf() */
#include <string.h>
#include <stdlib.h>

#define MAX_M_LEN 41
#define HMAC_LEN 20
#define MAX_ROUNDS 10
#define MSG_LEN 20
#define MAX_LEN 20

/*---------------------------------------------------------------------------*/
PROCESS(bob_process, "Alice process");
PROCESS(startup_process, "Statup Process");
AUTOSTART_PROCESSES(&startup_process);
/*---------------------------------------------------------------------------*/

static Params param;
Point pbkey_alice;
Point paraInverse;
NN_DIGIT prKey_alice[NUMWORDS];

NN_DIGIT prKey_alice1[NUMWORDS];
NN_DIGIT prKey_alice2[NUMWORDS];
/*NN_DIGIT prKey_alice3[NUMWORDS];*/
/*NN_DIGIT prKey_alice4[NUMWORDS];*/
/*NN_DIGIT prKey_alice5[NUMWORDS];*/

/*---------------------------------------------------------------------------*/
PROCESS_THREAD(startup_process, ev, data)
{
  PROCESS_BEGIN();

  memset(prKey_alice, 0, NUMWORDS*NN_DIGIT_LEN);
  memset(pbkey_alice.x, 0, NUMWORDS*NN_DIGIT_LEN);
  memset(pbkey_alice.y, 0, NUMWORDS*NN_DIGIT_LEN);

  /* set public key for Alice */
 
/*  pbkey_alice.x[5] = 0x00000000;*/
/*  pbkey_alice.x[4] = 0x21961f69;*/
/*  pbkey_alice.x[3] = 0xf02d202b;*/
/*  pbkey_alice.x[2] = 0xa4b41f1a;*/
/*  pbkey_alice.x[1] = 0x0aa08a86;*/
/*  pbkey_alice.x[0] = 0xdf27908d;*/
/*    */
/*  pbkey_alice.y[5] = 0x00000000;*/
/*  pbkey_alice.y[4] = 0x378e1278;*/
/*  pbkey_alice.y[3] = 0x62836d75;*/
/*  pbkey_alice.y[2] = 0x7acb7ca4;*/
/*  pbkey_alice.y[1] = 0x0dc0ad13;*/
/*  pbkey_alice.y[0] = 0x741e287c;*/

/*  prKey_alice[5] = 0x00000000;*/
/*  prKey_alice[4] = 0xc36e3e96;*/
/*  prKey_alice[3] = 0xc26c3d91;*/
/*  prKey_alice[2] = 0xc7ec7db1;*/
/*  prKey_alice[1] = 0xd7e47933;*/
/*  prKey_alice[0] = 0x16020c0d;*/

  /* Initialize ecc. */
  ecc_init();

  ecies_init();
  
  get_curve_param(&param);  
  
  //NN_ModInv (a, b, c, digits)     Computes a = 1/b mod c
	//NN_Mult (a, b, c, digits)       Computes a = b * c.
	//NN_Add (a, b, c, digits) Computes a = b + c.
	 
  //NN_Add(prKey_alice4, prKey_alice1, prKey_alice2, NUMWORDS);
	//NN_Add(prKey_alice5, prKey_alice4, prKey_alice3, NUMWORDS);
  
  NN_Mult(prKey_alice2, prKey_alice1, param.G, NUMWORDS);
  NN_ModInv(paraInverse, param.G, param.p, NUMWORDS);
  NN_Mult(prKey_alice, paraInverse, prKey_alice2, NUMWORDS);
  
  ecc_gen_public_key(&pbkey_alice, prKey_alice);
  
  //ecc_gen_private_key(prKey_alice);
  //ecc_gen_public_key(&pbkey_alice, prKey_alice);
  
  button_sensor.configure(SENSORS_ACTIVE, 1); 
  process_start(&bob_process, NULL);

  PROCESS_END();
}
/*---------------------------------------------------------------------------*/
static void
random_data(void *ptr, uint16_t len)
{
  uint16_t i;
  for(i=0; i<len; i++) {
    srand(100);
    ((uint8_t *)(ptr))[i] = rand() % 100; 
  }
}
/*---------------------------------------------------------------------------*/
static void enc_dec_message()
{
  uint8_t *M = malloc(MAX_M_LEN*sizeof(uint8_t));
  int M_len = MAX_M_LEN;
  
  uint8_t *C = malloc((2*KEYDIGITS*NN_DIGIT_LEN + 1 + MAX_M_LEN + HMAC_LEN)*sizeof(uint8_t));	
  int C_len; 
  
  uint8_t *dM = malloc(MAX_M_LEN*sizeof(uint8_t));
  int dM_len = MAX_M_LEN;
  
  random_data(M, MAX_M_LEN);

	C_len = encrypt(C, (2*KEYDIGITS*NN_DIGIT_LEN + 1 + M_len + HMAC_LEN), M, M_len, &pbkey_alice);
	
	dM_len = decrypt(dM, dM_len, C, C_len, prKey_alice); 
	
	/*-----------------------------------------------print texts-----------------------------*/  
  printf("id1:Plain text: ");
  int i = 0;
	for (i = 0; i != 20; ++i) 
	{
		if (i%32 == 0 && i != 0) printf("\n");
		printf("%02x", M[i]);
	}
	printf("\n");
/*-----------------------------------------------print texts ends------------------------*/  
	
/*-----------------------------------------------print texts-----------------------------*/  	
	printf("id2:Cipher text: ");
	int ii = 0;
	for (ii = 0; ii != 20; ++ii) 
	{
		if (ii%32 == 0 && ii != 0) printf("\n");
		printf("%02x", C[ii]);
	}
	printf("\n");
/*-----------------------------------------------print texts ends------------------------*/
/*-----------------------------------------------print texts-----------------------------*/  	
	printf("id3:plain text: ");
	int iii = 0;
	for (iii = 0; iii != 20; ++iii) 
	{
		if (iii%32 == 0 && iii != 0) printf("\n");
		printf("%02x", dM[iii]);
	}
	printf("\n");
/*-----------------------------------------------print texts ends------------------------*/
		
	if (*dM == *M) printf("works\n");
}

PROCESS_THREAD(bob_process, ev, data)
{
  PROCESS_BEGIN();
  while(1) {
    PROCESS_WAIT_EVENT_UNTIL(ev == sensors_event && data == &button_sensor);
    enc_dec_message();
    printf("message sent.\n");
  }  
  PROCESS_END();
}
/*---------------------------------------------------------------------------*/
