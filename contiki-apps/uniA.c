#include "ecc.h"
#include "ecies.h"
#include "contiki.h"
#include "curve_param.h"
#include "lib/random.h"
#include "net/rime.h"
#include "dev/button-sensor.h"

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
PROCESS(unicast_process, "Example unicast");
PROCESS(startup_process, "Statup Process");
AUTOSTART_PROCESSES(&startup_process);
/*---------------------------------------------------------------------------*/

Point pbkey_alice;
NN_DIGIT prKey_alice[NUMWORDS];
static Params param;
struct fulMsg *fmsg = &fulMsg_t;

static void 
recv_uc(struct unicast_conn *c, const rimeaddr_t *from)
{
  //printf("broadcast message received from %d.%d: '%s'\n",from->u8[0], from->u8[1], (char *)packetbuf_dataptr());
  printf("unicast message received from %d.%d: \n",from->u8[0], from->u8[1]);
  //uint8_t *rM= malloc(MAX_M_LEN*sizeof(uint8_t));
  //rM = (uint8_t *)(packetbuf_dataptr());  
  //if(*M == *rM) printf("Works\n");
}

static const struct unicast_callbacks unicast_callbacks = {recv_uc};

static struct unicast_conn uc;

/*---------------------------------------------------------------------------*/
static void random_data(void *ptr, uint16_t len)
{
  uint16_t i;
  for(i=0; i<len; i++) {
    srand(100);
   ((uint8_t *)(ptr))[i] = rand() % 100; 
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static void unicast_message()
{

	uint8_t *M = malloc(MAX_M_LEN*sizeof(uint8_t));
  int M_len = MAX_M_LEN;
  
  uint8_t *C = malloc((2*KEYDIGITS*NN_DIGIT_LEN + 1 + MAX_M_LEN + HMAC_LEN)*sizeof(uint8_t));	
  int C_len;
  
  random_data(M, MAX_M_LEN);
  
  C_len = encrypt(C, (2*KEYDIGITS*NN_DIGIT_LEN + 1 + M_len + HMAC_LEN), M, M_len, &pbkey_alice);
  printf("C len: %d\n",C_len);
  
/*-----------------------------------------------print texts-----------------------------  */
  printf("id1:Plain text: ");
  int i = 0;
	for (i = 0; i != 20; ++i) 
	{
		if (i%32 == 0 && i != 0) printf("\n");
		printf("%02x", M[i]);
	}
	printf("\n");
/*-----------------------------------------------print texts ends------------------------*/  

	fmsg->msgType = 1;
	fmsg->C = C;
	fmsg->whoEncrypts = 2;
	fmsg->whoDecrypts = 3;
	fmsg->subGroupNodeOne = 4;
	fmsg->subGroupNodeTwo = 5;
	fmsg->subGroupNodeThree = 6;
	fmsg->subGroupNodeFour= 7; 
	int a = 0;
	for(a = 0; a< 6; a++)
	fmsg->sintops[a] = prKey_alice[a];
	
	printf("id1:Cipher text: ");
	int ii = 0;
	for (ii = 0; ii != 20; ++ii) 
	{
		if (ii%32 == 0 && ii != 0) printf("\n");
		printf("%02x", fmsg->C[ii]);
	}
	printf("\n");
	
  packetbuf_clear();

	rimeaddr_t addr;
  packetbuf_copyfrom(fmsg, 110);
  addr.u8[0] = 2;
  addr.u8[1] = 0;
  if(!rimeaddr_cmp(&addr, &rimeaddr_node_addr))
  {
  	unicast_send(&uc, &addr);
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
PROCESS_THREAD(unicast_process, ev, data)
{
  PROCESS_EXITHANDLER(unicast_close(&uc);)   
  PROCESS_BEGIN();  
  while(1) {
    PROCESS_WAIT_EVENT_UNTIL(ev == sensors_event && data == &button_sensor);
    unicast_message();
    printf("message sent.\n");
  }     
  PROCESS_END();
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
PROCESS_THREAD(startup_process, ev, data)
{
  PROCESS_BEGIN();

  memset(prKey_alice, 0, NUMWORDS*NN_DIGIT_LEN);
  memset(pbkey_alice.x, 0, NUMWORDS*NN_DIGIT_LEN);
  memset(pbkey_alice.y, 0, NUMWORDS*NN_DIGIT_LEN);

  /* set public key for Alice */
  pbkey_alice.x[5] = 0x00000000;
  pbkey_alice.x[4] = 0x21961f69;
  pbkey_alice.x[3] = 0xf02d202b;
  pbkey_alice.x[2] = 0xa4b41f1a;
  pbkey_alice.x[1] = 0x0aa08a86;
  pbkey_alice.x[0] = 0xdf27908d;
    
  pbkey_alice.y[5] = 0x00000000;
  pbkey_alice.y[4] = 0x378e1278;
  pbkey_alice.y[3] = 0x62836d75;
  pbkey_alice.y[2] = 0x7acb7ca4;
  pbkey_alice.y[1] = 0x0dc0ad13;
  pbkey_alice.y[0] = 0x741e287c;

  prKey_alice[5] = 0x00000000;
  prKey_alice[4] = 0xc36e3e96;
  prKey_alice[3] = 0xc26c3d91;
  prKey_alice[2] = 0xc7ec7db1;
  prKey_alice[1] = 0xd7e47933;
  prKey_alice[0] = 0x16020c0d;
  
/*  fmsg->sintops[5] = 0x00000000;*/
/*  fmsg->sintops[4] = 0xc36e3e96;*/
/*  fmsg->sintops[3] = 0xc26c3d91;*/
/*  fmsg->sintops[2] = 0xc7ec7db1;*/
/*  fmsg->sintops[1] = 0xd7e47933;*/
/*  fmsg->sintops[0] = 0x16020c0d;*/

  /* Initialize ecc. */
  ecc_init();
  ecies_init();
  get_curve_param(&param);  
	unicast_open(&uc, 146, &unicast_callbacks);  
  button_sensor.configure(SENSORS_ACTIVE, 1); 
  process_start(&unicast_process, NULL);

  PROCESS_END();
}
