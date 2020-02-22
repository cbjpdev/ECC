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

NN_DIGIT prKey_alice[NUMWORDS];
Point pbkey_alice;
static Params param;
struct fulMsg *fmsg = &fulMsg_t;

static void
recv_uc(struct unicast_conn *c, const rimeaddr_t *from)
{ 
  //printf("unicast message received from %d.%d: %s\n",from->u8[0], from->u8[1], (uint8_t *)packetbuf_dataptr());
  printf("unicast message received from %d.%d: \n",from->u8[0], from->u8[1]);
  
  fmsg = (fulMsg_tt *)(packetbuf_dataptr());

  printf("msg int 1:%d \n",fmsg->msgType);
  printf("msg int 2:%d \n",fmsg->whoEncrypts);
	printf("msg int 3:%d \n",fmsg->whoDecrypts);
	printf("msg int 4:%d \n",fmsg->subGroupNodeOne);
	printf("msg int 5:%d \n",fmsg->subGroupNodeTwo);
	printf("msg int 6:%d \n",fmsg->subGroupNodeThree);
	printf("msg int 7:%d \n",fmsg->subGroupNodeFour);
	
	printf("id 2 Recieved Cipher text: ");
  int ii = 0;
	for (ii = 0; ii != 20; ++ii) 
	{
		if (ii%32 == 0 && ii != 0) printf("\n");
		printf("%02x", fmsg->C[ii]);
	}
	printf("\n");
	
  
  uint8_t *dM = malloc(MAX_M_LEN*sizeof(uint8_t));
  int dM_len = MAX_M_LEN;
  
  uint8_t *C = malloc((2*KEYDIGITS*NN_DIGIT_LEN + 1 + MAX_M_LEN + HMAC_LEN)*sizeof(uint8_t));	
  int C_len =102;
  
  C = fmsg->C;
  
  dM_len = decrypt(dM, dM_len, C, C_len, fmsg->sintops); 

/*----------------------print texts----------------------------- */ 
  printf("id 2 dm text: ");
  int i = 0;
	for (i = 0; i != 20; ++i) 
	{
		if (i%32 == 0 && i != 0) printf("\n");
		printf("%02x", dM[i]);
	}
	printf("\n");
/*----------------------print texts ends------------------------  */
	
}
/*---------------------------------------------------------------------------*/

static const struct unicast_callbacks unicast_callbacks = {recv_uc};

static struct unicast_conn uc;

/*---------------------------------------------------------------------------*/
static void unicast_message()
{
	
  rimeaddr_t addr;
  //packetbuf_copyfrom(dM, MAX_M_LEN);
  addr.u8[0] = 1;
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
  unicast_open(&uc, 146, &unicast_callbacks);
  
  button_sensor.configure(SENSORS_ACTIVE, 1); 
  process_start(&unicast_process, NULL);

  PROCESS_END();
}
