#include "ecc.h"
#include "ecies.h"
#include "contiki.h"
#include "curve_param.h"
#include "lib/random.h"
#include "net/rime.h"
#include "dev/button-sensor.h"
#include "node-id.h"
#include "messages.h"
#include "timer.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAX_M_LEN 41
#define HMAC_LEN 20
/*---------------------------------------------------------------------------*/
PROCESS(startup_process, "Statup Process");
PROCESS(unicast_key_eNode, "unicast key eNode");//msgType 3 
PROCESS(unicast_key_dNode, "unicast key dNode");//msgType 4
AUTOSTART_PROCESSES(&startup_process);
/*---------------------------------------------------------------------------*/

static void node_unicast_leder_encrypt_request(); //msgType 1 - encrypt request from a node - unicast to leader
static void leader_unicast_node_encypt_request(int whoEncrypts); //msgType 2 - 
static void encryptAndSend(); //msgType 4

static void node_unicast_leder_decrypt_request(); //msgType 5
static void leader_unicast_node_decypt_request(); //msgType 6
static void decryptMessage();

NN_DIGIT prKey_self[NUMWORDS];

NN_DIGIT prKey_selfM[NUMWORDS];

NN_DIGIT prKey_Eshare4[NUMWORDS];
NN_DIGIT prKey_Eshare5[NUMWORDS];
NN_DIGIT prKey_Eshare6[NUMWORDS];
NN_DIGIT prKey_Eshare7[NUMWORDS];

NN_DIGIT prKey_Dshare4[NUMWORDS];
NN_DIGIT prKey_Dshare5[NUMWORDS];
NN_DIGIT prKey_Dshare6[NUMWORDS];
NN_DIGIT prKey_Dshare7[NUMWORDS];

NN_DIGIT prKey_group[NUMWORDS];

int flag5e = 0,flag6e = 0,flag7e = 0,flag4d = 0,flag5d = 0,flag6d = 0;

NN_DIGIT s_s[NUMWORDS];
Point pbkey_group;
static Params param;
int WhoIsLeader;
//uint8_t *rC = malloc((2*KEYDIGITS*NN_DIGIT_LEN + 1 + MAX_M_LEN + HMAC_LEN)*sizeof(uint8_t));	
uint8_t *rC ;
static void random_data(void *ptr, uint16_t len)
{
  uint16_t i;
  for(i=0; i<len; i++) {
    srand(100);
   ((uint8_t *)(ptr))[i] = rand() % 100; 
  }
}

static void 
recv_uc(struct unicast_conn *c, const rimeaddr_t *from)
{
  struct fulMsg *fmsg = &fulMsg_t;
  fmsg = (fulMsg_tt *)(packetbuf_dataptr());  
  
  if (fmsg->msgType == 1)
  {
  	printf("I got a message from %d leader_unicast_node_encypt_request %d \n ",from->u8[0],fmsg->msgType);
  	leader_unicast_node_encypt_request(fmsg->whoEncrypts); // E - method 2
  }
  
  if (fmsg->msgType == 2)
  {
   	printf("I got a message from %d , leader_unicast_node_encypt_request %d \n ",from->u8[0],fmsg->msgType);
  	process_start(&unicast_key_eNode, NULL);
  }
  
  if (fmsg->msgType == 3)
  {
			if(from->u8[0] == 5)
			{
				printf("I got a message from %d , I got a key share from %d \n ",from->u8[0],fmsg->msgType);
				int a = 0;
				for(a = 0; a< 6; a++)
				prKey_Eshare5[a] = fmsg->sintops[a];			
				
				flag5e = 1;
				if ((flag5e == 1) & (flag6e == 1) & (flag7e == 1))
				{
					encryptAndSend();
				}
			}
		
			if(from->u8[0] == 6)
			{
				printf("I got a message from %d , I got a key share from %d \n ",from->u8[0],fmsg->msgType);
				int a = 0;
				for(a = 0; a< 6; a++)
				prKey_Eshare6[a] = fmsg->sintops[a];
				flag6e = 1;
				if ((flag5e == 1) & (flag6e == 1) & (flag7e == 1))
				{
					encryptAndSend();
				}
			}
			if(from->u8[0] == 7)
			{
				printf("I got a message from %d , I got a key share from %d \n ",from->u8[0],fmsg->msgType);
				int a = 0;
				for(a = 0; a< 6; a++)
				prKey_Eshare7[a] = fmsg->sintops[a];
				flag7e = 1;
				if ((flag5e == 1) & (flag6e == 1) & (flag7e == 1))
				{
					encryptAndSend();
				}
			}
  }
  
  if (fmsg->msgType == 4)
  {
  	printf("I got a message from %d , cipher text %d \n ",from->u8[0],fmsg->msgType); 
  	rC = fmsg->C;
  }
  
  if (fmsg->msgType == 5)
  {
  	if (from->u8[0] == 2)
  	{
  		printf("bad request\n");
  	}
  	else
  	{
  		printf("I got a message from %d , decrypt_request %d \n ",from->u8[0],fmsg->msgType);
  		leader_unicast_node_decypt_request(fmsg->whoDecrypts);
  	}
 	}
  if (fmsg->msgType == 6)
  {
  	printf("I got a message from %d , leader_unicast_node_decrypt_request %d \n ",from->u8[0],fmsg->msgType);
  	process_start(&unicast_key_dNode, NULL);
  }
  if (fmsg->msgType == 7)
  {
		if(from->u8[0] == 4)
				{
					printf("I got a message from %d , I got a key share from %d \n ",from->u8[0],fmsg->msgType);
					int a = 0;
					for(a = 0; a< 6; a++)
					prKey_Dshare4[a] = fmsg->sintops[a];
					flag4d = 1;
					if ((flag4d == 1) & (flag5d == 1) & (flag6d == 1))
					{
						decryptMessage();
					}
				}
		
				if(from->u8[0] == 5)
				{
					printf("I got a message from %d , I got a key share from %d \n ",from->u8[0],fmsg->msgType);
					int a = 0;
					for(a = 0; a< 6; a++)
					prKey_Dshare5[a] = fmsg->sintops[a];
					flag5d = 1;
					if ((flag4d == 1) & (flag5d == 1) & (flag6d == 1))
					{
						decryptMessage();
					}
				}
				if(from->u8[0] == 6)
				{
					printf("I got a message from %d , I got a key share from %d \n ",from->u8[0],fmsg->msgType);
					int a = 0;
					for(a = 0; a< 6; a++)
					prKey_Dshare6[a] = fmsg->sintops[a];
					flag6d = 1;
					if ((flag4d == 1) & (flag5d == 1) & (flag6d == 1))
					{
						decryptMessage();
					}
					//else printf("Someone miss in the group. Please re initiate the process\n");
				}			
	}
}

static const struct unicast_callbacks unicast_callbacks = {recv_uc};

static struct unicast_conn uc;

static void node_unicast_leder_encrypt_request() //ask leader to recive the keys to encrypt  //E - method one 
{
	WhoIsLeader = 10;
	struct fulMsg *fmsgOne = &fulMsg_t;
	fmsgOne->msgType = 1;//encrypt request
	fmsgOne->whoEncrypts = 4;
	fmsgOne->whoDecrypts = 0;
	fmsgOne->subGroupNodeOne = 4;
	fmsgOne->subGroupNodeTwo = 5;
	fmsgOne->subGroupNodeThree = 6;
	fmsgOne->subGroupNodeFour= 7;
	
	//ecc_gen_private_key(prKey_self);
	//NN_Mult(prKey_selfM, prKey_self, s_s, NUMWORDS);
	
  packetbuf_clear();

  rimeaddr_t addr;
	packetbuf_copyfrom(fmsgOne, 110);
  addr.u8[0] = WhoIsLeader;
  addr.u8[1] = 0;
  if(!rimeaddr_cmp(&addr, &rimeaddr_node_addr))
  {
  	unicast_send(&uc, &addr);
  }
}

static void leader_unicast_node_encypt_request(int WhichNodeWantEncrypt) //E - method 2
{	
	struct fulMsg *fmsgLtNsE = &fulMsg_t;
	fmsgLtNsE->msgType = 2;//request
	fmsgLtNsE->whoEncrypts = WhichNodeWantEncrypt;
	fmsgLtNsE->whoDecrypts = 0;
	fmsgLtNsE->subGroupNodeOne = 4;
	fmsgLtNsE->subGroupNodeTwo = 5;
	fmsgLtNsE->subGroupNodeThree = 6;
	fmsgLtNsE->subGroupNodeFour= 7;
	
	int i = 0; 
	for(i = 5; i < 8; i++)
	{
		rimeaddr_t addr;
		packetbuf_clear();
		packetbuf_copyfrom(fmsgLtNsE, 110);
		addr.u8[0] = i;
		addr.u8[1] = 0;
		if(!rimeaddr_cmp(&addr, &rimeaddr_node_addr))
		{
			unicast_send(&uc, &addr);
		}	
	}
}

PROCESS_THREAD(unicast_key_eNode, ev, data)
{
	PROCESS_EXITHANDLER(unicast_close(&uc);)
	PROCESS_BEGIN();
		unicast_open(&uc, 146, &unicast_callbacks);
	  static struct etimer et;
		rimeaddr_t addr;
		etimer_set(&et, CLOCK_SECOND);	
		PROCESS_WAIT_EVENT_UNTIL(etimer_expired(&et));		
		//generate key part
		//unicast them to encryption nodeOne
		//NN_Mult (a, b, c, digits)       Computes a = b * c.
/*		ecc_gen_private_key(prKey_self);*/
/*		NN_Mult(prKey_selfM, prKey_self, s_s, NUMWORDS);*/
		//send this to encryption node

		//create the struct and then send

		struct fulMsg *fmsgThree = &fulMsg_t;
		fmsgThree->msgType = 3;//request
		fmsgThree->whoEncrypts = 4;
		fmsgThree->whoDecrypts = 0;
		fmsgThree->subGroupNodeOne = 4;
		fmsgThree->subGroupNodeTwo = 5;
		fmsgThree->subGroupNodeThree = 6;
		fmsgThree->subGroupNodeFour= 7;
		int a = 0;
		for(a = 0; a< 6; a++)
		fmsgThree->sintops[a] = prKey_selfM[a];
	
		packetbuf_clear();

		packetbuf_copyfrom(fmsgThree, 110);
		addr.u8[0] = 4;
		addr.u8[1] = 0;
		if(!rimeaddr_cmp(&addr, &rimeaddr_node_addr))
		{			
			unicast_send(&uc, &addr);
		}
	PROCESS_END();
}

static void encryptAndSend()
{
	//NN_Add (a, b, c, digits) Computes a = b + c.
	//NN_ModInv (a, b, c, digits)     Computes a = 1/b mod c
	//NN_Mult (a, b, c, digits)       Computes a = b * c.
	
	NN_DIGIT Etemp1[NUMWORDS];
	NN_DIGIT Etemp2[NUMWORDS];
	NN_DIGIT Etemp3[NUMWORDS];
	NN_DIGIT Etemp4[NUMWORDS];

	memset(Etemp1, 0, NUMWORDS*NN_DIGIT_LEN);
	memset(Etemp2, 0, NUMWORDS*NN_DIGIT_LEN);
	memset(Etemp3, 0, NUMWORDS*NN_DIGIT_LEN);
	memset(Etemp4, 0, NUMWORDS*NN_DIGIT_LEN);
	
	int a = 0;
	for(a = 0; a< 6; a++)
	prKey_Eshare4[a] = prKey_selfM[a];
	

	
	NN_Add (Etemp1, prKey_Eshare4, prKey_Eshare5, NUMWORDS);
	NN_Add (Etemp2, prKey_Eshare6, prKey_Eshare7, NUMWORDS);
	NN_Add (Etemp3, Etemp2, Etemp1, NUMWORDS);
	
	
/*	int b1 = 0;*/
/*  for(b1 = 0; b1 <6; b1++ )*/
/*  {*/
/*  	printf("on mote 4 pr4 %p\n",Etemp3[b1]);  */
/*  }*/

	NN_ModInv(Etemp4, s_s, param.p, NUMWORDS);
  NN_Mult(prKey_group, Etemp3, Etemp4, NUMWORDS);
  
	ecc_gen_public_key(&pbkey_group, prKey_group);
	
	uint8_t *M = malloc(MAX_M_LEN*sizeof(uint8_t));
  int M_len = MAX_M_LEN;
  
  uint8_t *C = malloc((2*KEYDIGITS*NN_DIGIT_LEN + 1 + MAX_M_LEN + HMAC_LEN)*sizeof(uint8_t));	
  int C_len = 102;
  
  random_data(M, MAX_M_LEN);
  
  C_len = encrypt(C, (2*KEYDIGITS*NN_DIGIT_LEN + 1 + M_len + HMAC_LEN), M, M_len, &pbkey_group);
  
  struct fulMsg *fmsgFour = &fulMsg_t;
  fmsgFour->msgType = 4;
	fmsgFour->C = C;
	fmsgFour->whoEncrypts = 2;
	fmsgFour->whoDecrypts = 3;
	fmsgFour->subGroupNodeOne = 4;
	fmsgFour->subGroupNodeTwo = 5;
	fmsgFour->subGroupNodeThree = 6;
	fmsgFour->subGroupNodeFour= 7; 	
  packetbuf_clear();
  
  printf("Temperature : ");
  int i = 0;
	for (i = 0; i != 1; ++i) 
	{
		if (i%32 == 0 && i != 0) printf("\n");
		printf("%02x 'C ", M[i]);
	}
	printf("\n");
	
	packetbuf_clear();
		
	unicast_open(&uc, 146, &unicast_callbacks);

	rimeaddr_t addr;
	packetbuf_copyfrom(fmsgFour, 110);
	addr.u8[0] = 7;
	addr.u8[1] = 0;
	if(!rimeaddr_cmp(&addr, &rimeaddr_node_addr))
	{
		unicast_send(&uc, &addr);
	}
	
	packetbuf_clear();
		
	unicast_open(&uc, 146, &unicast_callbacks);

	packetbuf_copyfrom(fmsgFour, 110);
	addr.u8[0] = 2;
	addr.u8[1] = 0;
	if(!rimeaddr_cmp(&addr, &rimeaddr_node_addr))
	{
		unicast_send(&uc, &addr);
	}		
	memset(prKey_group, 0, NUMWORDS*NN_DIGIT_LEN);
}
//------------------------------------------------------------------------------------------------------------------//
static void node_unicast_leder_decrypt_request() //ask leader to recive the keys to encrypt
{
	WhoIsLeader = 10;
	struct fulMsg *fmsgFive = &fulMsg_t;
	fmsgFive->msgType = 5;//decrypt request
	fmsgFive->whoEncrypts = 0;
	fmsgFive->whoDecrypts = 7;
	fmsgFive->subGroupNodeOne = 4;
	fmsgFive->subGroupNodeTwo = 5;
	fmsgFive->subGroupNodeThree = 6;
	fmsgFive->subGroupNodeFour= 7;
	
  packetbuf_clear();
  
	unicast_open(&uc, 146, &unicast_callbacks);
	
  rimeaddr_t addr;
	packetbuf_copyfrom(fmsgFive, 110);
  addr.u8[0] = WhoIsLeader;
  addr.u8[1] = 0;
  if(!rimeaddr_cmp(&addr, &rimeaddr_node_addr))
  {
  	unicast_send(&uc, &addr);
  }
}

static void leader_unicast_node_decypt_request(int WhichNodeWantDecrypt)
{
	struct fulMsg *fmsgSix = &fulMsg_t;
	fmsgSix->msgType = 6;//request
	fmsgSix->whoEncrypts = 0;
	fmsgSix->whoDecrypts = WhichNodeWantDecrypt;
	fmsgSix->subGroupNodeOne = 4;
	fmsgSix->subGroupNodeTwo = 5;
	fmsgSix->subGroupNodeThree = 6;
	fmsgSix->subGroupNodeFour= 7;
	int i = 0; 
	for(i = 6; i >3; i--)
	{
		packetbuf_clear();
		
		unicast_open(&uc, 146, &unicast_callbacks);
	
		rimeaddr_t addr;
		packetbuf_copyfrom(fmsgSix, 110);
		addr.u8[0] = i;
		addr.u8[1] = 0;
		if(!rimeaddr_cmp(&addr, &rimeaddr_node_addr))
		{
			unicast_send(&uc, &addr);
		}	
	}
}

PROCESS_THREAD(unicast_key_dNode, ev, data)
{
	PROCESS_EXITHANDLER(unicast_close(&uc);)
	PROCESS_BEGIN();
		unicast_open(&uc, 146, &unicast_callbacks);
	  static struct etimer et;
		rimeaddr_t addr;
		etimer_set(&et, CLOCK_SECOND);	
		PROCESS_WAIT_EVENT_UNTIL(etimer_expired(&et));		
		//generate key part
		//unicast them to encryption nodeOne
		//NN_Mult (a, b, c, digits)       Computes a = b * c.
		//ecc_gen_private_key(prKey_self);
		//NN_Mult(prKey_selfM, prKey_self, s_s, NUMWORDS);
		//send this to encryption node

		//create the struct and then send

		struct fulMsg *fmsgSeven = &fulMsg_t;
		fmsgSeven->msgType = 7;//request
		fmsgSeven->whoEncrypts = 4;
		fmsgSeven->whoDecrypts = 7;
		fmsgSeven->subGroupNodeOne = 4;
		fmsgSeven->subGroupNodeTwo = 5;
		fmsgSeven->subGroupNodeThree = 6;
		fmsgSeven->subGroupNodeFour= 7;
		int a = 0;
		for(a = 0; a< 6; a++)
		fmsgSeven->sintops[a] = prKey_selfM[a];
	
		packetbuf_clear();

		packetbuf_copyfrom(fmsgSeven, 110);
		addr.u8[0] = 7;
		addr.u8[1] = 0;
		if(!rimeaddr_cmp(&addr, &rimeaddr_node_addr))
		{			
			unicast_send(&uc, &addr);
		}
	PROCESS_END();
}

static void decryptMessage()
{
	int a = 0;
	for(a = 0; a< 6; a++)
	prKey_Dshare7[a] = prKey_selfM[a];
	
	NN_DIGIT Dtemp1[NUMWORDS];
	NN_DIGIT Dtemp2[NUMWORDS];
	NN_DIGIT Dtemp3[NUMWORDS];
	NN_DIGIT Dtemp4[NUMWORDS];

	memset(Dtemp1, 0, NUMWORDS*NN_DIGIT_LEN);
	memset(Dtemp2, 0, NUMWORDS*NN_DIGIT_LEN);
	memset(Dtemp3, 0, NUMWORDS*NN_DIGIT_LEN);
	memset(Dtemp4, 0, NUMWORDS*NN_DIGIT_LEN);

	NN_Add (Dtemp1, prKey_Dshare4, prKey_Dshare5, NUMWORDS);
	NN_Add (Dtemp2, prKey_Dshare6, prKey_Dshare7, NUMWORDS);
	NN_Add (Dtemp3, Dtemp2, Dtemp1, NUMWORDS);
	
	NN_ModInv(Dtemp4, s_s, param.p, NUMWORDS);
  NN_Mult(prKey_group, Dtemp3, Dtemp4, NUMWORDS);
	
/*	int b1 = 0;*/
/*  for(b1 = 0; b1 <6; b1++ )*/
/*  {*/
/*  	printf("on mote 7 pr4 %p\n",Dtemp3[b1]);  */
/*  }*/
/*   */
//NN_Add (a, b, c, digits) Computes a = b + c.
	//NN_ModInv (a, b, c, digits)     Computes a = 1/b mod c
	//NN_Mult (a, b, c, digits)       Computes a = b * c.

	uint8_t *dM = malloc(MAX_M_LEN*sizeof(uint8_t));
  int dM_len = MAX_M_LEN;
	
  //uint8_t *rC = malloc((2*KEYDIGITS*NN_DIGIT_LEN + 1 + MAX_M_LEN + HMAC_LEN)*sizeof(uint8_t));	
  int rC_len =102;
  
  dM_len = decrypt(dM, dM_len, rC, rC_len, prKey_group);
  
  printf("Recieved Temperature: ");
  int ii = 0;
	for (ii = 0; ii != 1; ++ii) 
	{
		if (ii%32 == 0 && ii != 0) printf("\n");
		printf("%02x 'C", dM[ii]);
	}
	printf("\n");
	memset(prKey_group, 0, NUMWORDS*NN_DIGIT_LEN);
}

PROCESS_THREAD(startup_process, ev, data)
{

PROCESS_EXITHANDLER(unicast_close(&uc);)
    
  PROCESS_BEGIN();

  unicast_open(&uc, 146, &unicast_callbacks);
  memset(prKey_group, 0, NUMWORDS*NN_DIGIT_LEN);
  memset(prKey_self, 0, NUMWORDS*NN_DIGIT_LEN);
  memset(prKey_selfM, 0, NUMWORDS*NN_DIGIT_LEN);
  
  memset(prKey_Eshare4, 0, NUMWORDS*NN_DIGIT_LEN);
  memset(prKey_Eshare5, 0, NUMWORDS*NN_DIGIT_LEN);
  memset(prKey_Eshare6, 0, NUMWORDS*NN_DIGIT_LEN);
  memset(prKey_Eshare7, 0, NUMWORDS*NN_DIGIT_LEN);
  
  memset(prKey_Dshare4, 0, NUMWORDS*NN_DIGIT_LEN);
  memset(prKey_Dshare5, 0, NUMWORDS*NN_DIGIT_LEN);
  memset(prKey_Dshare6, 0, NUMWORDS*NN_DIGIT_LEN);
  memset(prKey_Dshare7, 0, NUMWORDS*NN_DIGIT_LEN);
  
  memset(s_s, 0, NUMWORDS*NN_DIGIT_LEN);
  memset(pbkey_group.x, 0, NUMWORDS*NN_DIGIT_LEN);
  memset(pbkey_group.y, 0, NUMWORDS*NN_DIGIT_LEN);

	s_s[5] = 0x00000000;
	s_s[4] = 0xc36e3e96;
	s_s[3] = 0xc26c3d91;
	s_s[2] = 0xc7ec7db1;
	s_s[1] = 0xd7e47933;
	s_s[0] = 0x16020c0d;

/*  Initialize ecc. */
  ecc_init();
  ecies_init();
  get_curve_param(&param);   
  
  ecc_gen_private_key(prKey_self);
	NN_Mult(prKey_selfM, prKey_self, s_s, NUMWORDS); 

  while(1) 
  {
   	PROCESS_WAIT_EVENT_UNTIL(ev == sensors_event && data == &button_sensor);
    if(node_id == 4)
    {
		  node_unicast_leder_encrypt_request();
		  printf("message sent by node 4.\n");
    }
    if(node_id == 7)
    {
		  node_unicast_leder_decrypt_request();
		  printf("message recieved 4.\n");
    }
     if(node_id == 2)
    {
		  node_unicast_leder_decrypt_request();
    }
  }
  PROCESS_END();
}
