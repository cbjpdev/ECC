#include "contiki.h"
#include "net/rime.h"
#include "dev/button-sensor.h"
#include "dev/leds.h"
#include "sys/etimer.h"
#include <stdio.h>
#include <stdlib.h>
#include "node-id.h"

/*---------------------------------------------------------------------------*/
PROCESS(example_unicast_process, "Example unicast");
PROCESS(test_unicast_process, "test unicast");
AUTOSTART_PROCESSES(&example_unicast_process);
/*---------------------------------------------------------------------------*/

static void
recv_uc(struct unicast_conn *c, const rimeaddr_t *from)
{
  printf("---------------------unicast message received from %d.%d\n",from->u8[0], from->u8[1]);
}

static const struct unicast_callbacks unicast_callbacks = {recv_uc};

static struct unicast_conn uc;

/*---------------------------------------------------------------------------*/
static void testFunction()
{
	 process_start(&test_unicast_process, 4);
}

PROCESS_THREAD(test_unicast_process, ev, data)
{
	PROCESS_EXITHANDLER(unicast_close(&uc);)
	PROCESS_BEGIN();
		unicast_open(&uc, 146, &unicast_callbacks);
	//while(1)
	//{
	printf("%d\n",data);
	  static struct etimer et;
		rimeaddr_t addr;
		etimer_set(&et, CLOCK_SECOND);	
		PROCESS_WAIT_EVENT_UNTIL(etimer_expired(&et));
		packetbuf_copyfrom("Hello", 5);
		addr.u8[0] = 10;
		addr.u8[1] = 0;
		if(!rimeaddr_cmp(&addr, &rimeaddr_node_addr))
		{			
			if(node_id == 7 || node_id == 6 || node_id == 5)
			unicast_send(&uc, &addr);
		}
  //}
	PROCESS_END();
}
PROCESS_THREAD(example_unicast_process, ev, data)
{

/*		if (node_id == 5) */
/*		{*/
/*			etimer_set(&et, CLOCK_SECOND*1);*/
/*		}*/
/*		*/
/*		if (node_id == 6)*/
/*		{*/
/*			etimer_set(&et, CLOCK_SECOND*2);*/
/*		}*/
/*		*/
/*		if (node_id == 7)*/
/*		{*/
/*			etimer_set(&et, CLOCK_SECOND*3);*/
/*		}*/
	
	
  PROCESS_BEGIN();
		testFunction();
  PROCESS_END();
}
/*---------------------------------------------------------------------------*/
