#include "nn.h"

struct fulMsg {
	int msgType;
	uint8_t *C;	
	int whoEncrypts;
	int whoDecrypts;
	int subGroupNodeOne;
	int subGroupNodeTwo;
	int subGroupNodeThree;
	int subGroupNodeFour;
	NN_DIGIT sintops[NUMWORDS]; 
} fulMsg_t;

typedef struct fulMsg fulMsg_tt;
