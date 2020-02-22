CONTIKI = /home/chap/Documents/project/2contiki-2.6

ifndef TARGET
TARGET=cooja
endif

CONTIKI_ECC_DIR = ..

all: alice bob

PROJECT_SOURCEFILES += messages.c
include $(CONTIKI_ECC_DIR)/Makefile.contikiecc

include $(CONTIKI)/Makefile.include
