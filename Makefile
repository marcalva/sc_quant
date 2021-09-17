
CC = gcc
CFLAGS = -g -O2 -pthread -Wall -lhts

# SRCS = main.c ac.c a_count.c bc_umi.o bins.o gc.o g_count.o \
# 	   gtf_anno.c overlap.c sam_read.c str_util.c variants.c

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

sc_quant : main.o ac.o a_count.o bc_umi.o bins.o gc.o g_count.o \
	gtf_anno.o overlap.o sam_read.o str_util.o variants.o
	$(CC) -o sc_quant $^ $(CFLAGS)


