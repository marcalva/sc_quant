
CC = gcc
CFLAGS = -g -O2 -Wall

all : sc_quant

HTSDIR = htslib
include $(HTSDIR)/htslib.mk
include $(HTSDIR)/htslib_static.mk
HTSLIB = $(HTSDIR)/libhts.a
HTSLIB_LIB = $(HTSLIB) $(HTSLIB_static_LIBS)
SCLIBS = $(HTSLIB_LIB) -lpthread

HTSLIB_CPPFLAGS = -I$(HTSDIR)
CPPFLAGS = -I. $(HTSLIB_CPPFLAGS)

OBJS = main.o ac.o a_count.o bc_umi.o bins.o gc.o g_count.o \
	gtf_anno.o overlap.o sam_read.o str_util.o variants.o

LDFLAGS = $(HTSLIB_LIBS)

sc_quant : $(OBJS) $(HTSLIB)
	$(CC) $(LDFLAGS) -o $@ $(OBJS) $(SCLIBS)

%.o: %.c %.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $< 

# for specific configuration but htslib.mk seems to work well
# HTSCONF = --disable-lzma --without-curses
# $(HTSLIB) : 
# 	cd $(HTSDIR) && autoreconf -i && ./configure $(HTSCONF) && make && cd ../

clean:
	rm -f *o
	rm sc_quant
	cd $(HTSDIR) && make clean && cd ../

