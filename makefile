CC = gcc
MPICC = mpic++
CFLAGS = -O3
LDFLAGS = -Wl,-R,\.
LIBS = -lm
UNAME_S = $(shell uname -s)

ifneq (, $(findstring SunOS, $(UNAME_S)))
	LIBS += -lnsl -lsocket -lresolv
endif

compile-parallel:
	$(MPICC) $(CFLAGS) -o BorgParICOW.exe iCOW_BORG_2020_05_18.cpp borgms.c mt19937ar.c $(LIBS)


.PHONY: compile, compile-parallel

