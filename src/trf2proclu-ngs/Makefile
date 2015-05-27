CC=gcc
CFLAGS=-O2 -Wall -m64
CFLAGS_DEBUG=-Wall -m64 -ggdb
LD=$(CC)
LDFLAGS=-lm
PROGRAM_NAME=trf2proclu-ngs

$(PROGRAM_NAME): $(PROGRAM_NAME).o
	$(LD) $(LDFLAGS) -o $(PROGRAM_NAME) $(PROGRAM_NAME).o

$(PROGRAM_NAME).debug.o: $(PROGRAM_NAME).c
	$(CC) $(CFLAGS_DEBUG) -o $(PROGRAM_NAME).debug.o -c $(PROGRAM_NAME).c

$(PROGRAM_NAME).o: $(PROGRAM_NAME).c
	$(CC) $(CFLAGS) -c $(PROGRAM_NAME).c

debug: $(PROGRAM_NAME).debug.o
	$(LD) $(LDFLAGS) -o $(PROGRAM_NAME) $(PROGRAM_NAME).debug.o

all: $(PROGRAM_NAME)

install:
	cp $(PROGRAM_NAME) ~/programs/bin/

clean:
	rm $(PROGRAM_NAME) $(PROGRAM_NAME).o $(PROGRAM_NAME).debug.o
