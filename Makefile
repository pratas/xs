BIN = .

CC = gcc
CFLAGS = -g -O2 -Wall

LIBS = -lm
DEPS = defs.h init.h

PROGS = $(BIN)/XS

OBJS = mem.o misc.o msg.o s454.o

all:
	$(MAKE) progs

progs: $(PROGS)

$(BIN)/XS: XS.c $(DEPS) $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/XS XS.c $(OBJS) $(LIBS)

mem.o: mem.c mem.h $(DEPS)
	$(CC) -c $(CFLAGS) mem.c

misc.o: misc.c misc.h $(DEPS)
	$(CC) -c $(CFLAGS) misc.c

msg.o: msg.c msg.h $(DEPS)
	$(CC) -c $(CFLAGS) msg.c

s454.o: s454.c s454.h $(DEPS)
	$(CC) -c $(CFLAGS) s454.c

clean:
	/bin/rm -f *.o $(PROGS)


