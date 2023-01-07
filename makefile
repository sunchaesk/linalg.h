CC=gcc
WARNING= -Wall
CFLAGS=-g -std=c11
SRCS=$(wildcard *.c)
OBJS=$(SRCS:.c=.o)

c: $(SRCS)
	$(CC) $(WARNING) $(CFLAGS) $(SRCS) -o c
	./c

warn: $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) -o c
	./c

clean:
	rm c
