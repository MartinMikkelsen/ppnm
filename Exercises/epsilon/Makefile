CC = gcc
CFLAGS = -O #-std=gnull
CFLAGS += -Wall
LDFLAGS = -I.
LDLIBS = -lm

out.txt : epsilon
	time ./$< > $@
epsilon : epsilon.c epsilon2.c
	$(CC) epsilon.c epsilon2.c -o epsilon

clean:
	$(RM) out.txt epsilon

