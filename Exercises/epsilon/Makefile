CC = gcc
CFLAGS = -Wall -fwrapv -O

.PHONEY: default
default: out.txt
	cat $<

out.txt: main
	./$< > $@

main: main.c tau.c
	$(CC) main.c tau.c -o main 
.PHONEY: clean

clean:
	$(RM) main out.txt
