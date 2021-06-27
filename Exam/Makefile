CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

output.txt: main
	./$< > $@

main: main.o SSA.c

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png
