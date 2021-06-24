CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all: out.txt
	@echo
	@cat out.txt
	@echo 

out.txt: main
	./$< > $@

main: main.o qNewton.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png 
