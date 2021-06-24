CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all: output.txt
	@echo
	@cat output.txt
	@echo 

output.txt: main
	./$< > $@

main: main.o qNewton.o 

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png 
