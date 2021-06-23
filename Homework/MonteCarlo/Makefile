CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)


all: output.txt
	@echo
	@cat output.txt

output.txt: main
	./$< > $@

main: main.o MonteCarlo.o Recursive_Adaptive_Integrator.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png 
