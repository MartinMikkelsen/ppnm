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

main: main.o SSA.o Recursive_Adaptive_Integrator.o Clenshaw–Curtis.o pMonteCarlo.o qMonteCarlo.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png

