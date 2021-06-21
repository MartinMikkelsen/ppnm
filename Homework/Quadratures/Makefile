CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all: output.txt
	@echo
	@cat output.txt
	@echo 
	@echo 'Compare result number 2 to result number 5 to see the effect of the Clenshaw-Curtis transformation.'
	@echo 'The result with the Clenshaw-Curtis transformation is more accurate'
	@echo 'Compared to the GSL integration routines the accuracy is much smaller'
output.txt: main
	./$< > $@

main: main.o Recursive_Adaptive_Integrator.o Clenshaw–Curtis.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png 
