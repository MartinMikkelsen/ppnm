CFLAGS = -O -Wall -std=gnu11
LDLIBS = -lm


output.txt error.txt: input.txt 
	./cmdline $$(cat input.txt) 1>output.txt 2>error.txt 
	
input.txt: Makefile 
	cat /dev/null > $@
	for x in $$(seq 0 0.5 5); do echo "$$x" >> $@; done

clean: 
	$(RM) input.txt out* error.txt cmdline
