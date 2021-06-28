CFLAGS = -O -Wall -std=gnu11
LDLIBS = -lm

all: output.txt
	@echo '================================='
	@cat output.txt
	@echo '================================='
	@cat out.txt


output.txt error.txt: input.txt cmdline stdio read-file 
	./cmdline $$(cat input.txt) 1>output.txt 2>error.txt 
	./stdio < input.txt 1> out.file.txt 
	./read-file input.txt out.txt
input.txt: Makefile 
	cat /dev/null > $@
	for x in $$(seq 0 0.5 5); do echo "$$x" >> $@; done

clean: 
	$(RM) input.txt out.file.txt stdio out* error.txt cmdline out.txt
	
