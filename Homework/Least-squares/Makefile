CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
CFLAGS  =  -lm -Wall -std=gnu11

CC   	=  gcc
PIPE    =  cat
INPUT   =  data.txt
OUTPUT  =  plot.pdf
EXE   	=  main
OBJ  	=  main.o GramSchmidt.o utilities.o lsfit.o
HDRS    =  GramSchmidt.h utilities.h lsfit.h
IMGFL   =  plot.pdf
TXT     =  output.txt

.PHONEY: default

plot.pdf: plot.ppl
	pyxplot $< 1> /dev/null 

plot.ppl: data.txt Makefile
	@echo 'set terminal pdf dpi 600'						>  $@ 
	@echo 'set output "plot.pdf"'							>> $@           		
	@echo 'set key outside' 								>> $@           		
	@echo 'set grid'										>> $@           		
	@echo 'set xlabel "Time [s]"'						>> $@      		
	@echo 'set ylabel "Activity [percent]"'				>> $@     		
	@echo 'set title "Measurements with error bars"'			>> $@	
	@echo 'plot \\'											>> $@
	@echo ' "data.txt" using 1:2:3 with yerrorbars color black title "Measurements" 	\\'>> $@
	
Out.data.text:main
		./$< 1>$@ 2>log

logdata.txt: main
		./$< 1>$@ Out.data.text


main: main.o GramSchmidt.o utilities.o lsfit.o

clean:
	$(RM) $(TXT) $(OUTPUT) $(OBJ) $(IMGFL)
