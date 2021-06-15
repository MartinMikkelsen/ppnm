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

logplotfit.pdf: logplotfit.ppl
	pyxplot $< 1> /dev/null

logplotfit.ppl: Out.data.txt Makefile
	@echo 'set terminal pdf dpi 600'						>  $@
	@echo 'set output "logplotfit.pdf"'						>> $@
	@echo ' "Out.data.txt.txt" using 1:2 w l lw 1 lt 5 color black title "$$\ln(y) = (%.2f \pm %.2f) + (%.3f \pm %.3f)t$$"%(lna, dlna, lambda, dlambda) \\'>> $@
	@echo ',"Out.data.txt" using 1:3 w l lw 1 lt 3 color black title "" 					\\'>> $@
	@echo ',"Out.data.txt" using 1:4 w l lw 1 lt 3 color black title "" 					\\'>> $@
	@echo ',"Out.data..txt" using 1:2:3 with yerrorbars color blue title "Measurements" 	\\'>> $@




Out.data.text:main
		./$< 1>$@ 2>log

main: main.o GramSchmidt.o utilities.o lsfit.o

clean:
	$(RM) $(TXT) $(OUTPUT) $(OBJ) $(IMGFL)
