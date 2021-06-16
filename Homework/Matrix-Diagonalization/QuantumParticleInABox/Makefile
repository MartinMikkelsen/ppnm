CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all: eigenfunctions.pdf output.txt
 
	
eigenfunctions.txt output.txt: main
	./$< > $@

eigenfunctions.pdf: eigenfunctions.ppl
	pyxplot $< 1> /dev/null 

eigenfunctions.ppl: eigenfunctions.txt Makefile
	@echo 'set terminal pdf dpi 600'						>  $@ 
	@echo 'set output "eigenfunctions.pdf"'					>> $@           		
	@echo 'set key outside' 								>> $@           		
	@echo 'set grid'										>> $@           		
	@echo 'set yrange [-1.1:1.1]'							>> $@   
	@echo 'set xrange [0:1]'								>> $@
	@echo 'set xlabel "$$x$$"'							>> $@      		
	@echo 'set ylabel "$$\psi$$"'							>> $@     		
	@echo 'set title "Particle in a box"'					>> $@
	@echo 'set colormap (5-1)*c1+1'							>> $@
	@echo 'set palette brickRed, limeGreen, cadetBlue'		>> $@
	@echo 'set c1tics -1,0.25'								>> $@ 
	@echo 'plot \\'											>> $@
	@echo '"eigenfunctions.txt" using 1:3 w l lw 3 color red title "$$n=1$$; numerical "	 	\\'>> $@
	@echo ',"eigenfunctions.txt" using 1:5 w l lw 3 color green title "$$n=2$$; numerical" 		\\'>> $@
	@echo ',"eigenfunctions.txt" using 1:7 w l lw 3 color blue title "$$n=3$$; numerical" 		\\'>> $@
	@echo ',"eigenfunctions.txt" using 1:2 w l lw 2 color black title "$$n=1$$" 	\\'>> $@	
	@echo ',"eigenfunctions.txt" using 1:4 w l lw 2 color black title "$$n=2$$" 	\\'>> $@
	@echo ',"eigenfunctions.txt" using 1:6 w l lw 2 color black title "$$n=3$$" 	\\'>> $@

main: main.o GramSchmidt.o utilities.o Jacobi.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png 

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)
