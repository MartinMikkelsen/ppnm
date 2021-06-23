CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)


Hydrogen.pdf: hydrogen.ppl
	pyxplot $< 1> /dev/null

hydrogen.ppl: hydrogen.txt Makefile
	@echo 'set terminal pdf dpi 600'						>  $@
	@echo 'set output "Hydrogen.pdf"'						>> $@
	@echo 'set key outside' 								>> $@
	@echo 'set grid'											>> $@
	@echo 'set xlabel "$$x$$"'								>> $@
	@echo 'set ylabel "$$y$$"'								>> $@
	@echo 'set title r"""Numerical integration of s-wave radial Schr\"odinger equation"""'  >>$@
	@echo 'plot \\'																								>> $@
	@echo '"hydrogen.txt" using 1:2 w l lt 1 lw 2 color BrickRed   title "Numerical" 				\\'>> $@
	@echo ',"hydrogen.txt" using 1:3 w l lt 5 lw 2 color Cyan  title "Analytical" 	\\'>> $@

hydrogen.txt: main
	./$< > $@

main: main.o Runge-Kutta.o NewtonJacobian.o GramSchmidt.o utilities.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png

