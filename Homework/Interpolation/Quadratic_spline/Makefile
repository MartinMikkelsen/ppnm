CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all: qspline.pdf
	@echo '15 data points made with sin(i*i)'


output.txt: main
	./$< > $@

clean:
	$(RM) Out* main

qspline.pdf: qspline.ppl
	pyxplot $< 1> /dev/null

qspline.ppl: output.txt Makefile
	@echo 'set terminal pdf dpi 600'						>  $@
	@echo 'set output "qspline.pdf"'						>> $@
	@echo 'set key outside' 								>> $@
	@echo 'set grid'											>> $@
	@echo 'set xlabel "$$x$$"'								>> $@
	@echo 'set ylabel "$$y$$"'								>> $@
	@echo 'set title "Quadratic Interpolation "'	>> $@
	@echo 'plot \\'																								>> $@
	@echo '"output.txt" index 0 using 1:2 w p pt 31 ps 0.75 color grey title "Data points" 		\\'>> $@
	@echo ',"output.txt" index 1 using 1:2 w l lt 1 lw 1 ps 0.75 color black  title "Interpolation" 		\\'>> $@
	@echo ',"output.txt" index 2 using 1:2 w l lt 5 lw 1 ps 0.75 color Cyan  title "Integral" 		\\'>> $@
	@echo ',"output.txt" index 3 using 1:2 w l lt 6 lw 1 ps 0.75 color Mahogany  title "Derivative" 		\\'>> $@

main: main.o utilities.o qspline.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png
