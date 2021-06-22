CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)


all: linspline.pdf


output.txt: main
	./$< > $@

clean:
	$(RM) Out* main

linspline.pdf: linspline.ppl
	pyxplot $< 1> /dev/null

linspline.ppl: output.txt Makefile
	@echo 'set terminal pdf dpi 600'						>  $@
	@echo 'set output "lspline.pdf"'						>> $@
	@echo 'set key outside' 								>> $@
	@echo 'set grid'											>> $@
	@echo 'set xlabel "$$x$$"'								>> $@
	@echo 'set ylabel "$$y$$"'								>> $@
	@echo 'set title "Linear Interpolation "'	>> $@
	@echo 'plot \\'																								>> $@
	@echo '"output.txt" index 2 using 1:2 w l lt 10 lw 4 color BrickRed  title "Linear spline" 			\\'>> $@
	@echo ',"output.txt" index 3 using 1:2 w l lt 10 lw 4 color OliveGreen  title "Integral spline" 		\\'>> $@
	@echo ',"output.txt" index 4 using 1:2 w l lt 5 lw 2 color Sepia  title "GSL linear spline" 		\\'>> $@
	@echo ',"output.txt" index 5 using 1:2 w l lt 5 lw 2 color YellowGreen  title "GSL integral-spline" 	\\'>> $@
	@echo ',"output.txt" index 1 using 1:2 w p pt 31 ps 0.75 color grey title "Exact integral" 	\\'>> $@
	@echo ',"output.txt" index 0 using 1:2 w p pt 31 ps 0.75 color white  title "Data points" 		\\'>> $@

main: main.o utilities.o lspline.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png
