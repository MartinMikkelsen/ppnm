CFLAGS = `gsl-config --cflags`
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

default: data.txt erf.gnuplot.pdf gamma.gnuplot.pdf

data.txt: main
	./$< > $@

erf.gnuplot.pdf: erf.gnuplot.ppl
	pyxplot $< 

erf.gnuplot.ppl: data.txt tabulatedvalues.txt Makefile 
	@echo 'set terminal pdf dpi 600'						>  $@
	@echo 'set output "erf.gnuplot.pdf"'						>> $@
	@echo 'set key outside' 								>> $@
	@echo 'set grid'											>> $@
	@echo 'set xlabel "$$x$$"'								>> $@
	@echo 'set ylabel "$$y$$"'								>> $@
	@echo 'set title "Error function "'	>> $@
	@echo 'plot \\'																								>> $@
	@echo '"data.txt" using 1:2 w l lt 1 lw 3 color BrickRed  title "erf from math.h" 			\\'>> $@
	@echo ',"data.txt" using 1:3 w l lt 2 lw 3 color Cyan  title "gsl-sf-erf" 			\\'>> $@
	@echo ',"data.txt" using 1:4 w l lt 5 lw 3 color Green  title "Homemade erf" 			\\'>> $@
	@echo ',"tabulatedvalues.txt" using 1:2 w p pt 31 ps 0.4 color OliveGreen  title "Integral spline" 		\\'>> $@

gamma.gnuplot.pdf: gamma.gnuplot.ppl
	pyxplot $< 

gamma.gnuplot.ppl: data2.txt Makefile tabulated2.txt  
	@echo 'set terminal pdf dpi 600'						>  $@
	@echo 'set output "gamma.gnuplot.pdf"'						>> $@
	@echo 'set key outside' 								>> $@
	@echo 'set grid'											>> $@
	@echo 'set xlabel "$$x$$"'								>> $@
	@echo 'set ylabel "$$y$$"'								>> $@
	@echo 'set title "Gamma function"'	>> $@
	@echo 'plot \\'																								>> $@
	@echo '"data2.txt" using 1:2 w l lt 1 lw 3 color BrickRed  title "Gamma from math.h" 			\\'>> $@
	@echo ',"data2.txt" using 1:3 w l lt 2 lw 3 color Cyan  title "gsl-sf-gamma" 			\\'>> $@
	@echo ',"data2.txt" using 1:4 w l lt 5 lw 3 color Green  title "Stirling" 			\\'>> $@
	@echo ',"tabulated2.txt" using 1:2 w p pt 31 ps 0.4  color OliveGreen  title "Tabulated data" 		\\'>> $@


main: main.o erf.o gamma.o 

.PHONEY:clean
clean: 
	$(RM) *.o main *.png data.txt test

