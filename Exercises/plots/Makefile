#CFLAGS = `gsl-config --cflags`
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)
default: erf.gnuplot.png gamma.gnuplot.png

erf.gnuplot.png: data.txt tabulatedvalues.txt Makefile 
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "erf from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	,"$<" using 1:4 with line title "homemade Erf"\
	,"tabulatedvalues.txt" using 1:2 with points title "tabulated data"\
	' | gnuplot

gamma.gnuplot.png: data2.txt Makefile   
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Gamma function";\
	plot \
         "$<" using 1:2 with line title "gamma from math.h"\
         ,"$<" using 1:3 with line title "gsl-sf-gamma"\
	 ,"$<" using 1:4 with line title "Stirling"\
        ,"tabulated2.txt" using 1:2 with points title "tabulated data"\
        '| gnuplot
data.txt: main
	./$< > $@

main: main.o erf.o gamma.o 

.PHONEY:clean
clean: 
	$(RM) *.o main *.png data.txt test

