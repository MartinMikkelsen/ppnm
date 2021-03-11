CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY:default
default: expo.pdf ms.pdf

expo.pdf: data.txt Makefile 
	echo '\
	set terminal pdfcairo;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "$x$";\
	set ylabel "$y$";\
	set title "Exponential function";\
	plot \
        "$<" using 1:2 with line title "expo func"\
        '| gnuplot

data.txt: main 
	./$< >$@

ms.pdf: ms.tex
	pdflatex $<
	pdflatex $<
#	bibtex $<
#	pdflatex $<


.PHONEY:clean
clean:
	$(RM) main 
