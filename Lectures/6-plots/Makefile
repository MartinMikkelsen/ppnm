#CFLAGS = `gsl-config --cflags`
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

erf.gnuplot2.png: data.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot "$<" with line title "error function";\
	' | tee test | gnuplot

data.txt: main
	./$< > $@

main: main.o erf.o

.PHONEY:clean
clean:
	$(RM) *.o main *.png data.txt test
test:
	echo $(CFLAGS)
	echo $(LDLIBS)


