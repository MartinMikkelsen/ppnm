CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += -pthread
#CFLAGS += -fopenmp
LDLIBS = -lm #-lgomp

default: out.txt mcplot.png
out.txt: main
	time ./$< >$@

mcplot.png: data.txt Makefile circle.txt
	echo '\
	set terminal png;\
	set size square;\
	set output "$@";\
	set tics out;\
	set key outside top;\
	set xlabel "x-coordinate";\
	set ylabel "y-coordinate";\
	set xrange [0:1];\
	set yrange [0:1];\
	set title "Monte Carlo Method";\
	set style fill transparent solid 0.2 noborder;\
	plot \
	"circle.txt" using 1:2 with filledcurves x1 lc "green" title "Circle r=1"\
	,"$<" using 1:2 with points lc "purple" title "Random points"\
	'|gnuplot

.PHONEY: clean

clean:
	$(RM) main *.txt *.png
