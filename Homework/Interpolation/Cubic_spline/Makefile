CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all: cspline.pdf

cspline.pdf: cspline.ppl
	pyxplot $<

cspline.ppl: cspline.txt datapoints.txt Makefile
	@echo 'set terminal pdf dpi 600'				>  $@
	@echo 'set output "cspline.pdf"'				>> $@
	@echo 'set key outside' 						>> $@
	@echo 'set grid'									>> $@
	@echo 'set xlabel "$$x$$"'						>> $@
	@echo 'set ylabel "$$y$$"'						>> $@
	@echo 'set title "Qubic spline"'				>> $@
	@echo 'plot \\'									>> $@
	@echo ' "datapoints.txt"      using 1:2 w p pt 31 color grey  title "Data points"									\\'>> $@
	@echo ',"cspline.txt" using 1:3 w l lt 1 lw 2 color red   title "interpolation" 		\\'>> $@
	@echo ',"cspline.txt" using 1:2 w l lt 2 lw 4 color black title "gsl interpolation" 	\\'>> $@
	@echo ',"cspline.txt" using 1:4 w l lw 4 lt 2 color black title "gsl derivative" 		\\'>> $@
	@echo ',"cspline.txt" using 1:5 w l lw 2 lt 1 color green title "derivative" 			\\'>> $@
	@echo ',"cspline.txt" using 1:6 w l lw 4 lt 2 color black title "gsl integral" 			\\'>> $@
	@echo ',"cspline.txt" using 1:7 w l lw 2 lt 1 color cyan  title "integral" 			\\'>> $@

datapoints.txt:
	@echo 0	   0  >  $@
	@echo 1	   1  >> $@
	@echo 2	   4  >> $@
	@echo 3	   9  >> $@
	@echo 4	   16  >> $@
	@echo 5	   25 >> $@
	@echo 6	   36 >> $@
	@echo 7	   49 >> $@
	@echo 8	   64 >> $@
	@echo 9	   81 >> $@

cspline.txt: main xs.txt ys.txt
	./$<

main: main.o cspline.o utilities.o

xs.txt ys.txt: datapoints.txt
	cut -d " " -f 1 datapoints.txt > xs.txt  
	cut -d " " -f 2 datapoints.txt > ys.txt  

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.pdf *.png

