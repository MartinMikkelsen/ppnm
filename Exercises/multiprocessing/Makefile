CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += -pthread
CFLAGS += -fopenmp
LDLIBS = -lm #-lgomp

default: out.txt mcplot.pdf
	@echo 
	@cat out.txt
	@echo
	@echo 'The plot is msplot.pdf'

out.txt: main
	time ./$< >$@

mcplot.pdf: mcplot.ppl
	pyxplot $< 1> /dev/null

mcplot.ppl: data.txt circle.txt Makefile
	@echo 'set terminal pdf dpi 600'						>  $@
	@echo 'set output "mcplot.pdf"'						>> $@
	@echo 'set key outside' 								>> $@
	@echo 'set grid'											>> $@
	@echo 'set xlabel "$$x$$"'								>> $@
	@echo 'set ylabel "$$y$$"'								>> $@
	@echo 'set xrange [0:1]'								>> $@
	@echo 'set yrange [0:1]'								>> $@
	@echo 'set title "Monte Carlo Method "'	>> $@
	@echo 'plot \\'																								>> $@
	@echo '"circle.txt" using 1:2 w l lt 10 lw 4 color BrickRed  title "Circle $$r=1$$" 			\\'>> $@
	@echo ',"data.txt" using 1:2 w p pt 31 lw 4 color grey title "Random points" 			\\'>> $@

.PHONEY: clean

clean:
	$(RM) main *.txt *.png
