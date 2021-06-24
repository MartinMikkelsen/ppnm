CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)


fit.pdf: fit.ppl
	pyxplot $< 1> /dev/null

fit.ppl: out.txt Makefile
	@echo 'set terminal pdf dpi 600'						>  $@
	@echo 'set output "fit.pdf"'						>> $@
	@echo 'set key outside' 								>> $@
	@echo 'set grid'											>> $@
	@echo 'set xlabel "$$x$$"'								>> $@
	@echo 'set ylabel "$$y$$"'								>> $@
	@echo 'set title "Neural Network curve fitting "'	>> $@
	@echo 'plot \\'																								>> $@
	@echo '"out.txt"  index 1 using 1:2 w l lw 2 color Cyan  title "Fit" 	\\'>> $@
	@echo ',"out.txt" index 0 using 1:2 w p pt 31 color BrickRed title "Data points" 				\\'>> $@
	@echo ',"out.txt" index 1 using 1:3 w l lw 2 color Green title "Derivative" 				\\'>> $@
	@echo ',"out.txt" index 1 using 1:4 w l lw 2 color Black title "Integral" 				\\'>> $@




out.txt:main
	./$< 1> $@ 2>log

main: main.o qNewton.o ann.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.log
