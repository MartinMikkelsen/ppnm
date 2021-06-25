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
	@echo 'set title "Neural Network curve fitting, derivative and anti-derivative "'	>> $@
	@echo 'plot \\'																								>> $@
	@echo '"out.txt"  index 1 using 1:2 w l lw 2 color Cyan  title "Fit" 	\\'>> $@
	@echo ',"out.txt" index 1 using 1:3 w l lw 2 color PineGreen title "Derivative" 				\\'>> $@
	@echo ',"out.txt" index 1 using 1:4 w l lw 2 color Black title "Integral" 				\\'>> $@
	@echo ',exp(-x*x)*(-5*sin(5*x-1)-2*x*cos(5*x-1)) w l lt 2 lw 2 color SpringGreen title "Analytical derivative" \\' >> $@
	@echo ',-0.1877184794-0.045046858263+cos(1.0)*x+5.0/2.0*sin(1.0)*x**2-9.0/2.0*cos(1.0)*x**3-155.0/24.0*sin(1.0)*x**4+937.0/120.0*cos(1.0)*x**5 w l lt 5 lw 2 color grey title "Analytical Taylor expansion to 5th order " \\' >> $@
	@echo ',"out.txt" index 0 using 1:2 w p pt 31 color BrickRed title "Data points" 				\\'>> $@



out.txt:main
	./$< 1> $@ 2>log

main: main.o qNewton.o ann.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.log
