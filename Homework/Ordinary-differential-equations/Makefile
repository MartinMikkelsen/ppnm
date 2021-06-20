CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all: solution1.pdf SIR_model_1.pdf SIR_model_2.pdf SIR_model_3.pdf
	@echo
	@echo 'Changeing Tc will drastically impact the behavoir of the functions'


solution1.pdf: solution1.ppl
	pyxplot $< 1> /dev/null

solution1.ppl: solu1.txt Makefile
	@echo 'set terminal pdf dpi 600'						>  $@
	@echo 'set output "solutions1.pdf"'						>> $@
	@echo 'set key outside' 								>> $@
	@echo 'set grid'											>> $@
	@echo 'set yrange [-1.1:1.1]'							>> $@
	@echo 'set xlabel "$$x$$"'								>> $@
	@echo 'set ylabel "$$y$$"'								>> $@
	@echo 'set title "Numerical integration of $$y^{\prime\prime} = 2(x²/2-1/2)y$$ with $$y(0) = 1, y^{\prime}(0) = 0$$ "'	>> $@
	@echo 'plot \\'																								>> $@
	@echo '"solu1.txt" using 1:2 w l lw 1 color red   title "$$y_0 = u$$" 				\\'>> $@
	@echo ',"solu1.txt" using 1:3 w l lw 1 color cyan  title "$$y_1 = u^{\prime}$$" 	\\'>> $@


SIR_model_1.pdf: SIR_model_1.ppl
	pyxplot $< 1> /dev/null

SIR_model_1.ppl: sir_model_1.txt Makefile
	@echo 'set terminal pdf dpi 600'						>  $@
	@echo 'set output "SIR_model_1.pdf"'							>> $@
	@echo 'set key outside' 								>> $@
	@echo 'set grid'											>> $@
	@echo 'set xlabel "$$t \, \mathrm{[days]}$$"'			>> $@
	@echo 'set ylabel "$$N$$"'								>> $@
	@echo 'set title "Epidemic development with $$N = 5.806 \\times 10^6, \; T_r = 2.0, \; T_c = 1.0, \; N_{i,0} = 15$$ "'	>> $@
	@echo 'plot \\'																								>> $@
	@echo ' "sir_model_1.txt" using 1:2 w l lw 2 color red title "Susceptible $$(S)$$" 	\\'>> $@
	@echo ',"sir_model_1.txt" using 1:3 w l lw 2 color green title "Infectious $$(I)$$" 	\\'>> $@
	@echo ',"sir_model_1.txt" using 1:4 w l lw 2 color blue title "Removed $$(R)$$" 	 	\\'>> $@

SIR_model_2.pdf: SIR_model_2.ppl
		pyxplot $< 1> /dev/null

SIR_model_2.ppl: sir_model_2.txt Makefile
	@echo 'set terminal pdf dpi 600'						>  $@
	@echo 'set output "SIR_model_2.pdf"'							>> $@
	@echo 'set key outside' 								>> $@
	@echo 'set grid'											>> $@
	@echo 'set xlabel "$$t \, \mathrm{[days]}$$"'			>> $@
	@echo 'set ylabel "$$N$$"'								>> $@
	@echo 'set title "Epidemic development with $$N = 5.806 \\times 10^6, \; T_r = 6.0, \; T_c = 2.0, \; N_{i,0} = 15$$ "'	>> $@
	@echo 'plot \\'																								>> $@
	@echo ' "sir_model_2.txt" using 1:2 w l lw 2 color red title "Susceptible $$(S)$$" 	\\'>> $@
	@echo ',"sir_model_2.txt" using 1:3 w l lw 2 color green title "Infectious $$(I)$$" 	\\'>> $@
	@echo ',"sir_model_2.txt" using 1:4 w l lw 2 color blue title "Removed $$(R)$$" 	 	\\'>> $@

SIR_model_3.pdf: SIR_model_3.ppl
			pyxplot $< 1> /dev/null

SIR_model_3.ppl: sir_model_3.txt Makefile
	@echo 'set terminal pdf dpi 600'						>  $@
	@echo 'set output "SIR_model_3.pdf"'							>> $@
	@echo 'set key outside' 								>> $@
	@echo 'set grid'											>> $@
	@echo 'set xlabel "$$t \, \mathrm{[days]}$$"'			>> $@
	@echo 'set ylabel "$$N$$"'								>> $@
	@echo 'set title "Epidemic development with $$N = 5.806 \\times 10^6, \; T_r = 8.0, \; T_c = 2.0, \; N_{i,0} = 15$$ "'	>> $@
	@echo 'plot \\'																								>> $@
	@echo ' "sir_model_3.txt" using 1:2 w l lw 2 color red title "Susceptible $$(S)$$" 	\\'>> $@
	@echo ',"sir_model_3.txt" using 1:3 w l lw 2 color green title "Infectious $$(I)$$" 	\\'>> $@
	@echo ',"sir_model_3.txt" using 1:4 w l lw 2 color blue title "Removed $$(R)$$" 	 	\\'>> $@

solu1.txt sir_model_1.txt sir_model_2.txt sir_model_3.txt: main
	./$<  > $@


main: main.o utilities.o Runge-Kutta.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png ../util.o
