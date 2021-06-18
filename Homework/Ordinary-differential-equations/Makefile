CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)


all: solution1.pdf
	@echo

	
solu1.txt: main
	./$< > $@

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

main: main.o utilities.o Runge-Kutta.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png ../util.o

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)
