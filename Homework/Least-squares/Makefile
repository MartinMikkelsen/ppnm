CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
CFLAGS  =  -lm -Wall -std=gnu11

CC   	=  gcc

default: output.txt fit.pdf
	cat $<


fit.pdf: fit.ppl
	pyxplot $< 1> /dev/null

fit.ppl: fitparams.txt Makefile
	@echo 'set terminal pdf dpi 600'						>  $@
	@echo 'set output "fit.pdf"'						>> $@
	@echo 'set key outside' 								>> $@
	@echo 'set grid'											>> $@
	@echo 'set xlabel "$$\log(t)$$"'								>> $@
	@echo 'set ylabel "$$\log(y)$$"'								>> $@
	@echo 'set title "Fit with errors"'	>> $@
	@echo 'plot \\'																								>> $@
	@echo '"fitparams.txt" index 1 using 1:2 w l lt 1 lw 2 color OliveGreen  title "Fit" 		\\'>> $@
	@echo ',"fitparams.txt" index 2 using 1:2 w l lt 5 lw 2 color Cyan  title "Fit$$+5\%$$ " 		\\'>> $@
	@echo ',"fitparams.txt" index 3 using 1:2 w l lt 5 lw 2 color Cyan  title "Fit$$-5\%$$ " 		\\'>> $@
	@echo ',"fitparams.txt" index 0 using 1:2:3 with yerrorbars color BrickRed  title "Data" 			\\'>> $@

output.txt: main
	./$< > $@

main: main.o GramSchmidt.o utilities.o lsfit.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png
