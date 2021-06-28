CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY:default

default: expo.pdf ms.pdf


expo.pdf: expo.ppl
	pyxplot $< 1> /dev/null

expo.ppl: data.txt Makefile
	@echo 'set terminal pdf dpi 600'						>  $@
	@echo 'set output "expo.pdf"'						>> $@
	@echo 'set key outside' 								>> $@
	@echo 'set grid'											>> $@
	@echo 'set xlabel "$$x$$"'								>> $@
	@echo 'set ylabel "$$y$$"'								>> $@
	@echo 'set title "Exponential function "'	>> $@
	@echo 'plot \\'																								>> $@
	@echo '"data.txt"  using 1:2 w l lt 10 lw 4 color BrickRed  title "Exponential function" 			\\'>> $@

data.txt: main 
	./$< >$@

ms.pdf: ms.tex
	pdflatex $<
	pdflatex $<
#	bibtex $<
#	pdflatex $<


.PHONEY:clean
clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png
 
