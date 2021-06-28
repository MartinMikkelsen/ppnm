CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc
.PHONEY: default

defualt: out.txt erfplot.pdf Besselplot.pdf

out.txt: main
	./$< >$@


erfplot.pdf: erfplot.ppl
	pyxplot $< 1> /dev/null

erfplot.ppl: DataB.txt Makefile
	@echo 'set terminal pdf dpi 600'						>  $@
	@echo 'set output "erfplot.pdf"'						>> $@
	@echo 'set key outside' 								>> $@
	@echo 'set grid'											>> $@
	@echo 'set xlabel "$$x$$"'								>> $@
	@echo 'set ylabel "$$y$$"'								>> $@
	@echo 'set title "Error function "'	>> $@
	@echo 'plot \\'																								>> $@
	@echo '"DataB.txt" using 1:2 w l lw 2 color red  title "Error function" 				\\'>> $@

Besselplot.pdf: Besselplot.ppl
	pyxplot $< 1> /dev/null

Besselplot.ppl: DataBessel.txt Makefile
	@echo 'set terminal pdf dpi 600'						>  $@
	@echo 'set output "Besselplot.pdf"'							>> $@
	@echo 'set key outside' 								>> $@
	@echo 'set grid'											>> $@
	@echo 'set xlabel "$$x$$"'			>> $@
	@echo 'set ylabel "$$y$$"'								>> $@
	@echo 'set title "Bessel "'	>> $@
	@echo 'plot \\'																								>> $@
	@echo ' "DataBessel.txt" using 1:2 w l lw 2 color red title "$$J_0$$" 	\\'>> $@



.PHONEY: clean

clean:
	$(RM) main *.txt *png
