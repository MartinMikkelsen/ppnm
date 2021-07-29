CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all:fft.pdf time.pdf


fft.pdf: fft.ppl
	pyxplot $< 1> /dev/null

fft.ppl: out.noise.data Makefile
	@echo 'set terminal pdf dpi 600'						>  $@
	@echo 'set output "fft.pdf"'						>> $@
	@echo 'set key outside' 								>> $@
	@echo 'set grid'											>> $@
	@echo 'set xlabel "time"'								>> $@
	@echo 'set ylabel "$$y$$"'								>> $@
	@echo 'set title "Noise filtering using DFT"'	>> $@
	@echo 'plot \\'																								>> $@
	@echo '"out.noise.data"  index 1 w l lw 2 color cyan title "Signal with random noise" 	\\'>> $@
	@echo ',"out.noise.data"  index 0 w l lw 2 color black title "Original signal" 	\\'>> $@
	@echo ',"out.noise.data"  index 2 w l lt 5 lw 2 color PineGreen title "Amplitude filtered" 	\\'>> $@
	@echo ',"out.noise.data"  index 3 w l lt 5 lw 2 color Brickred title "Frequency filtered" 	\\'>> $@

time.pdf: time.ppl
	pyxplot $< 1> /dev/null

time.ppl: out.timing.data Makefile
	@echo 'set terminal pdf dpi 600'						>  $@
	@echo 'set output "time.pdf"'						>> $@
	@echo 'set key outside' 								>> $@
	@echo 'set grid'											>> $@
	@echo 'set xlabel "data set size $$N$$"'								>> $@
	@echo 'set ylabel "time to perform FFT"'								>> $@
	@echo 'set title "FFT time as function of data set size"'	>> $@
	@echo 'plot \\'																								>> $@
	@echo '"out.timing.data" w p pt 2 color cyan title "Data" 	\\'>> $@
	@echo ',1*x*log(x) w l lt 2 lw 2 color SpringGreen title "N*log(N)"  	\\'>> $@



N = $(shell echo '3*2^12; 3*2^13; 3*2^14; 3*2^15'|bc -l)

out.timing.data:main Makefile
	>$@
	for n in $(N); do time -f "$$n %U" -ao $@ ./$< $$n; done

out.noise.data:main
	./$< > $@

main: main.o fft.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png
