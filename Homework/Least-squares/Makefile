CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
LDLIBS +=  -lm -Wall

CC   	=  gcc
PIPE    =  cat

INPUT   =  data.txt
OUTPUT  =  fit.plt
EXE   	=  main
OBJ  	=  main.o gramSchmidt.o backsub.o input2Array.o leastSquares.o utilities.o
HDRS    =  gramSchmidt.h backsub.h input2Array.h leastSquares.h utilities.h
IMGFL   =  fit.png
TXT     =  output.txt

.PHONEY: default

Out.data.txt:main
		./$< 1>$@ 2>log
main: main.o
qr.o: qr.c
	$(CC) $(CFLAGS) -c -o $@ $<
clean:
	$(RM) $(TXT) $(OUTPUT) $(OBJ) $(IMGFL)
