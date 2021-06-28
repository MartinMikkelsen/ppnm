CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY:default

default: out.txt

out.txt:main
	./$< > $@

.PHONEY: clean

clean:
	$(RM) main out.txt
