CFLAGS = -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

Out.data.txt: main
	./$< >$@

clean:
	$(RM) main
