FLAGS = -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)


out.myfunc.txt:myfunc
	./$< > $@
#	./myfunc > out.myfunc.txt

out.myfunc.txt: myerror
	./$< > $@
#	./myeror > out.myerror.txt

.PHONEY:clean
clean: $(RM) out* myfunc myerror
