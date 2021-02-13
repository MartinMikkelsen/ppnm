CFLAGS = -Wall -std=gnu99
LDLIBS = -lm
.PHONEY: all clean
all              : out.txt           ; cat $<
out.txt          : main              ; ./$< > $@
main             : main.o komplex.o  # built-in linking rule is used here
main.o komplex.o : komplex.h         # built-in compiling rule is used here
clean            :                   ; $(RM) main main.o komplex.o out.txt        
