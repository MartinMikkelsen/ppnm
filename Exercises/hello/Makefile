CC = gcc # the C compiler
CFLAGS = -O -std=gnu11 # options for the C compiler
LDLIBS = -lm # libraries to link

default: out.txt    # default target is to make out.txt
	cat out.txt # show out.txt on the screen

out.txt: hello            # out.txt needs hello
	./hello > out.txt # run hello, output goes into out.txt

hello: hello.o                           # hello needs hello.o
	$(CC) -o hello hello.o $(LDLIBS) # link hello.o into hello

hello.o: hello.c                   # hello.o needs hello.c
	$(CC) $(CFLAGS) -c hello.c # compile hello.c

.PHONEY: clean
clean:                              # clean is a phoney target
	$(RM) hello.o hello out.txt # clean up the directory

.PHONEY: test
test:                      # test target used for debugging
	echo $(LDLIBS)
	echo $(CC)
	echo $(RM)
