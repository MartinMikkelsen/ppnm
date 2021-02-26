CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY:default
default:out.main.txt out.params.txt out.mygamma.png
out.main.txt:main
	./$< > $@
out.params.txt:param
	./$< > $@

out.mygamma.png: out.mygamma.txt Makefile
	graph \
	--output-format png --bitmap-size 1024x768 \
	--x-label "x" --y-label "y" \
	--top-label "Gamma function" \
	--line-mode 1 $< \
	> $@

out.mygamma.txt:mygamma
	./$< > $@

.PHONEY:clean
clean:
	$(RM) out* main params mygamma


.PHONEY:test
test:
	echo CFLAGS
	echo $(CFLAGS)
	echo $(LDLIBS)
