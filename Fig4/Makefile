MISC= ../predpow.rda

all: Fig4a.pdf info_Fig4a.txt

%.pdf info_%.txt : $(MISC)
	R -f plot_Fig4a.R

clean:
	rm -f Fig4a.pdf info_Fig4a.txt
