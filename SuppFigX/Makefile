MISCPATH  = ../dmel_r5.57_FB2014_03/
MISCFILES = act_exons act_introns exons introns
MISC      = $(addprefix $(MISCPATH), $(addsuffix _r5.57.rda, $(MISCFILES)))

all: SuppFigXa.pdf SuppFigXb.pdf

SuppFigXa.pdf: $(MISC)
	R -f plot_SuppFigXa.R

SuppFigXb.pdf: $(MISC)
	R -f plot_SuppFigXb.R

$(MISC):
	$(MAKE) -C $(MISCPATH)

clean:
	rm SuppFigXa.pdf SuppFigXb.pdf
