MISCPATH  = ../dmel_r5.57_FB2014_03/
MISCFILES = act_exons_r5.57.rda act_introns_r5.57.rda
MISC      = $(addprefix $(MISCPATH), $(MISCFILES))

SuppFigX.pdf: $(MISC)
	R -f plot_SuppFigX.R

$(MISC):
	$(MAKE) -C $(MISCPATH)

clean:
	rm SuppFigX.pdf
