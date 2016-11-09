#MISCPATH  = ../dmel_r5.57_FB2014_03/
#MISCFILES = exons introns act_promoters
#MISC      = $(addprefix $(MISCPATH), $(addsuffix _r5.57.rda, $(MISCFILES)))
#ACTGENES  = act_genes_r5.57.txt
#SCORES4C  = scores4C.txt

all: gdom0.rda gdom1.rda Fig3a.pdf Fig3b.pdf Fig3c.pdf Fig3d.pdf

gdom%.rda:
	R -f make_domains.R

Fig3a.pdf: 
	R -f plot_Fig3a.R

Fig3b.pdf: gdom0.rda gdom1.rda
	R -f plot_Fig3b.R

Fig3c.pdf: gdom0.rda gdom1.rda
	R -f plot_Fig3c.R

Fig3d.pdf: gdom0.rda gdom1.rda
	R -f plot_Fig3d.R

clean:
	rm -f gdom0.rda gdom1.rda
