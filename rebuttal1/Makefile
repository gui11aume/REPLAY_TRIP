MISCPATH  = ../dmel_r5.57_FB2014_03/
MISCFILES = exons introns act_promoters
MISC      = $(addprefix $(MISCPATH), $(addsuffix _r5.57.rda, $(MISCFILES)))
ACTGENES  = act_genes_r5.57.txt
SCORES4C  = scores4C.txt
MODELS    = ../Fig4/models.rda

all: info_rev1_point2.txt info_rev1_point4.txt info_rev2_point2.txt \
	info_rev2_point3.txt info_rev2_point5.txt figure_rev3_point3.pdf 

info_rev1_point2.txt: $(MISC)
	R -f get_info_rev1_point2.R

info_rev1_point4.txt: $(MISC)
	R -f get_info_rev1_point4.R

info_rev2_point2.txt: $(MODELS)
	R -f get_info_rev2_point2.R

info_rev2_point3.txt: $(SCORES4C)
	R -f get_info_rev2_point3.R

info_rev2_point5.txt:
	R -f get_info_rev2_point5.R

figure_rev3_point3.pdf: $(MODELS)
	R -f plot_figure_rev3_point3.R

$(MODELS):
	$(MAKE) -C ../Fig4

$(ACTGENES): $(MISC)
	R -f make_act_genes.R

$(SCORES4C): $(ACTGENES)
	python create_2kb_term_vectors.py > $(SCORES4C)

$(MISC):
	$(MAKE) -C $(MISCPATH)

clean:
	rm -f info_rev*.txt figure_rev*.pdf
