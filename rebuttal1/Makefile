MISCPATH  = ../dmel_r5.57_FB2014_03/
MISCFILES = exons introns act_promoters
MISC      = $(addprefix $(MISCPATH), $(addsuffix _r5.57.rda, $(MISCFILES)))
ACTGENES  = act_genes_r5.57.txt
SCORES4C  = scores4C.txt

all: info_rev1_point2.txt info_rev1_point4.txt info_rev2_point3.txt \
	info_rev2_point5.txt

info_rev1_point2.txt: $(MISC)
	R -f get_info_rev1_point2.R

info_rev1_point4.txt: $(MISC)
	R -f get_info_rev1_point4.R

info_rev2_point3.txt: $(SCORES4C)
	R -f get_info_rev2_point3.R

info_rev2_point5.txt:
	R -f get_info_rev2_point5.R

$(ACTGENES): $(MISC)
	R -f make_act_genes.R

$(SCORES4C): $(ACTGENES)
	python create_2kb_term_vectors.py > $(SCORES4C)

$(MISC):
	$(MAKE) -C $(MISCPATH)

clean:
	rm -f info_rev1_point2.txt info_rev1_point4.txt info_rev2_point5.txt
