MISCPATH  = ../dmel_r5.57_FB2014_03/
MISCFILES = exons introns act_promoters
MISC      = $(addprefix $(MISCPATH), $(addsuffix _r5.57.rda, $(MISCFILES)))

all: info_rev1_point2.txt

info_rev1_point2.txt: $(MISC)
	R -f get_info_rev1_point2.R

$(MISC):
	$(MAKE) -C $(MISCPATH)

clean:
	rm -f info_rev1_point2.txt
