FLYBASE = ftp://flybase.org
GFFPATH = genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/gff
GFFALL  = dmel-all-no-analysis-r5.57.gff.gz

STEMS = \
exons     act_exons     inact_exons     \
genes     act_genes     inact_genes     \
introns   act_introns   inact_introns   \
promoters act_promoters inact_promoters \

TARGETS = $(addsuffix _r5.57.rda, $(STEMS))
PATTARG = $(addsuffix %r5.57.rda, $(STEMS))

.INTERMEDIATE: data/$(GFFALL)

all: $(TARGETS)

data/$(GFFALL):
	wget $(FLYBASE)/$(GFFPATH)/$(GFFALL) -P data

data/dmel-all-exons-r5.57.gff: data/$(GFFALL)
	zcat data/$(GFFALL) | \
		awk '{if ($$2 == "FlyBase" && $$3 == "exon") print}' \
		> data/dmel-all-exons-r5.57.gff

data/dmel-all-genes-r5.57.gff: data/$(GFFALL)
	zcat data/$(GFFALL) | \
		awk '{if ($$2 == "FlyBase" && $$3 == "gene") print}' \
		> data/dmel-all-genes-r5.57.gff

$(PATTARG): data/dmel-all-exons-r5.57.gff data/dmel-all-genes-r5.57.gff
	R -f generate_dmel_r5.57_objects.R

clean:
	rm -f *.rda data/*.gff
