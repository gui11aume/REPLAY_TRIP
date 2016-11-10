MISC    = ../predpow.rda
MODELS  = models.rda
POTFEAT = potprom.rda potterm.rda potEF1.rda potGSE.rda potenh.rda

all: Fig4a.pdf Fig4b.pdf


$(subst pot, %, $(POTFEAT)):
	R -f makepot.R

$(MODELS): $(POTFEAT)
	R -f fit_models.R

Fig4a.pdf: $(MISC) $(MODELS)
	R -f plot_Fig4a.R

Fig4b.pdf: $(MODELS)
	R -f plot_Fig4b.R

clean:
	rm -f info_Fig4a.txt Fig4*.pdf
