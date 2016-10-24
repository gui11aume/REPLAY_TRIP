all: allprom_nochromP.txt allprom.txt gdom0.rda gdom1.rda\
	violin_het.pdf domains.pdf violin_domains.pdf \
	Venn.pdf violin_color_expression.pdf PCA.pdf \
	outliers_heatmap.pdf prdx_domains.txt prdx_idx.txt \
	bun_locus.pdf prdx_vs_RED.pdf modENCODE_spie.pdf \
	barplot_enhancers.pdf violin_short_range.pdf \
	anticorrelation.pdf decay.pdf potEF1.rda potGSE.rda \
	stripsize.pdf


##  ------  DATA  ------  ##
allprom_nochromP.txt allprom.txt:
	# XXX Broken during development XXX 
	R -f generate_tables.R
gdom0.rda gdom1.rda: allprom_nochromP.txt
	R -f make_domains.R


##  ------  DEPENDENCIES  ------  ##
dmel_r5.57_FB2014_03/%.rda:
	make -C dmel_r5.57_FB2014_03


##  ------  FIGURES  ------  ##

# Figure 2a
violin_het.pdf: allprom.txt
	R -f plot_violin_het.R
# Figure 2b
domains.pdf: allprom_nochromP.txt gdom0.rda
	R -f plot_domains.R
# Figure 2c
violin_domains.pdf: allprom.txt gdom0.rda gdom1.rda
	R -f plot_violin_domains.R
# Figure 2d
Venn.pdf: allprom.txt
	R -f plot_Venn_BLACK.R
# Figure 2e
violin_color_expression.pdf: allprom_nochromP.txt
	R -f plot_violin_color_expression.R

# Figure 3a and 3b
PCA.pdf outliers_heatmap.pdf prdx_domains.txt: allprom_nochromP.txt gdom0.rda
	R -f plot_PCA.R
# Genes covered by paradoxical domains.
prdx_idx.txt stripsize.pdf: dmel_r5.57_FB2014_03/genes_r5.57.rda \
		prdx_domains.txt
	R -f prdx_overlap.R
# Figure 3c
bun_locus.pdf:
	R -f plot_bun_locus.R
# Figure 3d
prdx_vs_RED.pdf: prdx_domains.txt
	R -f prdx_vs_RED.R
# Figure 3e
modENCODE_spie.pdf: prdx_domains.txt
	R -f plot_modENCODE_spie.R
# Figure 3f
barplot_enhancers.pdf: prdx_domains.txt
	R -f plot_barplot_enhancers.R

# Figure 4b
violin_short_range.pdf anticorrelation.pdf decay.pdf: prdx_domains.txt
	R -f plot_violin_short_range.R
# Figure 4e
potGSE.rda potEF1.rda:
	R -f makepot.R
