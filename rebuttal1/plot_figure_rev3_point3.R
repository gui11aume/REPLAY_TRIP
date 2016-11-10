load("../Fig4/models.rda")

pdf("figure_rev3_point3.pdf", width=8, height=4, useDingbats=FALSE)
par(mfrow=c(1,2))
acf(models[[4]]$residuals, lag=200, main="Terminators")
acf(fullmodel$residuals, lag=200, main="Full model")
dev.off()
