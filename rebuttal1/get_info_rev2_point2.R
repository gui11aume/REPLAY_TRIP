load("../Fig4/models.rda")
modenh = models[[5]]

# Start writing to file.
sink(file("info_rev2_point2.txt"))

cat("Predictive power H3K27ac/H3K4me1\n")
var(modenh$fitted.values) / (var(modenh$residuals + modenh$fitted.values))

# End writing to file.
sink(NULL)
