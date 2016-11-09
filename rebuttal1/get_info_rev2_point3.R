o = read.delim("scores4C.txt")

y = o$nexp

x1 = log(o$prom + .01); z1 = x1^2
x2 = log(o$term + .01); z2 = x2^2

mod1 = lm(y ~ x1+z1)
predpow1 = var(mod1$fitted.values) /
   (var(mod1$residuals + mod1$fitted.values))
AIC1 = AIC(mod1)

mod2 = lm(y ~ x2+z2)
predpow2 = var(mod2$fitted.values) /
   (var(mod2$residuals + mod2$fitted.values))
AIC2 = AIC(mod2)

# Start writing to file.
sink(file("info_rev2_point3.txt"))

cat("Predictive power promoters on 4C\n")
predpow1
cat(paste("AIC: ", AIC1, "\n", sep=""))
cat("Predictive power terminators on 4C\n")
predpow2
cat(paste("AIC: ", AIC2, "\n", sep=""))

# End writing to file.
sink(NULL)
