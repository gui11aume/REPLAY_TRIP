# Read in the insertion data.
a = read.delim("../allprom_nochromP.txt", as.is=TRUE)
a = subset(a, prom != "p0")

# Use HMM to segment the genome in up/down domains.
library(HMMt)
b = bridge(a[,c(2,4,4,5)])
fit = BaumWelchT(b$x, b$series.length)

refstates = fit@ViterbiPath[b$nonvirtuals]

# Function used to fit models separately.
fit_models_separately = function() {
   indstates = rep(NA, length(refstates))
   promoters = c("pI", "pII", "pIII", "pIV")
   for (i in 1:4) {
      thisprom = promoters[i]
      thisa = subset(a, prom == thisprom)
      b = bridge(thisa[,c(2,4,4,5)])
      thisfit = BaumWelchT(b$x, b$series.length)
      thisstates = thisfit@ViterbiPath[b$nonvirtuals]
      indstates[a$prom == thisprom] = thisstates
   }

   return(mean(refstates == indstates))

}

refscore = fit_models_separately()

if (file.exists("newscores.rda")) {
   load("newscores.rda")
} else {
   # Resample.
   set.seed(123)
   newscores = rep(NA, 1000)

   for (iter in 1:1000) {
      print (iter)
      a$prom = sample(a$prom)
      newscores[iter] = fit_models_separately()
   }

   save(newscores, file="newscores.rda")
}

# Start writing to file.
sink(file("info_rev2_point5.txt"))

cat("Agreement between global and individual calls:\n")
mean(refscore)

cat("Average agreement when suffling the labels:\n")
mean(newscores)

# End writing to file.
sink(NULL)
