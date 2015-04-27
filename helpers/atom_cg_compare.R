# Import packages
library(sm)

# Read in files
setwd("~")
atom = read.table("aa_bonds.csv")
cg   = read.table("cg_bonds.csv")

# Check there's the same number of bonds in each file
stopifnot(ncol(atom) == ncol(cg))

# Do t-test and f-test to check mean and variance separately
for(i in 1:ncol(atom)){
  print(i)
  print(t.test(atom[,i+3], cg[,i]))
  print(var.test(atom[,i+3], cg[,i]))
  #sm.density.compare(atom[,i], cg[,i])
  plot(density(atom[,i]), col="blue", main=i, xlim=c(0, 1))
  lines(density(cg[,i]), col="red", main=i, xlim=c(0, 1))
}
