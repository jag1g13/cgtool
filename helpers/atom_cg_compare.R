# Import packages
library(sm)

# Read in files
setwd("~/projects/cgtool/build")
atom = read.table("aa_bonds.csv")
cg   = read.table("cg_bonds.csv")

# Check there's the same number of bonds in each file
stopifnot(ncol(atom) == ncol(cg))

# Do t-test and f-test to check mean and variance seperately
for(i in 1:ncol(atom)){
  print(i)
  print(t.test(atom[,i], cg[,i]))
  print(var.test(atom[,i], cg[,i]))
  #sm.density.compare(atom[,i], cg[,i])
  plot(density(atom[,i]), col="blue", main=i, xlim=c(0, 0.6))
  lines(density(cg[,i]*1.1), col="red")
}
