# Import packages
library(sm)

# Read in files
setwd("~/projects/cgtool/build")
atom = read.table("bonds.csv")
cg   = read.table("bonds.csv")

# Check there's the same number of bonds in each file
stopifnot(ncol(atom) == ncol(cg))

# Do t-test and f-test to check mean and variance seperately
for(i in 1:ncol(atom)){
  print(i)
  print(t.test(atom[,i], cg[,i]))
  print(var.test(atom[,i], cg[,i]))
  #sm.density.compare(atom[,i], cg[,i])
  plot(density(atom[,i]))
  plot(density(cg[,i]))
}
