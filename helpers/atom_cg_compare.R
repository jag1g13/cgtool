# Import packages
library(sm)

# Read in files
atom = read.table("../aa_bonds.dat")
cg   = read.table("cg_bonds.dat")

# Check there's the same number of bonds in each file
stopifnot(ncol(atom) == ncol(cg))

# Do t-test and f-test to check mean and variance separately
par(xaxs='i',yaxs='i')
for(i in 1:ncol(atom)){
  #print(i)
  #print(t.test(atom[,i], cg[,i]))
  #print(var.test(atom[,i], cg[,i]))
  #sm.density.compare(atom[,i], cg[,i])
  plot(density(atom[,i]), col="blue", main=i, xlim=c(0,1), xlab="Bond length (nm)")
  lines(density(cg[,i]), col="red", main=i, xlim=c(0,1))
  print(mean(atom[,i]))
  print(mean(cg[,i]))
  print(sd(atom[,i]))
  print(sd(cg[,i]))
}

atom = read.table("../aa_angles.dat")
cg   = read.table("cg_angles.dat")
stopifnot(ncol(atom) == ncol(cg))
par(xaxs='i',yaxs='i')
for(i in 1:ncol(atom)){
  plot(density(atom[,i]), col="blue", main=i, xlim=c(-180,180), xlab="Bond angle (deg)")
  lines(density(cg[,i]), col="red", main=i, xlim=c(-180,180))
  print(mean(atom[,i]))
  print(mean(cg[,i]))
  print(sd(atom[,i]))
  print(sd(cg[,i]))
}

atom = read.table("../aa_dihedrals.dat")
cg   = read.table("cg_dihedrals.dat")
stopifnot(ncol(atom) == ncol(cg))
par(xaxs='i',yaxs='i')
for(i in 1:ncol(atom)){
  plot(density(atom[,i]), col="blue", main=i, xlim=c(-180,180), xlab="Dihedral angle (deg)")
  lines(density(cg[,i]), col="red", main=i, xlim=c(-180,180))
  print(mean(atom[,i]))
  print(mean(cg[,i]))
  print(sd(atom[,i]))
  print(sd(cg[,i]))
}

