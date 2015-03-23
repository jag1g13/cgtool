setwd("~/projects/cgtool/build")
tab = read.table("bonds.csv")

for(i in 1:ncol(tab)){
  print(i)
  tryCatch({
    print(shapiro.test(tab[,i]))
  }, error=function(e){
    print("Too much data for Shapiro")
  })
  hist(tab[,i], main=i, n=15)
  qqnorm(tab[,i], main=i)
}

