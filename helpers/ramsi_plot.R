# Plot membrane thickness
dat = read.table("thickness_avg.dat")
names(dat) = seq(1,200)
library(ggplot2)
library(reshape2)
m = as.matrix(dat)
mel = melt(m)
gg = ggplot(mel, aes(Var1, Var2, fill=value))
gg = gg + theme(axis.title.x = element_text(size=30))
gg = gg + labs(x=expression(phi), y=expression(psi))
gg = gg + geom_raster()
gg = gg + coord_equal()
gg = gg + scale_fill_gradient(low="red", high="yellow")
gg = gg + scale_x_continuous(expand = c(0, 0))
gg = gg + scale_y_continuous(expand = c(0, 0))
gg = gg + theme_bw()
gg

# Plot area per lipid over time
apl = read.table("APL.dat")
plot(apl$V2~apl$V1, ylim=c(0,2), type='n')
cols = c("black", "blue", "red", "yellow")
for(i in 2:ncol(apl)){
  lines(apl[,i] ~ apl[,1], col=cols[i])
}