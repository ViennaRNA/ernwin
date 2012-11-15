library(ggplot2)

t <- read.csv('../stats/stem_stem_orientations.csv', head=F, sep=' ')
ts <- read.csv('../stats/stem_stem_orientations_sampled.csv', head=F, sep=' ')

t1 <- t[t$V1 < 20,]
t2 <- t[t$V2 < 20,]

tnew <- data.frame(val=t1$V3, sampled="N")
tsnew <- data.frame(val=ts1$V3, sampled="Y")
anew <- rbind(tnew, tsnew)

ggplot(anew, aes(x=val, fill=sampled)) + geom_density(alpha=.3)

