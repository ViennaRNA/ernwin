library(ggplot2)

t <- read.table('../stats/cylinder_intersection_fractions.csv', head=F)
ts <- read.table('../stats/cylinder_intersection_fractions_sampled.csv', head=F)
anew <- rbind(data.frame(val=t$V2, sampled="N"), data.frame(val=ts$V2, sampled="Y"))
  
ggplot(anew, aes(x=val, fill=sampled)) + geom_density(alpha=.3)

hist(t$V2)
hist(ts$V2)
