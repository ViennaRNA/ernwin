library(ggplot2)

t <- read.table('../stats/cylinder_intersection_fractions.csv', head=F)
ts <- read.table('../stats/cylinder_intersection_fractions_sampled.csv', head=F)
ts1 <- read.table('../stats/cylinder_intersection_fractions_1mfq.csv', head=F)
anew <- rbind(data.frame(val=t$V2, sampled="N"), data.frame(val=ts$V2, sampled="Y"), data.frame(val=ts1$V2, sampled='X'))
anew <- rbind(data.frame(val=t$V2, sampled="N"), data.frame(val=ts$V2, sampled="Y"))
  
ggplot(anew, aes(x=val, fill=sampled)) + geom_density(alpha=.3)

hist(ts1$V2)
hist(ts$V2)
hist(t$V2)

