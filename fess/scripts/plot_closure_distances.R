#t <- read.csv('../stats/closure_distances_2.csv', head=F, sep=' ')
t <- read.csv('../stats/all_closure_distances.csv', head=F, sep=' ')
t <- read.csv('../stats/closure_distances_all.csv', head=F, sep=' ')
head(t)
t1 <- t[t$V1 == 10,]

xlim=c(min(t1$V3, t1$V4), max(t1$V3, t1$V4))

plot(t1$V2, t1$V5, col='red', xlim=xlim)
plot(t1$V3, t1$V5, col='red', xlim=xlim)
par(new=T)
plot(t1$V4, t1$V5, col='blue', xlim=xlim)

t2 <- t[t$V5 < .2,]
head(t2)

a1 <- aggregate(t2$V4, by=list(t2$V1), FUN=max)
head(a1)
plot(a1$Group.1, a1$x)
lm(a1$x ~ a1$Group.1)

a2 <- aggregate(t2$V3, by=list(t2$V1), FUN=max)
head(a2)
plot(a2$Group.1, a2$x)
