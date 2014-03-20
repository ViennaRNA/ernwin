ta <- read.csv('../stats/temp.stats', head=F, sep=' ')
tb <- read.csv('~/data/ernwin/processed/1jj2/stats/temp.angles', head=F, sep=' ')
head(tb)

ta1 <- ta[ta$V1=='angle',]
tb1 <- tb[tb$V1=='angle',]

head(ta1)

s1=0
s2=9

ta2 <- ta1[ta1$V3 == s1 & ta1$V4 == s2,]
tb2 <- tb1[tb1$V3 == s1 & tb1$V4 == s2,]

library(ggplot2)


ggplot(ta2, aes(x=V10, y=V9, color=V8, fill=V8)) + geom_point(colour="grey50", size=4) + geom_point(aes(colour=V8))

plot(ta2$V9 ~ ta2$V10, ylim=c(0, 3.14), xlim=c(-3.14, 3.14), col='red')
par(new=T)
plot(tb2$V5 ~ tb2$V6, ylim=c(0, 3.14), xlim=c(-3.14, 3.14), col='blue')

foo <- t2$V8
max(foo)
max(t1[t1$V3 == 0 & t1$V4 == 6,]$V8)

hist(foo, prob=T)

lines(density(foo))
lines(density(foo, adjust=2), lty='dotted')
curve(dnorm(x, mean=mean(foo), sd=sd(foo)), add=T)

#### Lengths

t2 <- t1[t1$V3 == 0 & t1$V4 == 5,]
hist(t2$V8, n=8)


qqnorm(t2$V8)

t1[t1$V2 == '1gid',]



library(sn)
t2 <- t1[t1$V3 == 0 & t1$V4 == 5,]
r = (t2$V8)
png('junction_closure_distribution.png')
sa <- sn.mle(y=r)
dev.off()
sa
log(dsn(5, sa$cp[1], sa$cp[2], sa$cp[3]))

num <- 10
log(dsn(num, 13.35, 5.76, 3.34)) - log(dsn(num, 16.9, 13.13, 0.76))


### Long Range Interaction Lengths
s <- read.csv('../stats/temp.longrange.contact', head=F, sep=' ')

head(s)

r = s$V3
png('skew_normal_native_interactions.png')
sa <- sn.mle(y=r)
dev.off()

log(dsn(115, sa$cp[1], sa$cp[2], sa$cp[3]))

### Plot all interior loop angles
library(ggplot2)

t <- read.csv('../stats/temp.stats', head=F, sep=' ')
t <- read.csv('~/data/ernwin/processed/1jj2/stats/temp.angles', head=F, sep=' ')

t1 <- t[t$V1=='angle',]
t1 <- t1[t1$V3 > 0,]
t1 <- t1[t1$V11 == 0,]
t2 <- t1[t1$V3 > 0,]
head(t1)




png("rosetta_interior_loop_stem_angles.png")
c <- ggplot(data=t2)
c + geom_point(aes(x=V6, y=V5, size=V4/V3), alpha=0.5) + xlab("azimuth") + ylab("altitude")
dev.off()

plot(t2$V6, t2$V5)
