t <- read.csv('../stats/stem_nt.stats', sep=' ', head=F)
t1 <- read.csv('../stats/stem_nt_sampled.stats', sep=' ', head=F)

head(t)

t$dist = sqrt(t$V4 ** 2 + t$V5 ** 2)
t1$dist = sqrt(t1$V4 ** 2 + t1$V5 ** 2)

t$ang = atan2(t$V4, t$V5)
t1$ang = atan2(t1$V4, t1$V)


require(CircStats)
library(circular)

d <- circular(t$ang)
radial.plot(d, t$ang)

t$ang

b <- head(t$ang, n=100)


b <- circular(b)

a <- density(b, bw=200)

plot(b,stack=T,cex=0, ticks=T)
lines(a)

h <- hist(t$ang)

radial.pie(h$counts, h$mids)
radial.plot(h$counts, h$mids, rp.type="r", lwd=10)


t1 <- t[abs(t$V4) < 1,]
t2 <- t[abs(t$V4) < 1 & abs(t$V3) < 1,]

head(t1)
hist(t$dist)
hist(t1$dist)
hist(t2$dist)

print(h)

library(lattice)
wireframe(t$V3 ~ t$V4 * t$V5)

data = data.frame(
    x = rep( c(0.1, 0.2, 0.3, 0.4, 0.5), each=5),
    y = rep( c(1, 2, 3, 4, 5), 5)
)

data$z = runif(
    25,
    min = (data$x*data$y - 0.1 * (data$x*data$y)),
    max = (data$x*data$y + 0.1 * (data$x*data$y))
)

data
str(data)
wireframe(z ~ x * y, data=data)

t <- read.csv('../stats/vres_distances.csv', sep=' ', head=F)
min(t$V1)
?rose.diag
