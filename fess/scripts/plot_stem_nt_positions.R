t <- read.csv('../stats/stem_nt.stats', sep=' ', head=F)
head(t)
t$dist = sqrt(t$V4 ** 2 + t$V5 ** 2)

t1 <- t[abs(t$V4) < 1,]
t2 <- t[abs(t$V4) < 1 & abs(t$V3) < 1,]

head(t1)
hist(t$dist)
hist(t1$dist)
hist(t2$dist)

print(h)
