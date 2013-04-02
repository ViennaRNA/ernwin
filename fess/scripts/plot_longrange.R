options(width=150)
t <- read.table('../stats/temp.longrange.stats', head=F)
ts <- read.table('../stats/temp.longrange.stats.sampled', head=F)

head(t)
head(ts)
colnames(t) <- c('key1', 'type1', 'len1', 'key2', 'type2', 'len2', 'dist', 'seq1', 'seq2', 'longrange')
colnames(ts) <- c('key1', 'type1', 'len1', 'key2', 'type2', 'len2', 'dist', 'seq1', 'seq2', 'longrange')

t1 <- t[t$dist > 0 & t$dist < 150,]
ts1 <- ts[ts$dist > 0 & ts$dist < 150,]
#t1 <- t[t$dist > 0,]
#ts1 <- ts[ts$dist > 0,]

loop_loop <- t1[t1$type1 == 'l' & t1$type2 == 'l',]
ts_loop_loop <- ts1[ts1$type1 == 'l' & ts1$type2 == 'l',]

g1 <- ggplot(loop_loop, aes(x=dist, fill=longrange)) + geom_density(alpha=0.3)
g2 <- ggplot(ts_loop_loop, aes(x=dist, fill=longrange)) + geom_density(alpha=0.3)

grid.arrange(g1, g2)


loop <- t1[t1$type1 == 'l',]
loop_stem <- t1[t1$type1 == 'l' & t1$type2 == 's',]
loop_bulge <- t1[t1$type1 == 'l' & t1$type2 == 'i',]


ls_y <- loop_stem[loop_stem$longrange == "Y",]
lb_y <- loop_bulge[loop_bulge$longrange == "Y",]

hist(ls_y$len1)
hist(lb_y$len1)

hist(ls_y[ls_y$len1 == 5,]$dist)
ls_y
lb_y

library(ggplot2)
library(gridExtra)
                

a <- loop_stem
a <- loop_bulge
a <- loop
a <- loop_loop
ay <- a[a$longrange=="Y",]
an <- a[a$longrange=="N",]

hist(ay$len1)
hist(an$len1)
ggplot(a, aes(x=len1, fill=longrange)) + geom_density(alpha=0.3)
ggplot(a, aes(x=dist, fill=longrange)) + geom_density(alpha=0.3)



interacting <- unique(ay$key1)
noninteracting <- setdiff(unique(an$key1), unique(ay$key1))

hist(a[match(noninteracting, a$key1),]$len1)
a[match(interacting, a$key1),]

dfy <- data.frame(len1 = a[match(interacting, a$key1),]$len1, longrange="Y")
dfn <- data.frame(len1 = a[match(noninteracting, a$key1),]$len1, longrange="N")
df <- rbind(dfy, dfn)
df
ggplot(df, aes(x=len1, fill=longrange, y=..count../sum(..count..))) + geom_histogram(position="dodge")
