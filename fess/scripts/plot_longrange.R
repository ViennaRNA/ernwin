library(ggplot2)
library(gridExtra)

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

# Distance given interaction
loop_loop_y <- loop_loop[loop_loop$longrange == 'Y',]
ts_loop_loop_y <- ts_loop_loop[ts_loop_loop$longrange == 'Y',]

g1 <- ggplot(loop_loop, aes(x=dist, fill=longrange)) + geom_density(alpha=0.3)
g2 <- ggplot(ts_loop_loop, aes(x=dist, fill=longrange)) + geom_density(alpha=0.3)

grid.arrange(g1, g2)

# Interaction given size
head(loop_loop)
test <- 'b263'
loop_loop[loop_loop$key1 == test,]

ll_score <- function(loop1, loop2) {
  size_density_y <- density(loop_loop_y$len1)
  size_density_y$y
  
  size_density_n <- density(loop_loop_n$len1)

  dist_density_y <- density(loop_loop_y$dist)
  dist_density_n <- density(loop_loop_n$dist)

  pl_size_l1 <- 
}

ggplot(loop_loop, aes(x=len1, fill=longrange)) + geom_density(alpha=0.3)

## Difference of distance densities
d <- density(loop_loop$dist)
ds <- density(ts_loop_loop$dist)

a <- seq(0, max(d$x), .5)
da <- approx(d$x, d$y, a)
dsa <- approx(ds$x, ds$y, a)

plot(a, log(da$y) - log(dsa$y))

## Shortest loop loop distances
head(loop_loop)
a <- aggregate(dist ~ key1, data=loop_loop, FUN=min)
b <- merge(a, loop_loop)

a_ts <- aggregate(dist ~ key1, data=ts_loop_loop, FUN=min)
b_ts <- merge(a_ts,  ts_loop_loop)

g1 <- ggplot(b, aes(x=dist, fill=longrange)) + geom_density(alpha=0.3)
g2 <- ggplot(b_ts, aes(x=dist, fill=longrange)) + geom_density(alpha=0.3)

grid.arrange(g1, g2)

# All loop distances
all_loops <- rbind(data.frame(dist=loop_loop$dist, sampled='N'),
                   data.frame(dist=ts_loop_loop$dist, sampled='Y'))

head(all_loops)
ggplot(all_loops, aes(x=dist, fill=sampled)) + geom_density(alpha=0.3)


## Loop stem interactions
loop_stem <- t1[t1$type1 == 'l' & (t1$type2 == 's' | t1$type2 == 'i'),]
ls_y <- loop_stem[loop_stem$longrange == "Y",]
hist(ls_y$len1)
hist(ls_y$dist)

## Loop junction interactions
loop_junction <- t1[t1$type1 == 'l' & t1$type2 == 'm',]
loop_junction_y <- loop_junction[loop_junction$longrange == 'Y',]

hist(loop_junction_y$len2)

loop_bulge <- t1[t1$type1 == 'l' & t1$type2 == 'i',]
lb_y <- loop_bulge[loop_bulge$longrange == "Y",]
hist(lb_y$len1)

hist(ls_y$len1)
hist(lb_y$len1)

# All loop junction
all_juncs = rbind

head(b_ts)

hist(ls_y[ls_y$len1 == 5,]$dist)
ls_y
lb_y

library(ggplot2)
library(gridExtra)
                
loop <- t1[t1$type1 == 'l',]
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
