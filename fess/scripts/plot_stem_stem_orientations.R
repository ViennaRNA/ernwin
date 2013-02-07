library(ggplot2)
library(gridExtra)


plotorientations <- function(t, ts) {
  length(t[t$V1 < 20,]$V1)
  min_dist = 0
  max_dist = 30

  t1 <- t[t$V1 < max_dist & t$V1 > min_dist,]
  ts1 <- ts[ts$V1 < max_dist & ts$V1 > min_dist,]

  min(t1$V3, abs(pi - t1$V3))

  t1$offset3 <- sapply(t1$V3,function(x) min(x, abs(pi - x)))
  ts1$offset3 <- sapply(ts1$V3,function(x) min(x, abs(pi - x)))

  t1$offset2 <- sapply(t1$V2,function(x) min(x, abs(pi - x)))
  ts1$offset2 <- sapply(ts1$V2,function(x) min(x, abs(pi - x)))

head(t1)
  
  length(t1$V1)
  length(ts1$V1)
  #sample_len <- max(length(t1$V3), length(ts1$V3))
  sample_len <- 10000
  
  t1_V3 <- sample(t1$offset3, sample_len, replace=T)
  ts1_V3 <- sample(ts1$offset3, sample_len, replace=T)

  t1_V2 <- sample(t1$offset2, sample_len, replace=T)
  ts1_V2 <- sample(ts1$offset2, sample_len, replace=T)

  anew <- rbind(data.frame(val=t1_V3, sampled="N"),
                data.frame(val=ts1_V3, sampled="Y"))

  anew1 <- rbind(data.frame(val=t1_V2, sampled="N"),
                 data.frame(val=ts1_V2, sampled="Y"))

  length(t1_V3)
  length(ts1_V3)

  #plot(t1$V2, t1$V3)
  #plot(ts1$V2, ts1$V3)

  g1 <- ggplot(anew, aes(x=val, fill=sampled)) + geom_density(alpha=.3)
  g2 <- ggplot(anew1, aes(x=val, fill=sampled)) + geom_density(alpha=.3)

  grid.arrange(g1,g2)
}
plotorientations(t, ts)
ts <- read.csv('../stats/stem_stem_orientations_sampled.csv', head=F, sep=' ')
t <- read.csv('../stats/stem_stem_orientations.csv', head=F, sep=' ')

png("stem_stem_orientations_new.png")
plotorientations(t, ts)
dev.off()

ts <- read.csv('../stats/stem_stem_orientations_sampled.csv.old', head=F, sep=' ')
t <- read.csv('../stats/stem_stem_orientations.csv.old', head=F, sep=' ')
png("stem_stem_orientations_old.png")
plotorientations(t, ts)
dev.off()

ts <- read.csv('../stats/stem_stem_orientations_sampled.csv.old.1', head=F, sep=' ')
t <- read.csv('../stats/stem_stem_orientations.csv.old.1', head=F, sep=' ')
png("stem_stem_orientations_old_1.png")
plotorientations(t, ts)
dev.off()



