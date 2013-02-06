library(ggplot2)

plotorientations <- function(t, ts) {
  length(t[t$V1 < 20,]$V1)
  min_dist = 0
  max_dist = 22

  t1 <- t[t$V1 < max_dist & t$V1 > min_dist,]
  ts1 <- ts[ts$V1 < max_dist & ts$V1 > min_dist,]

  tnew <- data.frame(val=t1$V3, sampled="N")
  tsnew <- data.frame(val=ts1$V3, sampled="Y")
  anew <- rbind(tnew, tsnew)

  length(t1$V3)
  length(ts1$V3)

  ggplot(anew, aes(x=val, fill=sampled)) + geom_density(alpha=.3)
}

ts <- read.csv('../stats/stem_stem_orientations_sampled.csv.old', head=F, sep=' ')
t <- read.csv('../stats/stem_stem_orientations.csv.old', head=F, sep=' ')

png("stem_stem_orientations_old.png")
plotorientations(t, ts)
dev.off()

ts <- read.csv('../stats/stem_stem_orientations_sampled.csv.old', head=F, sep=' ')
t <- read.csv('../stats/stem_stem_orientations.csv.old', head=F, sep=' ')

png("stem_stem_orientations_old.png")
plotorientations(t, ts)
dev.off()


ts <- read.csv('../stats/stem_stem_orientations_sampled.csv', head=F, sep=' ')
t <- read.csv('../stats/stem_stem_orientations.csv', head=F, sep=' ')
png("stem_stem_orientations_new.png")
plotorientations(t, ts)
dev.off()
