t <- read.csv('../stats/temp.stats', sep=' ', head=F)
t <- read.csv('../stats/closure_distances_all.csv', sep=' ', head=F)
head(t)

plot_length <- function(t, l) {
  t1 <- t[t$V1 == l,]
  head(t1)

  x_min = min(t1$V2, t1$V3)
  x_max = max(t1$V2, t1$V3)

  y_min = min(t1$V4)
  y_max = max(t1$V4)

  plot(t1$V3, t1$V4, col='blue', xlim=c(x_min, x_max), ylim=c(y_min, y_max))
  par(new=T)
  plot(t1$V2, t1$V4, col='red', xlim=c(x_min, x_max), ylim=c(y_min, y_max))
}

t1 <- t[t$V4 < 0.15,]
a <- aggregate(t1, by=list(t1$V1), FUN=max)
plot(a$V1, a$V3)
plot(a$V1, a$V2)

cor(a$V1, a$V3)
cor(a$V1, a$V2)

f <- lm(a$V2 ~ a$V1)
abline(f)
1 * 5.94 + 13


max(t[t$V1 == 9 & t$V4 < 0.15,])
