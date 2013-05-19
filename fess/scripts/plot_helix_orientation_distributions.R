library(plotrix)
library(ggplot2)
library(MASS)
library(gridExtra)
library(RColorBrewer)
library(grid)

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
      r = diameter / 2
          tt <- seq(0,2*pi,length.out = npoints)
          xx <- center[1] + r * cos(tt)
          yy <- center[2] + r * sin(tt)
          return(data.frame(x = xx, y = yy))
    }

plot_sphere_distributions <- function(t) {
  colnames(t) <- c('type', 'pdb', 's1', 's2', 'u', 'v', 't', 'r', 'u1', 'v1', 'atype', 'something1', 'something2')
  tx <- t[t$type == 'angle',]
  ta <- tx[sample(nrow(tx), 800, replace=T),]
  
  print(str(length(ta$u)))

  ta$x <- sin(ta$u) * cos(ta$v)
  ta$y <- sin(ta$u) * sin(ta$v)
  ta$z <- cos(ta$u)

  taf <- ta[ta$x > 0,]
  tab <- ta[ta$x < 0,]

  circDat <- circleFun(c(0,0), 2, npoints=100)
  circDat$z <- 1

  g <- ggplot(ta, aes(x=y, y=z, colour=x))
  g + geom_point(size=3, alpha=0.4) + scale_colour_gradient(low="red")

  gf <- ggplot(taf, aes(x=y, y=z, colour=x))
  gf <- gf + geom_point(size=5, alpha=0.4) +scale_colour_gradient(low="#DEEBF7", high="#3182BD")
  gf <- gf + geom_path(data=circDat, aes(x,y,colour=z))

  gb <- ggplot(tab, aes(x=y, y=z, colour=-x))
  gb <- gb + geom_point(size=5, alpha=0.4) + scale_colour_gradient(low="#FEE0D2", high="#DE2D26")
  gb <- gb + geom_path(data=circDat, aes(x,y,colour=z))
  g <- arrangeGrob(gf, gb, nrow=1)
  
  return(g)
}


plot_2d_distributions <- function(t) {
  colnames(t) <- c('type', 'pdb', 's1', 's2', 'u', 'v', 't', 'r', 'u1', 'v1', 'atype', 'something1', 'something2')
  tx <- t[t$type == 'angle',]
  ta <- tx
  print(str(length(ta$u)))

  g <- ggplot(ta, aes(x=v, y=u))
  g + geom_point(size=3, alpha=0.4) + scale_colour_gradient(low="red")
}

t <- read.csv('../stats/temp.real.stats', sep=' ', head=F)


gx <- plot_sphere_distributions(read.csv('../stats/temp.real.stats', sep=' ', head=F))
ggsave('helix_distribution_sphere_real.png', gx, width=12, height=6)
plot_2d_distributions(read.csv('../stats/temp.real.stats', sep=' ', head=F))
ggsave('helix_distribution_plane_real.png')

gx <- plot_sphere_distributions(read.csv('../stats/temp.sampled.stats', sep=' ', head=F))
ggsave('helix_distribution_sphere_sampled.png', gx, width=12, height=6)
plot_2d_distributions(read.csv('../stats/temp.sampled.stats', sep=' ', head=F))
ggsave('helix_distribution_plane_sampled.png')

tr <- read.csv('../stats/temp.real.stats', sep=' ', head=F)
ts <- read.csv('../stats/temp.sampled.stats', sep=' ', head=F)

augment_angles <- function(t) {
  colnames(t) <- c('type', 'pdb', 's1', 's2', 'u', 'v', 't', 'r', 'u1', 'v1', 'atype', 'something1', 'something2')
  
  tx <- t[t$type == 'angle',]
  ta <- tx[sample(nrow(tx), 800, replace=T),]
  
  print(str(length(ta$u)))

  ta$x <- sin(ta$u) * cos(ta$v)
  ta$y <- sin(ta$u) * sin(ta$v)
  ta$z <- cos(ta$u)

  return(ta)
}

tra <- augment_angles(tr)
tsa <- augment_angles(ts)


kr <- kde2d(tsa$u, tsa$v, lims=c(0, pi, -pi, pi), n=25)
ks <- kde2d(tra$u, tra$v, lims=c(0, pi, -pi, pi), n=25)

kr$z - ks$z
kd <- kr
kd$z <- log(kr$z) - log(ks$z)

colors <- colorRampPalette(brewer.pal(11, "RdYlBu"))
filled.contour(kd, color.palette=colors)
par(new=T)
plot(tra$u, tra$v, col="red")
par(new=T)
plot(tsa$u, tsa$v, col="blue")

install.packages('rgl')
