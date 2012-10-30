t <- read.csv('../stats/temp.stats', sep=' ', head=F)
t <- t[t$V1 == 'angle',]
t <- t[t$V3 == 0,]
head(t)



max(t[t$V4 == 3,]$V8) - max(t[t$V4 == 2,]$V8)
max(t[t$V4 == 4,]$V8) - max(t[t$V4 == 3,]$V8)
max(t[t$V4 == 5,]$V8) - max(t[t$V4 == 4,]$V8)
max(t[t$V4 == 6,]$V8) - max(t[t$V4 == 5,]$V8)

max(t[t$V4 == 2,]$V8)
max(t[t$V4 == 3,]$V8)
max(t[t$V4 == 4,]$V8)
max(t[t$V4 == 5,]$V8)
max(t[t$V4 == 6,]$V8)
