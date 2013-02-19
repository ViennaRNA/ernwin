te <- read.csv('../../results/coarse_grain_ernwin_rmsds.csv', head=F, sep=' ', stringsAsFactors=F)
tr <- read.csv('../../results/coarse_grain_rosetta_rmsds.csv', head=F, sep=' ', stringsAsFactors=F)

te
tr

commonIds <- intersect(te$V1, tr$V1)

tr1 <- tr[match(commonIds, tr$V1),]
te1 <- te[match(commonIds, te$V1),]

mean(tr1$V3)
mean(te1$V3)

te1$from = "ernwin"
tr1$from = "rosetta"


t1 <- rbind(te1, tr1)
t2 <- transform(t1, V1=reorder(V1, V3))
head(t2)

mean(te1$V3)
mean(tr1$V3)

library(ggplot2)
head(t2)
png("rosetta_ernwin_rmsds.png")
c <- ggplot(data=t2, aes(x=V1,y=V3,fill=from))
#c + geom_bar(stat="identity", position=position_dodge())
#c + geom_dotplot(binaxis = "y", alpha=0.3) + opts(axis.text.x=theme_text(angle=90, hjust = 1)) + xlab("pdb_file") + ylab("rmsd")
#c + geom_bar(stat="identity", position="dodge") + opts(axis.text.x=theme_text(angle=90, hjust = 1)) + xlab("pdb_file") + ylab("rmsd")
ggplot(data=t2, aes(x=V1, y=V3, color=from, group=from)) + geom_line() + opts(axis.text.x=theme_text(angle=90, hjust = 1)) + xlab("pdb_file") + ylab("rmsd") + geom_point() + ylim(c(0, max(t2$V3)))
dev.off()

