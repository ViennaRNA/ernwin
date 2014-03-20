te <- read.csv('../../results/coarse_grain_ernwin_rmsds.csv', head=F, sep=' ', stringsAsFactors=F)
tr <- read.csv('../../results/coarse_grain_rosetta_rmsds.csv', head=F, sep=' ', stringsAsFactors=F)

colnames(te) <- c('pdb', 'n', 'rmsd')
colnames(tr) <- c('pdb', 'n', 'rmsd')

############################3
# All-atom
############################
te <- read.csv('../../results/all_atom_ernwin_rmsds.csv', head=F, sep=' ', stringsAsFactors=F)
tr <- read.csv('../../results/all_atom_rosetta_rmsds.csv', head=F, sep=' ', stringsAsFactors=F)

te
colnames(tr) <- c('pdb', 'n', 'num_atoms', 'rmsd')
colnames(te) <- c('pdb', 'n', 'num_atoms', 'rmsd')

head(tr)
commonIds <- intersect(te$pdb, tr$pdb)

commonIds

tr1 <- tr[match(commonIds, tr$pdb),]
te1 <- te[match(commonIds, te$pdb),]

mean(tr1$rmsd)
mean(te1$rmsd)

te1
te1$from = "ernwin"
tr1$from = "rosetta"


t1 <- rbind(te1, tr1)
t2 <- transform(t1, pdb=reorder(pdb, num_atoms))
head(t2)
t2
library(ggplot2)

#c + geom_bar(stat="identity", position=position_dodge())
#c + geom_dotplot(binaxis = "y", alpha=0.3) + opts(axis.text.x=theme_text(angle=90, hjust = 1)) + xlab("pdb_file") + ylab("rmsd")
#c + geom_bar(stat="identity", position="dodge") + opts(axis.text.x=theme_text(angle=90, hjust = 1)) + xlab("pdb_file") + ylab("rmsd")
g <- ggplot(data=t2, aes(x=pdb, y=rmsd, color=from, group=from)) + geom_line() + opts(axis.text.x=theme_text(angle=90, hjust = 1, size=15), axis.title.x=theme_text(face="bold", size=15), axis.title.y=theme_text(angle=90, face="bold", size=15)) + xlab("pdb_file") + ylab("rmsd") + geom_point() + ylim(c(0, max(t2$rmsd))) + coord_fixed(ratio=.2)
g
ggsave('rosetta_ernwin_rmsds.png', width=9, height=4)


