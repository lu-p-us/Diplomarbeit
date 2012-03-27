#Heatmap
#Aflatoxin cluster
#Aspergillus flavus
library("fields")

#Read matrix from file and define x,y vectors
mat = as.matrix(read.table(file="/home/lupus/Desktop/daten/motifsearch/aflatoxin_aspergillus_flavus/fimo/density/scoring_function/heatmap.csv"))
x = 1:nrow(mat)
y = 1:ncol(mat)
maxZ<-max(mat, na.rm=TRUE)+1

#Draw heatmap of complete contig
image.plot(y, x, t(mat), xlab="Frame-Nummer", ylab="Frame-Länge", col=rainbow(maxZ, start = 4/6, end = 2/6), nlevel=maxZ, legend.shrink=1, legend.lab="Frame-Score")

# Draw heatmap of small x-range
image.plot(y, x, t(mat), xlab="Frame-Nummer", ylab="Frame-Länge", xlim=c(565,585), ylim=c(1,25), col=rainbow(maxZ, start = 4/6, end = 2/6), nlevel=maxZ, legend.shrink=1, legend.lab="Frame-Score")