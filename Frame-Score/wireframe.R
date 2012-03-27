#Wireframe
#3D
#Aflatoxin cluster
#Aspergillus flavus
library("lattice")

#Read matrix from file and define x,y,z vectors
mat = as.matrix(read.table(file="/home/lupus/Desktop/daten/motifsearch/aflatoxin_aspergillus_flavus/fimo/density/scoring_function/heatmap.csv"))
x = 1:nrow(mat)
y = 1:ncol(mat)
maxZ<-max(mat, na.rm=TRUE)+1
z = 1:maxZ-1

#Draw wireframe of complete contig
wireframe(t(mat), xlab=list("Frame-Nummer",just="left"), ylab=list("Frame-Länge",just="right"), zlab=list("Frame-Score",just="right"), aspect=c(0.5,0.5), drape=TRUE, colorkey=TRUE, at=z, col.regions=rainbow(maxZ, start = 4/6, end = 2/6), zoom=0.65)

#Draw wireframe of small x-range
wireframe(t(mat), xlab=list("Frame-Nummer",just="left"), ylab=list("Frame-Länge",just="right"), zlab=list("Frame-Score",just="right"), aspect=c(0.5,0.5), drape=TRUE, colorkey=TRUE, at=z, col.regions=rainbow(maxZ, start = 4/6, end = 2/6), zoom=0.65,  xlim=c(564,588), scales=list(arrows=FALSE))