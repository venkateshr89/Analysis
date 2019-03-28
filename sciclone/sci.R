#a=1
#b=1
#c=a+b
#write.table(c,"/home/ubuntu/sciclone/outtmp.txt",quote=F,row.names=F,col.names=F)

library(bmm)
library(sciClone)
library(devtools)

v1 = read.table("/home/ubuntu/sciclone/s1_2caller.dat",header=T)
v2 = read.table("/home/ubuntu/sciclone/s2_2caller.dat",header=T)
regions=read.table("/home/ubuntu/sciclone/exclude.LOH")
head(regions)
cn1=read.table("/home/ubuntu/sciclone/s1/cnv.1718613181_uniq.dat")
cn2=read.table("/home/ubuntu/sciclone/s2/cnv.1718613179_uniq.dat")
names = c("Sample1","Sample2")

sc = sciClone(vafs=list(v1,v2), copyNumberCalls=list(cn1,cn2), sampleNames=names[1:2], regionsToExclude=regions)

#create output
writeClusterTable(sc, "/home/ubuntu/sciclone/results/WES_clusters2")
sc.plot1d(sc,"/home/ubuntu/sciclone/results/WES_clusters2.1d.pdf")
sc.plot2d(sc,"/home/ubuntu/sciclone/results/WES_clusters2.2d.pdf")
