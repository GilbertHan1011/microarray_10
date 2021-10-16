#== setting environment-----------------
options("BioC_mirror" = "http://mirrors.tuna.tsinghua.edu.cn/bioconductor/") 
options(repos=c(CRAN="http://mirror.tuna.tsinghua.edu.cn/CRAN/")) 
library(oligo)
library(ggplot2)
library(affy)
library(limma)
celpath = "data/"
CELFiles <- list.files(celpath,full.names=TRUE, pattern = "CEL")
rawData <- read.celfiles(CELFiles)

#== quality control-----------------------
ph = rawData@phenoData
ph@data[,1]<- c("control1", "control2", "control3","control4", "treatment1","treatment2","treatment3","treatment4")
#histogram plot
for (i in 1:8)
{
  name = paste("histogram",i,".jpg",sep="")
  jpeg(name)
  hist(rawData[,i],lwd=2,which='pm',ylab='Density',xlab='Log2 intensities',main=ph@data$index[i],target="core")
  dev.off()
}
color=c('green','green','green',"green","red","red",'red','red')
hist(rawData[,1:8],lwd=2,which='pm',col=color,ylab='Density',xlab='Log2 intensities',main='Histogram of raw data',target="core")

# box plot
pmexp = pm(rawData)
sampleNames = vector()
logs = vector()
for (i in 1:8)
{
  sampleNames = c(sampleNames,rep(ph@data[i,1],dim(pmexp)[1]))
  logs = c(logs,log2(pmexp[,i]))
}
logData = data.frame(logInt=logs,sampleName=sampleNames)
dataBox = ggplot(logData,aes(sampleName,logInt))
dataBox + geom_boxplot()



#== normalization-----------------------------
data.rma = oligo::rma(rawData)
data.matrix = exprs(data.rma)
color=c('green','green','green',"green",'red','red','red','red')
data.PC = prcomp(t(data.matrix),scale.=TRUE)
plot(data.PC$x[,1:2],col=color)

#== DE analysis-----------------------------------------
ph@data[ ,2] = c("control","control","control","control","mutant","mutant","mutant","mutant")
colnames(ph@data)[2]="source"
groups = ph@data$source
f = factor(groups,levels=c("control","mutant"))
design = model.matrix(~ 0 + f)
colnames(design) = c("control","mutant")
data.fit = lmFit(data.matrix,design)

contrast.matrix = makeContrasts(mutant-control,levels=design)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
data.fit.eb = eBayes(data.fit.con)
options(digits=2)
tab = topTable(data.fit.eb,coef=1,number = 2000,adjust.method="BH")

#== gene ID transform-------------------------
library(mogene20sttranscriptcluster.db)
library("AnnotationDbi")
probeId <- rownames(data.matrix)
OUT <- select(mogene20sttranscriptcluster.db,keys= probeId, columns=c("SYMBOL", "ENTREZID", "GENENAME"),keytype="PROBEID")
dim(na.omit(OUT))
dim(OUT)


#== 