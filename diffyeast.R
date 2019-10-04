diffbatch <- read.delim("~/Documents/salk/diffbatch.txt")
diffnobatch <- read.delim("~/Documents/salk/diffnobatch.txt")

#This is the differential expression of the batch and nobatch files from 
#getDiffExpression, not the right valuessee fpkmheat

library(gplots)

#diffbatch fold change and adj p-value
names(diffbatch)[1]<-"Transcript"
batchsig <-subset(diffbatch, abs(diffbatch$WT.vs..abf2.Log2.Fold.Change) >= 0.5
                  & diffbatch$WT.vs..abf2.adj..p.value <= 0.05)
batchname <- data.frame(batchsig$Transcript)

#Get numeric values of batch 
batchsignum <- batchsig[,9:14]
rownames(batchsignum) <- batchsig$Transcript
#Rename samplenames
comb2batch <- sapply(colnames(batchsignum), function(x) paste(strsplit(x,"_", fixed = T)[[1]][2],
                                                         strsplit(x,"_",fixed = T)[[1]][3],
                                                         sep = "_"
                                                          ))
colnames(batchsignum) <- comb2batch
#heatmap.2(as.matrix(batchsignum), col=topo.colors(100))

#heatmap.2(exprs(esetSel), col=topo.colors(75), scale="none", ColSideColors=patientcolors,
#          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
heatmap.2(as.matrix(batchsignum), col=redgreen, labRow = "", scale="row",
          key=TRUE, symkey=FALSE, density.info="none", trace="none",
          cexRow=0.5)





#diffnobatch fold change and adj p-value
names(diffnobatch)[1]<-"Transcript"
nobatchsig <-subset(diffnobatch, abs(diffnobatch$WT.vs..abf2.Log2.Fold.Change) >= 0.5
                  & diffnobatch$WT.vs..abf2.adj..p.value <= 0.05)
nobatchname <- data.frame(nobatchsig$Transcript)

#Get numeric values of batch 
batchsignum <- batchsig[,9:14]
rownames(batchsignum) <- batchsig$Transcript
#Rename samplenames
comb2 <- sapply(colnames(batchsignum), function(x) paste(strsplit(x,"_", fixed = T)[[1]][2],
                                                         strsplit(x,"_",fixed = T)[[1]][3],
                                                         sep = "_"
))
colnames(batchsignum) <- comb2
heatmap.2(as.matrix(batchsignum), col=topo.colors(100))


