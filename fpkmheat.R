##########################
# Creates heatmaps of differentially expressed genes from 
# fpkm.txt file which are generated through getDiffExpression.pl on raw.txt
# which outputted diffbatch.txt and diffnobatch.txt

# diffnobatch is about 448 genes
# separated into upregulated genes upregulatedgenes.txt
# and downregulatedgenes.txt
#########################


fpkmyeast <- read.delim("~/Documents/salk/fpkmyeast.txt")
diffbatch <- read.delim("~/Documents/salk/diffbatch.txt")
diffnobatch <- read.delim("~/Documents/salk/diffnobatch.txt")
dev.off()

library(gplots)

#Change rownames of fpkm to standard gene symbol
names(fpkmyeast)[1]<-"Transcript"
fpkmname <- as.character(fpkmyeast$Annotation.Divergence)
fpkmsym<- sapply(fpkmname, function(x) strsplit(x,"|", fixed = T)[[1]][1])
rownames(fpkmyeast) <- fpkmsym

#transform original fpkm numeric values only for heatmap
transformed <- log(fpkmyeast[,9:14]+5) #genes by sample

################# DIFF BATCH ###############

#Change rownames of diffbatch to standard gene symbol
names(diffbatch)[1]<-"Transcript"
batchnames <- as.character(diffbatch$Annotation.Divergence)
batchsym <- sapply(batchnames, function(x) strsplit(x,"|", fixed = T)[[1]][1])
rownames(diffbatch) <- batchsym

#diffbatch fold change and adj p-value find genes in diffbatch
batchsig <-subset(diffbatch, abs(diffbatch$WT.vs..abf2.Log2.Fold.Change) >= 0.5
                  & diffbatch$WT.vs..abf2.adj..p.value <= 0.05)

#Write to table for ORA
batchsigtrans <- data.frame(batchsig$Transcript)
write.table(batchsigtrans,"batchsigtrans.txt",sep="\t",row.names=FALSE, col.names = FALSE, quote=FALSE)

#Find genes in fpkmtransformed and subset if found in diffbatch
fpkmbatchsig <- subset(transformed, rownames(transformed) %in% rownames(batchsig))


#Rename samplenames for map
comb2batch <- sapply(colnames(fpkmbatchsig), function(x) paste(strsplit(x,"_", fixed = T)[[1]][2],
                                                                  strsplit(x,"_",fixed = T)[[1]][3],
                                                                  sep = "_"
))
colnames(fpkmbatchsig) <- comb2batch
#heatmap.2(as.matrix(fpkmbatchsignum), col=topo.colors(100))
yb<-colorRampPalette(c("blue","white","red"))(100)




############# LARGE HEATMAP BATCH GENES ################
pdf("Batch.pdf", width=8, height=30)
par(mar=c(2,2,2,2), cex=1.0)
heatmap.2(as.matrix(fpkmbatchsig), col=yb, scale="row",
                    key=TRUE, symkey=FALSE, density.info="none", trace="none",
                    cexRow=0.3,lhei=c(0.05,0.95), lwid = c(0.5,0.5),
                    reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
                    distfun=function(x) as.dist(1-cor(t(x))),
                    hclustfun=function(x) hclust(x, method="ward.D2"))
dev.off()


############# SMALL HEATMAP BATCH GENES ################
pdf("BatchSmall.pdf", width=7, height=7)
par(mar=c(2,2,2,2), cex=1.0)
heatmap.2(as.matrix(fpkmbatchsig), col=yb, scale="row", labRow = "",
          key=TRUE, symkey=FALSE, density.info="none", trace="none",
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x, method="ward.D2"))
dev.off()


################### NOBATCH #############


#Change rownames of diffbatch to standard gene symbol
names(diffnobatch)[1]<-"Transcript"
nobatchnames <- as.character(diffnobatch$Annotation.Divergence)
nobatchsym <- sapply(nobatchnames, function(x) strsplit(x,"|", fixed = T)[[1]][1])
rownames(diffnobatch) <- nobatchsym

#diffbatch fold change and adj p-value find genes in diffbatch
nobatchsig <-subset(diffnobatch, abs(diffnobatch$WT.vs..abf2.Log2.Fold.Change) >= 0.5
                  & diffnobatch$WT.vs..abf2.adj..p.value <= 0.05)

#Find genes in fpkmtransformed and subset if found in diffbatch
fpkmnobatchsig <- subset(transformed, rownames(transformed) %in% rownames(nobatchsig))

#Write to table for ORA
nobatchsigtrans <- data.frame(nobatchsig$Transcript)
write.table(nobatchsigtrans,"nobatchsigtrans.txt",sep="\t",row.names=FALSE, col.names = FALSE, quote=FALSE)


#Rename samplenames for map
comb2nobatch <- sapply(colnames(fpkmnobatchsig), function(x) paste(strsplit(x,"_", fixed = T)[[1]][2],
                                                               strsplit(x,"_",fixed = T)[[1]][3],
                                                               sep = "_"
))
colnames(fpkmnobatchsig) <- comb2nobatch
write.table(diffnobatch,"allgenes.txt",sep="\t",row.names=TRUE, col.names = TRUE, quote=FALSE)

############# LARGE HEATMAP NOBATCH GENES ################
pdf("NoBatchHeat.pdf", width=8, height=30)
par(mar=c(2,2,2,2), cex=1.0)
heatmap.2(as.matrix(fpkmnobatchsig), col=yb, scale="row",
          key=TRUE, symkey=FALSE, density.info="none", trace="none",
          cexRow=0.3,lhei=c(0.05,0.95), lwid = c(0.5,0.5),
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x, method="ward.D2"))
dev.off()

############# SMALL HEATMAP NOBATCH GENES ################
pdf("NoBatchHeatSmall.pdf", width=7, height=7)
par(mar=c(2,2,2,2), cex=1.0)
heatmap.2(as.matrix(fpkmnobatchsig), col=yb, scale="row", labRow = "",
          key=TRUE, symkey=FALSE, density.info="none", trace="none",
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x, method="ward.D2"))
dev.off()


########## CHECK WHICH GENES ARE IN BATCH BUT NOT IN NOBATCH #####
diffgenes <- setdiff(rownames(fpkmbatchsig), rownames(fpkmnobatchsig))
write.table(diffgenes,"diffbatchnobatch.txt",sep="\t",row.names=FALSE, col.names = FALSE, quote=FALSE)


#### getDiffexpression outputs diffbatch and diffnobatch.txt which has pvalue fold change

############  #########################





############# Upregulated GENES #############

#diffbatch fold change and adj p-value find genes in diffbatch

#161 0.5-FC 22 1-FC
upnobatch <-subset(diffnobatch, diffnobatch$WT.vs..abf2.Log2.Fold.Change >= 0.5
                  & diffnobatch$WT.vs..abf2.adj..p.value <= 0.05)
write.table(upnobatch,"upregulatedgenes.txt",sep="\t",row.names=TRUE, col.names = TRUE, quote=FALSE)



#### down regulated genes######
#287 0.5-FC 50 1-FC
downnobatch <-subset(diffnobatch, diffnobatch$WT.vs..abf2.Log2.Fold.Change <= -0.5
                   & diffnobatch$WT.vs..abf2.adj..p.value <= 0.05)
write.table(downnobatch,"downregulatedgenes.txt",sep="\t",row.names=TRUE, col.names = TRUE, quote=FALSE)


#Write to table for ORA
batchsigtrans <- data.frame(batchsig$Transcript)
write.table(batchsigtrans,"batchsigtrans.txt",sep="\t",row.names=TRUE, col.names = TRUE, quote=FALSE)

#Find genes in fpkmtransformed and subset if found in diffbatch
fpkmbatchsig <- subset(transformed, rownames(transformed) %in% rownames(batchsig))



