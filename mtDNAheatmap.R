library(gplots)
yb<-colorRampPalette(c("blue","white","red"))(100)

fpkmyeast <- read.delim("~/Documents/salk/fpkmyeast.txt")


#Change rownames of fpkm to standard gene symbol
names(fpkmyeast)[1]<-"Transcript"
fpkmname <- as.character(fpkmyeast$Annotation.Divergence)
fpkmsym<- sapply(fpkmname, function(x) strsplit(x,"|", fixed = T)[[1]][1])
rownames(fpkmyeast) <- fpkmsym

#Read in curated and high throughput
curatedgenes <- read.table("~/Downloads/mitochondrion_annotations_Manually Curated.txt", sep = "\t", skip = 1, header = T, comment.char = "!")
highthrugenes <- read.table("~/Downloads/mitochondrion_annotations_High-throughput.txt", sep = "\t", skip = 1, header = T, comment.char = "!")


#Write curated genes for IPA pathway analysis
transcur <- as.character(curatedgenes$Gene.Systematic.Name)
curinfpkm <- subset(fpkmyeast, fpkmyeast$Transcript %in% transcur)
curtransf <- log(curinfpkm[,9:14]+5) #genes by sample

#Remove genes that haven't changed between WT and KO
curvar <- curtransf[ apply(curtransf, 1, var, na.rm = TRUE) != 0 , ]

#Renamed curated col names for heatmap
comb2cur <- sapply(colnames(curvar), function(x) paste(strsplit(x,"_", fixed = T)[[1]][2],
                                                            strsplit(x,"_",fixed = T)[[1]][3],
                                                            sep = "_"
))
colnames(curvar) <- comb2cur


#############Curated heatmap ##############
pdf("Curated.pdf", width=7, height=7)
par(mar=c(2,2,2,2), cex=1.0)
heatmap.2(as.matrix(curvar), col=yb, scale="row", dendrogram = "row", labRow = "",
          key=TRUE, symkey=FALSE, density.info="none", trace="none",
          cexRow=0.3,lhei=c(0.20,0.70), lwid = c(0.25,0.5),
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x, method="ward.D2"))

dev.off()

#####High throughput genes
transhighthru <- as.character(highthrugenes$Gene.Systematic.Name)
highinfpkm <- subset(fpkmyeast, fpkmyeast$Transcript %in% transhighthru)
hightransf <- log(highinfpkm[,9:14]+5) #genes by sample


#Remove genes that haven't changed between WT and KO
highvar <- hightransf[ apply(hightransf, 1, var, na.rm = TRUE) != 0 , ]

#Renamed curated col names for heatmap
comb2high <- sapply(colnames(highvar), function(x) paste(strsplit(x,"_", fixed = T)[[1]][2],
                                                         strsplit(x,"_",fixed = T)[[1]][3],
                                                         sep = "_"
))
colnames(highvar) <- comb2high


#############High throughput heatmap ##############
pdf("Highthroughput.pdf", width=7, height=7)
par(mar=c(2,2,2,2), cex=1.0)
heatmap.2(as.matrix(highvar), col=yb, scale="row", dendrogram = "row", labRow = "",
          key=TRUE, symkey=FALSE, density.info="none", trace="none",
          cexRow=0.3,lhei=c(0.20,0.70), lwid = c(0.25,0.5),
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x, method="ward.D2"))

dev.off()
