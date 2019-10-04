diff <- read.delim("~/Documents/salk/diff.txt")
names(diff)[1]<-"Transcript"

#Overenrichment Analysis List

#changed fold change to 1 and pvalue to 0.05
ORAfc1 <- subset(diff, diff$WT.vs..KO.Log2.Fold.Change >= 1 & diff$WT.vs..KO.adj..p.value <= 0.05 )
ORAnamesfc1 <- data.frame(ORAfc1$Transcript)

#changed fold change to 0.5 and pvalue to 0.05
ORAfc5 <- subset(diff, diff$WT.vs..KO.Log2.Fold.Change >= 0.5 & diff$WT.vs..KO.adj..p.value <= 0.05 )
ORAnamesfc5 <- data.frame(ORAfc5$Transcript)

#both 0.05 and 1 fold change cutoffs are the same
fc <- setequal(ORAnamesfc1, ORAnamesfc5)

write.table(ORAnamesfc1,"ORAnamesfc.txt",sep="\t",row.names=FALSE, col.names = FALSE, quote=FALSE)

#fc 1 pvalue 0.01
ORAfc1p1 <- subset(diff, diff$WT.vs..KO.Log2.Fold.Change >= 1 & diff$WT.vs..KO.adj..p.value <= 0.01 )
ORAnamesfc1p1 <- data.frame(ORAfc1p1$Transcript)

#fc 0.5 pvalue 0.01
ORAfc5p1 <- subset(diff, diff$WT.vs..KO.Log2.Fold.Change >= 0.5 & diff$WT.vs..KO.adj..p.value <= 0.01 )
ORAnamesfc5p1 <- data.frame(ORAfc5p1$Transcript)

pval <- setequal(ORAnamesfc1p1, ORAnamesfc5p1)

write.table(ORAnamesfc1p1,"ORAnamespval.txt",sep="\t",row.names=FALSE, col.names = FALSE, quote=FALSE)



#Gene Set Analysis
diff$sign <- ifelse(diff$WT.vs..KO.Log2.Fold.Change > 0, 1, -1)
diff$logpvalue <- log10(diff$WT.vs..KO.adj..p.value)
diff$signtimeslog <- diff$sign * diff$logpvalue


diff <- diff[order(-diff$signtimeslog),]

sortedIDs <- data.frame(diff$"Transcript/RepeatID", diff$signtimeslog)
write.table(sortedIDs,"genenamesGSA.txt",sep="\t",row.names=FALSE, quote=FALSE)
