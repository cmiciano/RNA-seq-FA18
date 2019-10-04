##########################
# Creates table of DE genes
# from diffno3.txt (which was generated from rawno3.txt in getDiff
# rawno3 was made through analyzerepeats of tag directory)

# creates a table called tablediff which has 1050 DE genes
# across two conditions, one in which the third abf3 replicate is present
# and another in which the third replicate is removed
# the table contains the list of genes that are expressed in either of the 
# conditions. zerodiff contains the genes that are expressed in 3v3 but then
# aren't in 2v3


#Get all DE 3v3 subsetted from diffnobatch
#Get DE 2v3 from diffno3
#Get transcripts of both
#union both 2v3 and 3v3 
# get annotation of both
#determine if the uniongenelist gene is in the DE frame for certain cond
#zerodiff if certa


#########################
#FPKM values of 3v3 
fpkmyeast <- read.delim("~/Documents/salk/fpkmyeast.txt")
diffnobatch <- read.delim("~/Documents/salk/diffnobatch.txt")
names(diffnobatch)[1]<-"Transcript"


#FPKM and diff of no3
fpkmno3 <- read.delim("~/Documents/salk/no3/fpkmno3.txt")
diffno3 <- read.delim("~/Documents/salk/no3/diffno3.txt")

#Upreg and downreg in 3v3
#448 total DE genes
#161 UP 287 down
up3v3 <- read.delim("~/Documents/salk/no3/upregulatedgenes.txt")
down3v3 <- read.delim("~/Documents/salk/no3/downregulatedgenes.txt")


#subsetting from origi upreg
diffup3v3 <-subset(diffnobatch, diffnobatch$WT.vs..abf2.Log2.Fold.Change >= 0.5
                   & diffnobatch$WT.vs..abf2.adj..p.value <= 0.05)


diffdown3v3 <-subset(diffnobatch, diffnobatch$WT.vs..abf2.Log2.Fold.Change <= -0.5
                   & diffnobatch$WT.vs..abf2.adj..p.value <= 0.05)

#subsetting from original nobatch file
sigall3v3 <- subset(diffnobatch, abs(diffnobatch$WT.vs..abf2.Log2.Fold.Change) >= 0.5
                    & diffnobatch$WT.vs..abf2.adj..p.value <= 0.05)

#448 total DE genes
all3v3 <- rbind(up3v3, down3v3)

#subsetting from original nobatch file
#sigall3v3 <- subset(diffnobatch, abs(diffnobatch$WT.vs..abf2.Log2.Fold.Change) >= 0.5
#                    & diffnobatch$WT.vs..abf2.adj..p.value <= 0.05)

#subsetting from origi upreg
#diffup3v3 <-subset(diffnobatch, diffnobatch$WT.vs..abf2.Log2.Fold.Change >= 0.5
#                   & diffnobatch$WT.vs..abf2.adj..p.value <= 0.05)

#### 2v3 ##########

#Change rownames to gene symbol
names(diffno3)[1]<-"Transcript"
ann2v3 <- as.character(diffno3$Annotation.Divergence)
sym2v3 <- sapply(ann2v3, function(x) strsplit(x,"|", fixed = T)[[1]][1])
onlysym2v3 <- as.character(sym2v3)

diffno3 <- cbind(onlysym2v3, diffno3)
colnames(diffno3)[1] <- "Symbol"

#Number of sig genes when abf3 is removed is 1007 if FC cutoff 0.5
#403 0.5-FC 33 1-FC 
diffupsig <-subset(diffno3, diffno3$WT.vs..abf2.Log2.Fold.Change >= 0.5
                  & diffno3$WT.vs..abf2.adj..p.value <= 0.05)

############# Upregulated GENES #############

#diffbatch fold change and adj p-value find genes in diffbatch

#403 0.5-FC 33 1-FC
#write.table(diffupsig,"upregulatedgenes2v3.txt",sep="\t",row.names=F, col.names = TRUE, quote=FALSE)


#604 -0.5-FC 83 1-FC
diffdownsig <-subset(diffno3, diffno3$WT.vs..abf2.Log2.Fold.Change <= -0.5
                   & diffno3$WT.vs..abf2.adj..p.value <= 0.05)


#### down regulated genes######
#604 -0.5-FC 83 1-FC
#write.table(diffdownsig,"downregulatedgenes2v3.txt",sep="\t",row.names=F, col.names = TRUE, quote=FALSE)



all2v3comb<- subset(diffno3, abs(diffno3$WT.vs..abf2.Log2.Fold.Change) >= 0.5
                        & diffno3$WT.vs..abf2.adj..p.value <= 0.05)

#all 1007 genes
#write.table(all2v3comb,"allgenes2v3.txt",sep="\t",row.names=F, col.names = TRUE, quote=FALSE)


#all2v3 <- rbind(diffupsig, diffdownsig)
#Compare this with subsetting diffno3 based on abs fold change same 
#preserves order of diffno3 instead of changing it

#Differentially expressed genes found in 3v3
diff3v3 <- union(up3v3$Transcript, down3v3$Transcript)#get all up then down
diffan3v3 <- union(up3v3$Annotation.Divergence,
                   down3v3$Annotation.Divergence) #new

diff3v3 <- as.character(sigall3v3$Transcript)# new

#Differentially expressed genes found in 2v3
diff2v3 <- union(diffupsig$Transcript, diffdownsig$Transcript)
diffan2v3 <- union(diffupsig$Annotation.Divergence, 
              diffdownsig$Annotation.Divergence) #new


#diff2v3 <- as.character(all2v3comb$Transcript)# new
#Get total list of expressed genes across both conditions
#Transcript names 
un2and3 <- union(diff2v3, diff3v3)
#genesinsub <- subset(diff)


#Annotations of both conditions (might be wrong since only got diffno3 not both)
#ann2and3 <- subset(as.character(diffno3$Annotation.Divergence), diffno3$Transcript %in% un2and3)
ann2and3 <- union(diffan2v3, diffan3v3) #new

#If diff gene expressed in 2v3 (all of them should)
vec2v3 <- ifelse(un2and3 %in% diff2v3, "yes", "no")

#Check if the geneis differentially expressed in 3v3
vec3v3 <- ifelse(un2and3 %in% diff3v3, "yes", 'no')

#Table of genes that are differentially expressed across both conditions
tablediff <- data.frame(un2and3, ann2and3, vec2v3, vec3v3)

##43 genes in which they were significant when abf3 present and then
## not significant when abf3 was removed
zerodiff <- subset(tablediff, tablediff$vec2v3 == "no" & 
                     tablediff$vec3v3 == "yes")
#FC of 43 genes
#fc43nosig <- subset(all3v3$WT.vs..abf2.Log2.Fold.Change,
  #                   all3v3$Transcript %in% zerodiff$un2and3)

#Find corresponding FC of zerodiff transcripts in diffno3

#fc43nosig <- subset(diffno3,
   #                  diffno3$Transcript %in% zerodiff$un2and3)


# ##Changed fc should not be sig
 fc43nosig <- subset(diffno3,
                     diffno3$Transcript %in% zerodiff$un2and3)
 
transandfc <- fc43nosig[, c("Transcript", "WT.vs..abf2.Log2.Fold.Change")]
loc70 <- which(transandfc) 
# fc <- ifelse(zerodiff$un2and3 %in% transandfc$Transcript,
#              transandfc$WT.vs..abf2.Log2.Fold.Change, "no")
# 
# zerodiff$fold_change_in_orig <- fc


zerodiff$fold_change_in_orig <- fc43nosig

colnames(tablediff) <- c("Gene Symbol", "Gene Annotation", 
                         "Gene DE when abf3 is removed", 
                         "Gene DE when abf3 is present")
colnames(zerodiff) <- c("Gene Symbol", "Gene Annotation", 
                         "Gene DE when abf3 is removed", 
                         "Gene DE when abf3 is present",
                         "FC when abf3 is present")

#Write 1007 DE genes to table
#write.table(tablediff,"1007DEgenesnoabf3rev.txt",sep="\t",row.names=F, col.names = TRUE, quote=FALSE)
dev.off()

#Write 43 genes that will not be DE when removed
#write.table(zerodiff,"43genesnoDEnoabf3rev.txt",sep="\t",row.names=F, col.names = TRUE, quote=FALSE)









