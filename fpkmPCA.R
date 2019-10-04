##########################
# Creates PCA plot of genes from fpkm.txt file 
#
#########################



fpkmyeast <- read.delim("~/Documents/salk/fpkmyeast.txt")
# yeastmtDNA <- read.csv("~/Documents/salk/yeastmtDNA.txt")
# chr1yeast <- read.csv("~/Documents/salk/chr1yeast.txt")
# rawyeast <- read.delim("~/Documents/salk/rawyeast.txt")

#transpose fpkm file, take numeric values and rename rows as gene names
#numeric_fpkm <- t(fpkmyeast[,9:14]) #samples x genes
#View(numeric_fpkm[1:6, 1:10])
transformed <- t(log(fpkmyeast[,9:14]+5)) #samples by genes

genenames <- as.character(fpkmyeast$Annotation.Divergence)
singlename <- sapply(genenames, function(x) strsplit(x,"|", fixed = T)[[1]][1])
rownames(fpkmyeast) <- singlename
colnames(transformed) <- singlename
#View(transformed[1:6,1:5])



# #See if there are chr 1 genes and mtDNA genes in fpkm
# sum(yeastmtDNA$Gene.stable.ID %in% fpkmyeast[,1])
# which(yeastmtDNA$Gene.stable.ID %in% fpkmyeast[,1])
# 
# 
# #Get gene names with transformed data
# whichthings <- transformed[, which(colSums(transformed) != 0)]
# columnSum <- colSums(transformed)
# any(columnSum == 0)

#Remove columns that have zero variance #sample by gene
no_zeroes <- transformed[ , apply(transformed, 2, var, na.rm = TRUE) != 0]
tencol <- no_zeroes[1:6,1:10]
# 
# #Sum of expression top 500 genes
# no_zeroes_t <- t(no_zeroes) #Genes by samples
# sumgene <- rowSums(no_zeroes_t)
# sumgene <- sumgene[order(-sumgene)]
# sumgene_500 <- sumgene[1:500]
# sumnames500 <- names(sumgene_500)
# 
# trans_sum_500 <- no_zeroes[,sumnames500]
# 
# #variance top 500 genes with most variance
# #transformed_var <- no_zeroes[,apply(no_zeroes, 2, var)] 
# var_apply <- apply(no_zeroes, 2, var)
# var_apply <- var_apply[order(-var_apply)]
# var_500 <- var_apply[1:500]
# varnames500 <- names(var_500)
# trans_var_500 <- no_zeroes[,varnames500]
# 



#Original FPKM plot        
plot1 <- prcomp(no_zeroes, center = T, scale. = T)


#compute standard deviation of each principal component
std_dev <- plot1$sdev

#compute variance
pr_var <- std_dev^2

#proportion of variance explained
prop_varex <- pr_var/sum(pr_var)

pc1var <- prop_varex[1]
pc2var <- prop_varex[2]

#plot(plot1[["x"]][,1],plot1[["x"]][,2], xlab= "PC1", ylab = "PC2",
#     xlim = c(-100,100), ylim =c(-100,100))
#text(plot1[["x"]][,1],plot1[["x"]][,2], labels = rownames(plot1[["x"]]), 
#     xlim = c(-100,100), ylim =c(-100,100))
#plot1PC <- cbind(no_zeroes, plot1$x[,1:2])

library(ggplot2)

plot1DF <- data.frame(plot1[["x"]])
Grouping <- sapply(rownames(plot1DF), function(x) strsplit(x,"_", fixed = T)[[1]][2])
groupnum <- sapply(rownames(plot1DF), function(x) strsplit(x,"_", fixed = T)[[1]][3])
combined_grp <- paste(Grouping, groupnum, sep = "_")
plot1DF <- cbind(plot1DF,Grouping,groupnum,combined_grp)
ggplot(plot1DF, aes(PC1 ,PC2, col = Grouping)) +
  xlab(paste("PC1",round(pc1var * 100, 0), "%"))+
  ylab(paste("PC2",round(pc2var * 100, 0), "%")) +
  xlim(-100, 100) +
  geom_text(aes(label= combined_grp), size = 5) +
  scale_fill_discrete(name = "Grouping") + 
  theme_classic()
ggsave("ggplotyeast.jpeg", width = 10, height = 10, units = "cm", dpi = 600)

class(plot1[["x"]])
View(plot1[["x"]])

library(rgl)
plot3d(plot1[["x"]], col=plot1DF$groupnum)
text3d(plot1DF,texts=plot1DF$combined_grp)


# #500 genes with highest variance
# varsum <- prcomp(trans_var_500,center = T, scale. = T)
# plot(varsum[["x"]][,1],varsum[["x"]][,2], xlab= "PC1", ylab = "PC2",
#      xlim = c(-50,50), ylim =c(-50,50))
# #text(varsum[["x"]][,1],varsum[["x"]][,2], labels = rownames(varsum[["x"]]))
# 
# 
# #500 genes with highest expression
# varsum <- prcomp(trans_sum_500,center = T, scale. = T)
# plot(plot1[["x"]][,1],plot1[["x"]][,2], xlab= "PC1", ylab = "PC2",
#      xlim = c(-100,100), ylim =c(-100,100))
# text(plot1[["x"]][,1],plot1[["x"]][,2], labels = rownames(plot1[["x"]]), 
#      xlim = c(-100,100), ylim =c(-100,100))


