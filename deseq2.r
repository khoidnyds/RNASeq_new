# Try to describe the codes from the DESeq2 to generate 
# different expressed gene list, 
# MA-plot
# PCA plot, 
# heatmap for the DE genes and 
# ggplot to display any gene expression level (DESeq2 vignette will be helpful)


library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(rlist)
pvalue_cutoff = 0.01

# read features_matrix (49677 x 54) 
cts <- read.csv("features_matrix.txt", sep = '\t', header=T,  row.names = "Geneid")
cts <- as.matrix(cts)

# read conditions (54 x 1)
coldata <- read.csv("conditions.txt", row.names = 1)
coldata$condition <- factor(coldata$condition)
colnames(cts) <- rownames(coldata)

# create DESeqDataSet obj 
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "normal_colon")
# run DE 
dds <- DESeq(dds)
#pre-filtering
keep <- rowSums(counts(dds))>=10
dds <- dds[keep,]

resultsNames(dds)
res_pc <- results(dds, contrast=c("condition","primary_colorectal_cancer","normal_colon"))
res_mc <- results(dds, contrast=c("condition","metastasized_cancer","normal_colon"))

# remove row if pvalue is N/A
res_pc <- na.omit(res_pc)
res_mc <- na.omit(res_mc)

# plot MA, highlight the dot having pvalue smaller than the cutoff
generate_MA <- function(title, input, cutoff, plot_title="") {
  pdf(file = title, width = 12, height = 8)
  par(mar = c(5, 5, 5, 5))
  total_genes = nrow(input)
  sig_genes = nrow(input[input$padj<cutoff,])
  title = paste("MA plot. Total genes: ", total_genes,". Number of genes having the adjust pvalue < 0.01 
                (blue dots) is ", sig_genes," - ",format(sig_genes/total_genes*100, digits=4),"%", 
                plot_title, sep="")
  plotMA(input, alpha=cutoff, main=title, ylim=c(-4, 4))
  dev.off()
}
generate_MA("MA plot - primary colon cancer vs. normal colon.pdf",
                res_pc,
                pvalue_cutoff)
generate_MA("MA plot - metastasized cancer vs. normal colon.pdf",
                res_mc,
                pvalue_cutoff)

res_pc['abs_log2FC'] <- abs(res_pc$log2FoldChange)
res_mc['abs_log2FC'] <- abs(res_mc$log2FoldChange)

# select only significant gene based on pvalue threshold and log2FoldChange threshold
res_pcSig <- res_pc[res_pc$padj < pvalue_cutoff, ]
res_pcSig <- res_pcSig[res_pcSig$abs_log2FC > 4, ]
res_pcSig <- res_pcSig[order(res_pcSig$pvalue),]
write.csv(as.data.frame(res_pcSig), file="pc_results.csv")
write(rownames(res_pcSig), "DE_pc.txt", sep = "\n") 

res_mcSig <- res_mc[res_mc$padj < pvalue_cutoff, ]
res_mcSig <- res_mcSig[res_mcSig$abs_log2FC > 8, ]
res_mcSig <- res_mcSig[order(res_mcSig$pvalue),]
write.csv(as.data.frame(res_mcSig), file="mc_results.csv")
write(rownames(res_mcSig), "DE_mc.txt", sep = "\n") 


# Visualize heatmap of the count matrix
vsd <- vst(dds, blind=FALSE)

pdf(file = "Heatmap of count matrix (Primary colon cancer vs. Normal colon).pdf", width = 12, height = 8)
res_pcOrdered <- res_pc[order(res_pc$padj),]
select_genes <- order(rownames(subset(res_pcOrdered, padj < 0.01 & abs_log2FC > 4)),decreasing=TRUE)
df <- as.data.frame(colData(dds)[,c("condition")])
colnames(df) <- c("Conditions")
rownames(df) <- rownames(coldata)
select_cols <- rownames(subset(df, Conditions!="metastasized_cancer"))
df <- subset(df, Conditions!="metastasized_cancer")
pheatmap(assay(vsd)[select_genes,select_cols], 
         cluster_rows=FALSE, 
         cluster_cols=TRUE,
         show_rownames=FALSE, 
         show_colnames=FALSE, 
         annotation_col=df,
         main="Heatmap of count matrix (Primary colon cancer vs. Normal colon).")
dev.off()


pdf(file = "Heatmap of count matrix (Metastasized cancer vs. Normal colon).pdf", width = 12, height = 8)
res_pcOrdered <- res_mc[order(res_mc$padj),]
select_genes <- order(rownames(subset(res_pcOrdered, padj < 0.01 & abs_log2FC > 8)),decreasing=TRUE)
df <- as.data.frame(colData(dds)[,c("condition")])
colnames(df) <- c("Conditions")
rownames(df) <- rownames(coldata)
select_cols <- rownames(subset(df, Conditions!="primary_colorectal_cancer"))
df <- subset(df, Conditions!="primary_colorectal_cancer")
pheatmap(assay(vsd)[select_genes,select_cols], 
         cluster_rows=FALSE, 
         cluster_cols=TRUE,
         show_rownames=FALSE, 
         show_colnames=FALSE, 
         annotation_col=df,
         main="Heatmap of count matrix (Metastasized cancer vs. Normal colon).")
dev.off()


# Visualize heatmap of the sample-to-sample distances
pdf(file = "Heatmap sample-to-sample distances.pdf", width = 12, height = 8)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
vsd_label <- as.character(vsd$condition)
vsd_label <- replace(vsd_label, vsd_label=="primary_colorectal_cancer", "pcc")
vsd_label <- replace(vsd_label, vsd_label=="normal_colon", "norm")
vsd_label <- replace(vsd_label, vsd_label=="metastasized_cancer", "meta")
rownames(sampleDistMatrix) <- paste(colnames(vsd),vsd_label, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         main="Heatmap of sample-to-sample distances.")
dev.off()

# PCA of samples
pdf(file = "PCA.pdf", width = 12, height = 8)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) + 
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of the samples. \nPC1 explains 60% variance of the data. PC2 explains 15% variance of the data.")
coord_fixed() 
dev.off()

# Display gene expression level of any gene
generate_plot_counts <- function(gene_name) {
  jpeg(file = paste("plot_count_",gene_name,".jpg",sep=""),  width = 1000, height = 600, res=100)
  plotCounts(dds, gene=gene_name, intgroup="condition")
  dev.off()
}
generate_plot_counts('INHBA')
generate_plot_counts('COL11A1')
generate_plot_counts('ADAM12')
