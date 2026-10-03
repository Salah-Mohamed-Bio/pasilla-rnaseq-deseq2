

# load the library

library(DESeq2)

BiocManager::install(c("org.Dm.eg.db","AnnotationDbi"))

library(AnnotationDbi)
library(org.Dm.eg.db)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")



if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")

library(clusterProfiler)
library(org.Dm.eg.db)

install.packages("gdtools")                   

if (!requireNamespace("EnhancedVolcano", quietly = TRUE))
  BiocManager::install("EnhancedVolcano")

library(EnhancedVolcano)

if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")

library(pheatmap)



#set workign directory 


count_dir <- "~/RNA_analysis/featurcount.d"

#load the date 

sample_files <- c("GSM461176_untreat_single_featureCounts.counts", "GSM461178_untreat_paired_featureCounts.counts", "GSM461180_treat_paired_featureCounts.counts", "GSM461182_untreat_single_featureCounts.counts",
                  "GSM461177_untreat_paired_featureCounts.counts", "GSM461179_treat_single_featureCounts.counts", "GSM461181_treat_paired_featureCounts.counts")
length(sample_files)
sample_names <- c("GSM461176_untreat","GSM461178_untreat","GSM461180_treat","GSM461182_untreat","GSM461177_untreat","GSM461179_treat","GSM461181_treat")

length(sample_names)
condition <- factor(c("untreated","untreated","treated","untreated","untreated","treated","treated"))
length(condition)
count_list <- lapply(sample_files, function(f){read.table(file.path(count_dir,f), header = FALSE, col.names = c("gene_id", "count"))})

print(count_list)
head(count_list)

counts <- Reduce(function(x,y)merge(x,y, by= "gene_id"), count_list)

gene_id <- counts$gene_id

print(gene_id)

counts <- counts[,-1]

counts[] <- lapply(counts, as.numeric)

str(counts)

print(counts)

colnames(counts) <- sample_names

rownames(counts) <-gene_id

colData <- data.frame(row.names = sample_names, condition = condition)

ncol(counts)


nrow(colData)

all(colnames(counts) == rownames(colData))

counts[is.na(counts)] <- 0

dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts) , colData = colData, design = ~ condition)

dds <- DESeq(dds)

res <- results(dds)

summary(res)

# data preparation 

red_df <- as.data.frame(res)


res$symbol <- mapIds(org.Dm.eg.db, keys =rownames(res), column="SYMBOL", keytype="FLYBASE", multiVals="first")

res$genenames <- mapIds(org.Dm.eg.db, keys = rownames(res), column = "GENENAME", keytype = "FLYBASE", multiVals="first")

res$entrez <- mapIds( org.Dm.eg.db,keys = rownames(res), column = "ENTREZID", keytype ="FLYBASE" , multiVals = "first")
head(res)
columns(org.Dm.eg.db)
res_ordered <- res[order(res$padj),]

head(res_ordered) 


res_df <- as.data.frame(res_ordered)

write.csv(res_df,"DE_results_annotated.csv")



significant_genes <- subset(res_df,padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(significant_genes, "significant_genes.csv")

nrow(significant_genes)


# MA plot

plotMA(res, ylim = c(-5, 5), main = "MA Plot: Treated vs Untreated")

# enhanced volcano 

EnhancedVolcano(res,
                lab = res$symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Treated vs Untreated',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 3.0)



# heatmap drawing 

norm_counts <- counts(dds, normalized = TRUE)

top_genes <- rownames(significant_genes)[1:20]

pheatmap(log2(norm_counts[top_genes, ] + 1),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         annotation_col = colData,
         main = "Top 20 Differentially Expressed Genes")



#GO enrichment

sig_genes <- rownames(significant_genes)


length(sig_genes)


go_results <- enrichGO(gene = sig_genes,
                       OrgDb = org.Dm.eg.db,
                       keyType = "FLYBASE",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.3,      
                       qvalueCutoff = 1)

nrow(as.data.frame(go_results))
head(as.data.frame(go_results), 10)



dotplot(go_results, showCategory = 7, title = "GO Enrichment: Biological Process")

dir.create("results", showWarnings = FALSE)

write.csv(as.data.frame(res), "results/deseq2_full_results.csv")
write.csv(as.data.frame(significant_genes), "results/significant_genes.csv")
write.csv(as.data.frame(go_results), "results/go_enrichment_results.csv")











