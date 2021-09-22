library(edgeR)
library(EnhancedVolcano)

setwd("/Users/monikagrabowska/Downloads")

rna <- read.table("GSE179862_Raw_gene_counts_matrix.txt", 
               header = TRUE, sep = "\t", quote="\"", stringsAsFactors = FALSE)

rna_raw <- rna[,c(2:5)] # extract raw counts
gene_names <- paste(rna[,10], rna[,1], sep = "|") # HGNC gene symbol | Ensembl ID
row.names(rna_raw) <- make.unique(gene_names, sep = "|") 
rna_raw <- as.matrix(rna_raw)

group = factor(c("HT29_control","HT29_treat","NCM460_control","NCM460_treat"))

# filter - have data for 64962 genes before filtering --> reduced to 17125 after filtering 
keep <- rowSums(cpm(rna_raw)) > 0.5
rna_raw <- rna_raw[keep,]

d <- DGEList(counts = rna_raw, group = group)

# normalization 
d <- calcNormFactors(d)

# create MDS plot
png(filename = "/Users/monikagrabowska/Downloads/mdsplot_GSE179862.png")
mds_output <- plotMDS(d, pch = 19, labels = NULL, 
              col = c("darkgreen","blue","red", "orange")[d$samples$group], 
              xlim = c(-4,4), ylim = c(-2,2), xlab = "Dimension 1", ylab = "Dimension 2", cex.lab=1.25, cex.axis=1.5)
legend("topleft", 
       legend=levels(d$samples$group), 
       pch = 19, col= c("darkgreen","blue","red", "orange"), title = "Condition",  
       bty = 'n')
dev.off()

# cannot calculate dispersion (no replicates), so set bcv manually at 0.2 
bcv <- 0.2

# excluding unannoted genes
exclude <- grep("^\\|", rownames(d), value = T)
d <- d[which(!rownames(d) %in% exclude),]

# calculate differential expression for HT29 cells
de_HT29 <- exactTest(d, pair = c("HT29_control","HT29_treat"), dispersion=bcv^2)
tt_exact_test_HT29 <- topTags(de_HT29, n = nrow(d))
tt_HT29 <- tt_exact_test_HT29

# calculate differential expression for NCM460 cells
de_NCM460 <- exactTest(d, pair = c("NCM460_control", "NCM460_treat"), dispersion=bcv^2)
tt_exact_test_NCM460 <- topTags(de_NCM460,n = nrow(d))
tt_NCM460 <- tt_exact_test_NCM460

# select genes with fdr < 0.05 
select_genes_HT29 <- which(tt_HT29$table$FDR < 0.05)
length(select_genes_HT29) # 487 genes for HT-29 cells

select_genes_NCM460 <- which(tt_NCM460$table$FDR < 0.05)
length(select_genes_NCM460) # 1260 genes for NCM460 cells

# separate gene name into HGNC gene symbols and Ensembl ID
gene_symbols_HT29 <- unlist(lapply(rownames(tt_HT29$table), function(data) 
{unlist(strsplit(data,"\\|"))[1]}))
gene_ids_HT29 <- unlist(lapply(rownames(tt_HT29$table), function(data) 
{unlist(strsplit(data,"\\|"))[2]})) 

gene_symbols_NCM460 <- unlist(lapply(rownames(tt_NCM460$table), function(data) 
{unlist(strsplit(data,"\\|"))[1]}))
gene_ids_NCM460 <- unlist(lapply(rownames(tt_NCM460$table), function(data) 
{unlist(strsplit(data,"\\|"))[2]})) 

# volcano plot for visualizing differential expression results
res_HT29 <- cbind(tt_HT29$table$logFC, tt_HT29$table$PValue)
rownames(res_HT29) <- gene_symbols_HT29
colnames(res_HT29) <- c("logFC","pvalue")

volcano_HT29 <- EnhancedVolcano(res_HT29, lab = rownames(res_HT29), x = "logFC", y = "pvalue",
                title = "HT-29", selectLab = c("CLDN2", "CLEC3A", "SNORD3A"), 
                cutoffLineType = "blank", pCutoff = 0.05 / nrow(res_HT29), FCcutoff = 1,
                col = c("black", "black", "black", "red3"), colAlpha = 1, legendPosition = "none")

res_NCM460 <- cbind(tt_NCM460$table$logFC, tt_NCM460$table$PValue)
rownames(res_NCM460) <- gene_symbols_NCM460
colnames(res_NCM460) <- c("logFC","pvalue")
volcano_NCM460 <- EnhancedVolcano(res_NCM460, lab = rownames(res_NCM460), x = "logFC", y = "pvalue", 
                  title = "NCM460", selectLab = c("ARHGAP5", "CNOT7", "ARHGAP29"), 
                  cutoffLineType = "blank", pCutoff = 0.05 / nrow(res_NCM460), FCcutoff = 1,
                  col = c("black", "black", "black", "red3"), colAlpha = 1, legendPosition = "none")

# creating a ranked gene list for input into GSEA 
# HT-29
ranks_HT29 <- sign(tt_HT29$table$logFC) * -log10(tt_HT29$table$PValue)
ranks_HT29 <- cbind(gene_symbols_HT29, ranks_HT29) # need HGNC gene symbols for GSEA input

colnames(ranks_HT29) <- c("HGNCGeneSymbol","Rank")
ranks_HT29 <- ranks_HT29[order(as.numeric(ranks_HT29[,2]),decreasing = TRUE),]

write.table(ranks_HT29, file.path("/Users/monikagrabowska/Downloads/HT29_RNASeq_ranks.rnk"), 
            col.names = TRUE, sep="\t", row.names = FALSE, quote = FALSE) # write to .rnk file 

#NCM460
ranks_NCM460 = sign(tt_NCM460$table$logFC) * -log10(tt_NCM460$table$PValue)
ranks_NCM460 <- cbind(gene_symbols_NCM460, ranks_NCM460) # need HGNC gene symbols for GSEA input

colnames(ranks_NCM460) <- c("HGNCGeneSymbol","Rank")
ranks_NCM460 <- ranks_NCM460[order(as.numeric(ranks_NCM460[,2]),decreasing = TRUE),]

write.table(ranks_NCM460, file.path("/Users/monikagrabowska/Downloads/NCM460_RNASeq_ranks.rnk"), 
            col.names = TRUE, sep="\t", row.names = FALSE, quote = FALSE) # write to .rnk file 


