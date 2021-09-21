library(edgeR)

setwd("/Users/monikagrabowska/Downloads")

rna <- read.table("GSE179862_Raw_gene_counts_matrix.txt", 
               header = TRUE, sep = "\t", quote="\"", stringsAsFactors = FALSE)

rna_raw <- rna[,c(2:5)] # extract raw counts
gene_names <- paste(rna[,10], rna[,1], sep = "|") # HUGO gene symbol | Ensembl ID
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
mds_output <- plotMDS(d, pch = 1, labels = d$samples$group, col = c("darkgreen","blue","red", "orange")[d$samples$group], xlim = c(-4,4), ylim = c(-2,2))
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

# creating a ranked gene list for input into GSEA 
# HT-29
ranks_HT29 <- sign(tt_HT29$table$logFC) * -log10(tt_HT29$table$PValue)

gene_symbols_HT29 <- unlist(lapply(rownames(tt_HT29$table), function(data) 
                        {unlist(strsplit(data,"\\|"))[1]}))
gene_ids_HT29 <- unlist(lapply(rownames(tt_HT29$table), function(data) 
                        {unlist(strsplit(data,"\\|"))[2]})) 

ranks_HT29 <- cbind(gene_symbols_HT29, ranks_HT29) # need HUGO gene symbols for GSEA input

colnames(ranks_HT29) <- c("HUGOGeneSymbol","Rank")
ranks_HT29 <- ranks_HT29[order(as.numeric(ranks_HT29[,2]),decreasing = TRUE),]

write.table(ranks_HT29, file.path("/Users/monikagrabowska/Downloads/HT29_RNASeq_ranks.rnk"), 
            col.names = TRUE, sep="\t", row.names = FALSE, quote = FALSE) # write to .rnk file 

#NCM460
ranks_NCM460 = sign(tt_NCM460$table$logFC) * -log10(tt_NCM460$table$PValue)

gene_symbols_NCM460 <- unlist(lapply(rownames(tt_NCM460$table), function(data) 
                          {unlist(strsplit(data,"\\|"))[1]}))
gene_ids_NCM460 <- unlist(lapply(rownames(tt_NCM460$table), function(data) 
                          {unlist(strsplit(data,"\\|"))[2]})) 

ranks_NCM460 <- cbind(gene_symbols_NCM460, ranks_NCM460) # need HUGO gene symbols for GSEA input

colnames(ranks_NCM460) <- c("HUGOGeneSymbol","Rank")
ranks_NCM460 <- ranks_NCM460[order(as.numeric(ranks_NCM460[,2]),decreasing = TRUE),]

write.table(ranks_NCM460, file.path("/Users/monikagrabowska/Downloads/NCM460_RNASeq_ranks.rnk"), 
            col.names = TRUE, sep="\t", row.names = FALSE, quote = FALSE) # write to .rnk file 


