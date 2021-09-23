# CRC - KEGG

gsea-cli.sh GSEAPreranked -gmx ftp.broadinstitute.org://pub/gsea/gene_sets/c2.cp.kegg.v7.4.symbols.gmt \
-collapse No_Collapse 
-mode Max_probe \
-norm meandiv \
-nperm 1000 \
-rnk /Users/johnshelley/Dropbox/John/MD-PhD/PhD/G2/BMIF 6310/Project 1/HT29_RNASeq_ranks.rnk \
-scoring_scheme weighted \
-rpt_label CRC_KEGG \
-create_svgs false \
-help false \
-include_only_symbols true \
-make_sets true \
-plot_top_x 20 \
-rnd_seed timestamp \
-set_max 500 -set_min 15 \
-zip_report false \
-out /Users/johnshelley/gsea_home/output/sep22

# WT - KEGG

gsea-cli.sh GSEAPreranked -gmx ftp.broadinstitute.org://pub/gsea/gene_sets/c2.cp.kegg.v7.4.symbols.gmt \
-collapse No_Collapse 
-mode Max_probe \
-norm meandiv \
-nperm 1000 \
-rnk /Users/johnshelley/Dropbox/John/MD-PhD/PhD/G2/BMIF 6310/Project 1/NCM460_RNASeq_ranks.rnk \
-scoring_scheme weighted \
-rpt_label WT_KEGG \
-create_svgs false \
-help false \
-include_only_symbols true \
-make_sets true \
-plot_top_x 20 \
-rnd_seed timestamp \
-set_max 500 -set_min 15 \
-zip_report false \
-out /Users/johnshelley/gsea_home/output/sep22

# CRC - Oncogenic signatures

gsea-cli.sh GSEAPreranked -gmx ftp.broadinstitute.org://pub/gsea/gene_sets/c6.all.v7.4.symbols.gmt \
-collapse No_Collapse 
-mode Max_probe \
-norm meandiv \
-nperm 1000 \
-rnk /Users/johnshelley/Dropbox/John/MD-PhD/PhD/G2/BMIF 6310/Project 1/HT29_RNASeq_ranks.rnk \
-scoring_scheme weighted \
-rpt_label CRC_Oncogenic \
-create_svgs false \
-help false \
-include_only_symbols true \
-make_sets true \
-plot_top_x 20 \
-rnd_seed timestamp \
-set_max 500 -set_min 15 \
-zip_report false \
-out /Users/johnshelley/gsea_home/output/sep22

# WT - Oncogenic signatures

gsea-cli.sh GSEAPreranked -gmx ftp.broadinstitute.org://pub/gsea/gene_sets/c6.all.v7.4.symbols.gmt \
-collapse No_Collapse 
-mode Max_probe \
-norm meandiv \
-nperm 1000 \
-rnk /Users/johnshelley/Dropbox/John/MD-PhD/PhD/G2/BMIF 6310/Project 1/NCM460_RNASeq_ranks.rnk \
-scoring_scheme weighted \
-rpt_label WT_Oncogenic \
-create_svgs false \
-help false \
-include_only_symbols true \
-make_sets true \
-plot_top_x 20 \
-rnd_seed timestamp \
-set_max 500 -set_min 15 \
-zip_report false \
-out /Users/johnshelley/gsea_home/output/sep22