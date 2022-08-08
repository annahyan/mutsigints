# from https://github.com/CRI-iAtlas/ImmuneSubtypeClassifier

library(ImmuneSubtypeClassifier)
library(here)

data.dir = here("data", "expression")

out.dir = here("data/RDS/PCAWG/immune_states")

pcawg.data = read.delim(file.path(data.dir, 
                                  "tophat_star_fpkm_uq.v2_aliquot_gl_filtered.tsv"),
                        check.names = FALSE)

rownames(pcawg.data) = pcawg.data$feature
pcawg.data$feature = NULL

pcawg.data[1:4, 1:4]

pcawg.data.immgenes = pcawg.data[ebpp_genes_sig$Ensembl, ]
pcawg.data.immgenes

immgenes.symbol = ebpp_genes_sig[match(rownames(pcawg.data.immgenes), 
               ebpp_genes_sig$Ensembl), ]$Symbol

pcawg.data.immgenes.filtered = pcawg.data.immgenes[-which(is.na(immgenes.symbol)), ]
rownames(pcawg.data.immgenes.filtered) = immgenes.symbol[-which(is.na(immgenes.symbol))]
pcawg.data.immgenes = pcawg.data.immgenes.filtered
rm(pcawg.data.immgenes.filtered)

pcawg.data.immgenes = as.matrix(pcawg.data.immgenes)

calls <- callEnsemble(X=pcawg.data.immgenes, geneids='symbol')

saveRDS(calls, file = file.path(out.dir, "pcawg.xgboost.out.RDS"))
