library(msigdb)
library(ExperimentHub)
library(GSEABase)

eh = ExperimentHub()
query(eh , 'msigdb')


## Sources of data
## https://gdc.cancer.gov/about-data/publications/panimmune

## https://www.biorxiv.org/content/10.1101/2020.01.17.910950v1.full.pdf

library(ImmuneSubtypeClassifier)
data("geneSetSymbols")
names(genesetsymbols)

immune.gene.lists = genesetsymbols

data.dir = here("data", "expression")
out.dir = here("data/RDS/PCAWG/immune_states")

### The function borrowed from https://rpubs.com/pranali018/SSGSEA
ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
    row_names = rownames(X)
    num_genes = nrow(X)
    gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
    
    # Ranks for genes
    R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
    
    # Calculate enrichment score (es) for each sample (column)
    es = apply(R, 2, function(R_col) {
        gene_ranks = order(R_col, decreasing = TRUE)
        
        # Calc es for each gene set
        es_sample = sapply(gene_sets, function(gene_set_idx) {
            # pos: match (within the gene set)
            # neg: non-match (outside the gene set)
            indicator_pos = gene_ranks %in% gene_set_idx
            indicator_neg = !indicator_pos
            
            rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
            
            step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
            step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
            
            step_cdf_diff = step_cdf_pos - step_cdf_neg
            
            # Normalize by gene number
            if (scale) step_cdf_diff = step_cdf_diff / num_genes
            
            # Use ssGSEA or not
            if (single) {
                sum(step_cdf_diff)
            } else {
                step_cdf_diff[which.max(abs(step_cdf_diff))]
            }
        })
        unlist(es_sample)
    })
    
    if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
    
    # Normalize by absolute diff between max and min
    if (norm) es = es / diff(range(es))
    
    # Prepare output
    rownames(es) = names(gene_sets)
    colnames(es) = colnames(X)
    return(es)
}

pcawg.data = read.delim(file.path(data.dir, 
                                  "tophat_star_fpkm_uq.v2_aliquot_gl_filtered.tsv"),
                        check.names = FALSE)

rownames(pcawg.data) = pcawg.data$feature
pcawg.ensembl = pcawg.data$feature
pcawg.data$feature = NULL
pcawg.data = as.matrix(pcawg.data)
pcawg.data[1:4, 1:4]

library("org.Hs.eg.db") # remember to install it if you don't have it already
symbols <- mapIds(org.Hs.eg.db, keys = pcawg.data$feature, keytype = "ENSEMBL", column="SYMBOL", )


pcawg.median = apply(pcawg.data, MARGIN = 1, median)
head(pcawg.median)
pcawg.mean = apply(pcawg.data, MARGIN = 1, mean)
head(pcawg.mean)

medians.p = pcawg.median %>% 
    enframe() %>% 
    ggplot(aes(x = value)) + 
        geom_histogram(binwidth = 0.1) + 
        coord_cartesian(x = c(0, 30)) + 
    ggtitle("Ensembl gene medians")

means.p = pcawg.mean %>% 
    enframe() %>% 
    ggplot(aes(x = value)) + 
    geom_histogram(binwidth = 0.1) + 
    coord_cartesian(x = c(0, 30)) +
    ggtitle("Ensembl gene means")


pp = ggarrange(means.p, medians.p)

pcawg.exp.prop = rowSums(pcawg.data > 0) / ncol(pcawg.data)

pcawg.filtered = pcawg.data[pcawg.mean > 1, ]
filtered.symbols <- mapIds(org.Hs.eg.db, 
                           keys = rownames(pcawg.filtered), keytype = "ENSEMBL", column="SYMBOL", )

pcawg.filtered = pcawg.filtered[!is.na(filtered.symbols), ]
filtered.symbols = filtered.symbols[ !is.na(filtered.symbols)]

filtered.symbols = filtered.symbols[!duplicated(filtered.symbols)]

pcawg.filtered = pcawg.filtered[ names(filtered.symbols), ]
rownames(pcawg.filtered) = filtered.symbols

gsva_out = gsva(as.matrix(pcawg.filtered), immune.gene.lists, method = "ssgsea")


saveRDS(gsva_out, file = file.path(out.dir, "ssGSEA_pcawg.RDS") )
