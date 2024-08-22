library(here)
source(here("R/load_packages.R"))
library(openxlsx)

fig.dir = here("figures/sig_sig_analysis")
data.out.dir = here("data/RDS/PCAWG", "sig_sig_ints")

for (dir.path in c(fig.dir, data.out.dir) ) {
    if (! file.exists(dir.path)) {
        dir.create(dir.path)
    }
}

PCAWG.full.subset.ann = readRDS(here("data", "RDS", "PCAWG", "signatures",
                                     "PCAWG.full.subset.ann.RDS"))

PCAWG.null.dists = readRDS(file.path(here("data", "RDS", "PCAWG", "signatures"),
                                      "PCAWG.null.model.metrics.dist.RDS"))

tissues = names(PCAWG.null.dists)
### Cooccurrence

# PCAWG.cooccurrence = calculate_bootstrap_metric(PCAWG.full.subset.ann, cooccurrence,
#                  min.tissue.samples = 20,
#                  sample.rate = 0.9, sample.N = 100,
#                  N = 1, seed = 1, p.adjust = TRUE)
# 
# saveRDS(PCAWG.cooccurrence,
#         file = file.path(data.out.dir, "PCAWG.cooccurrence.RDS") )

# # Plotting individual experiments
# # pp = plot_all_experiments(PCAWG.cooccurrence$Biliary_AdenoCA, "title")
# 
# PCAWG.cooccurrence.summary.nets = lapply(PCAWG.cooccurrence, summarize_calcs )
# 
# saveRDS(PCAWG.cooccurrence.summary.nets,
#         file = file.path(data.out.dir, "PCAWG.cooccurrence.summary.nets.RDS") )
# 
# ### Mutual information(BCMI)
# 
# PCAWG.bcmi = calculate_bootstrap_metric(PCAWG.full.subset.ann, bcmi,
#                                       min.tissue.samples = 20,
#                                       sample.rate = 0.9, sample.N = 100,
#                                       N = 1, seed = 1, p.adjust = TRUE)
# 
# saveRDS(PCAWG.bcmi,
#         file = file.path(data.out.dir, "PCAWG.bcmi.RDS") )

# #Plotting individual experiments
# # pp = plot_all_experiments(PCAWG.bcmi$Biliary_AdenoCA, "title")
# 
# PCAWG.bcmi.summary.nets = lapply(PCAWG.bcmi, summarize_calcs )
# 
# saveRDS(PCAWG.bcmi.summary.nets,
#         file = file.path(data.out.dir, "PCAWG.bcmi.summary.nets.RDS") )

# ### Pearson correlation
# 
# PCAWG.pearson = calculate_bootstrap_metric(PCAWG.full.subset.ann, cor_sigs,
#                               min.tissue.samples = 20,
#                               sample.rate = 0.9, sample.N = 100,
#                               N = 1, seed = 1, p.adjust = TRUE,
#                               method = "pearson")
# 
# saveRDS(PCAWG.pearson,
#         file = file.path(data.out.dir, "PCAWG.pearson.RDS") )
# 
# #Plotting individual experiments
# # pp = plot_all_experiments(PCAWG.pearson$Biliary_AdenoCA, "title")
# 
# PCAWG.pearson.summary.nets = lapply(PCAWG.pearson, summarize_calcs )
# 
# saveRDS(PCAWG.pearson.summary.nets,
#         file = file.path(data.out.dir, "PCAWG.pearson.summary.nets.RDS") )
# 
# 
# ### Spearman correlation
# 
# PCAWG.spearman = calculate_bootstrap_metric(PCAWG.full.subset.ann, cor_sigs,
#                                  min.tissue.samples = 20,
#                                  sample.rate = 0.9, sample.N = 100,
#                                  N = 1, seed = 1, p.adjust = TRUE,
#                                  method = "spearman")
# 
# saveRDS(PCAWG.spearman,
#         file = file.path(data.out.dir, "PCAWG.spearman.RDS") )
# 
# #Plotting individual experiments
# # pp = plot_all_experiments(PCAWG.spearman$Biliary_AdenoCA, "title")
# 
# PCAWG.spearman.summary.nets = lapply(PCAWG.spearman, summarize_calcs )
# 
# saveRDS(PCAWG.spearman.summary.nets,
#         file = file.path(data.out.dir, "PCAWG.spearman.summary.nets.RDS") )
# 
# 
# ### CoDa correlation
# 
# PCAWG.coda = calculate_bootstrap_metric(PCAWG.full.subset.ann, cor_coda,
#                                   min.tissue.samples = 20,
#                                   sample.rate = 0.9, sample.N = 100,
#                                   N = 1, seed = 1, p.adjust = TRUE,
#                                   rand.add = FALSE)
# 
# saveRDS(PCAWG.coda,
#         file = file.path(data.out.dir, "PCAWG.coda.RDS") )
# 
# #Plotting individual experiments
# # pp = plot_all_experiments(PCAWG.coda$Biliary_AdenoCA, "title")
# 
# PCAWG.coda.summary.nets = lapply(PCAWG.coda, summarize_calcs )
# 
# saveRDS(PCAWG.coda.summary.nets,
#         file = file.path(data.out.dir, "PCAWG.coda.summary.nets.RDS") )


### Writing to supplements
PCAWG.cooccurrence.summary.nets = readRDS(file = file.path(data.out.dir, "PCAWG.cooccurrence.summary.nets.RDS") )

PCAWG.bcmi.summary.nets = readRDS(file = file.path(data.out.dir, "PCAWG.bcmi.summary.nets.RDS") )

PCAWG.spearman.summary.nets = readRDS(file = file.path(data.out.dir, "PCAWG.spearman.summary.nets.RDS") )

PCAWG.coda.summary.nets = readRDS(file = file.path(data.out.dir, "PCAWG.coda.summary.nets.RDS") )


### Removing the supplementary file if it exists
out.file = here("supp_data/pcawg_sig_sig_all_metrics_report.xlsx")

if (file.exists(out.file)) {
    file.remove(out.file)  
} 

all.interactions = list(MI = PCAWG.bcmi.summary.nets,
                        cooccurrence = PCAWG.cooccurrence.summary.nets,
                        CoDa = PCAWG.coda.summary.nets,
                        Spearman = PCAWG.spearman.summary.nets)


sig.threshold = 0.05

all.sig.interactions = all.interactions
for(tissue in tissues) {
    
    for (metric.name in names(all.interactions)) {
        cat(tissue, ":", metric.name, "\n")
        out =  all.interactions[[metric.name]][[tissue]]
        metric.values = out %>% rm_zeros()
        if (length(metric.values) == 0) {next} 
        
        p.values = matrix(1, nrow = nrow(metric.values), 
                          ncol = ncol(metric.values), 
                          dimnames = (dimnames(metric.values)))
        
        metric.values[lower.tri(metric.values, diag = TRUE) ] = 0
        
        non.zero.values = which(metric.values != 0, arr.ind = TRUE)
        
        for(i in 1:nrow(non.zero.values)) {
            sig1 = rownames(metric.values)[non.zero.values[i, 1]]
            sig2 = rownames(metric.values)[non.zero.values[i, 2]]
            
            p.values[sig1, sig2] = p.values[sig1, sig2] = min(percentile(PCAWG.null.dists[[tissue]][[metric.name]][sig1, sig2, ],
                                                                         abs(metric.values[sig1, sig2])),
                                                              1 - percentile(PCAWG.null.dists[[tissue]][[metric.name]][sig1, sig2, ],
                                                                             abs(metric.values[sig1, sig2]))  )
            
            if(p.values[sig1, sig2] > sig.threshold) {
                out[sig1, sig2] = 0
                out[sig2, sig1] = 0
            }
        }
        all.sig.interactions[[metric.name]][[tissue]] = out
    }
}


### Write the metric name
for (tissue in tissues) {
    for (metric.name in names(all.sig.interactions)) {
        
        cat(tissue, "\n")
        metric.values = all.sig.interactions[[metric.name]][[tissue]] %>% rm_zeros()
        
        if (length(metric.values) == 0) {next} 
        
        p.values = matrix(NA, nrow = nrow(metric.values), 
                          ncol = ncol(metric.values), 
                          dimnames = (dimnames(metric.values)))
        
        metric.values[lower.tri(metric.values, diag = TRUE) ] = NA
        
        non.zero.values = which(metric.values != 0, arr.ind = TRUE)
        
        for(i in 1:nrow(non.zero.values)) {
            sig1 = rownames(metric.values)[non.zero.values[i, 1]]
            sig2 = rownames(metric.values)[non.zero.values[i, 2]]
            
            p.values[sig1, sig2] = min(percentile(PCAWG.null.dists[[tissue]][[metric.name]][sig1, sig2, ],
                                                  abs(metric.values[sig1, sig2])),
                                       1 - percentile(PCAWG.null.dists[[tissue]][[metric.name]][sig1, sig2, ],
                                                      abs(metric.values[sig1, sig2]))  )
        }
        
        if (!file.exists(out.file)) {
            wb <- createWorkbook()
            cat("File doesn't exist.\n")
        } else { 
            wb = loadWorkbook(out.file)
        }
        
        if (!(tissue %in% names(wb))) {
            cat("Sheet doesn't exist.\n")
            addWorksheet(wb, tissue)
            writeData(wb = wb, 
                      sheet = tissue,
                      x = metric.name)
            start.row = 2
        } else {
            cat("Sheet exists.\n")
            start.row = nrow(readWorkbook(out.file, 
                                          sheet = tissue, 
                                          colNames = FALSE, 
                                          skipEmptyRows = FALSE))+5
        }
        
        writeData(wb = wb, 
                  sheet = tissue,
                  x = metric.name,
                  colNames = TRUE,
                  startRow = start.row - 1,
                  startCol = 1) 
        cat("File exists.\n")
        
        
        cat(metric.name, "\t", start.row, "\n") 
        
        ### Write the metric values
        writeData(wb = wb, 
                  sheet = tissue,
                  x = metric.values,
                  colNames = TRUE,
                  rowNames = TRUE,
                  startRow = start.row)
        
        ### Write the p-values annotation
        writeData(wb = wb, 
                  sheet = tissue,
                  x = "p-values",
                  colNames = TRUE,
                  startRow = start.row - 1,
                  startCol = ncol(metric.values)+4) 
        
        ### Write the p-values matrix
        
        p.values [upper.tri(p.values)] =  ## replace NA's with 1 in upper.tri
            pmin(1, p.values[upper.tri(p.values)], na.rm = TRUE)
        
        writeData(wb = wb, 
                  sheet = tissue,
                  x = p.values,
                  colNames = TRUE,
                  rowNames = TRUE,
                  startRow = start.row,
                  startCol = ncol(metric.values)+4)
        
        saveWorkbook(wb, file = out.file, overwrite = TRUE)
        
    }
}

saveRDS(all.sig.interactions, file = here("data/RDS/PCAWG/signatures/PCAWG.all.significant.interactions.RDS"))
