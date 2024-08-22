library(here)

rm(list = ls())

source(here("R/load_packages.R"))
library(openxlsx)

fig.dir = here("figures/sig_sig_analysis")
data.out.dir = here("data/RDS/TCGA", "sig_sig_ints")

for (dir.path in c(fig.dir, data.out.dir) ) {
    if (! file.exists(dir.path)) {
        dir.create(dir.path)
    }
}

TCGA.full.subset.ann = readRDS(here("data", "RDS", "TCGA", "signatures",
                                     "TCGA.full.subset.ann.RDS"))

TCGA.null.dists = readRDS(file.path(here("data", "RDS", "TCGA", "signatures"),
                                      "TCGA.null.model.metrics.dist.RDS"))
tissues = names(TCGA.null.dists)
# ### Cooccurrence
# 
# TCGA.cooccurrence = calculate_bootstrap_metric(TCGA.full.subset.ann, cooccurrence,
#                  min.tissue.samples = 20,
#                  sample.rate = 0.9, sample.N = 100,
#                  N = 1, seed = 1, p.adjust = TRUE)
# 
# saveRDS(TCGA.cooccurrence,
#         file = file.path(data.out.dir, "TCGA.cooccurrence.RDS") )
# 
# # Plotting individual experiments
# # pp = plot_all_experiments(TCGA.cooccurrence$Biliary_AdenoCA, "title")
# 
# TCGA.cooccurrence.summary.nets = lapply(TCGA.cooccurrence, summarize_calcs )
# 
# saveRDS(TCGA.cooccurrence.summary.nets,
#         file = file.path(data.out.dir, "TCGA.cooccurrence.summary.nets.RDS") )
# 
# ### Mutual information(BCMI)
# 
# TCGA.bcmi = calculate_bootstrap_metric(TCGA.full.subset.ann, bcmi,
#                                       min.tissue.samples = 20,
#                                       sample.rate = 0.9, sample.N = 100,
#                                       N = 1, seed = 1, p.adjust = TRUE)
# 
# saveRDS(TCGA.bcmi,
#         file = file.path(data.out.dir, "TCGA.bcmi.RDS") )
# 
# #Plotting individual experiments
# # pp = plot_all_experiments(TCGA.bcmi$Biliary_AdenoCA, "title")
# 
# TCGA.bcmi.summary.nets = lapply(TCGA.bcmi, summarize_calcs )
# 
# saveRDS(TCGA.bcmi.summary.nets,
#         file = file.path(data.out.dir, "TCGA.bcmi.summary.nets.RDS") )
# 
# ### Pearson correlation
# 
# TCGA.pearson = calculate_bootstrap_metric(TCGA.full.subset.ann, cor_sigs,
#                               min.tissue.samples = 20,
#                               sample.rate = 0.9, sample.N = 100,
#                               N = 1, seed = 1, p.adjust = TRUE,
#                               method = "pearson")
# 
# saveRDS(TCGA.pearson,
#         file = file.path(data.out.dir, "TCGA.pearson.RDS") )
# 
# #Plotting individual experiments
# # pp = plot_all_experiments(TCGA.pearson$Biliary_AdenoCA, "title")
# 
# TCGA.pearson.summary.nets = lapply(TCGA.pearson, summarize_calcs )
# 
# saveRDS(TCGA.pearson.summary.nets,
#         file = file.path(data.out.dir, "TCGA.pearson.summary.nets.RDS") )
# 
# 
# ### Spearman correlation
# 
# TCGA.spearman = calculate_bootstrap_metric(TCGA.full.subset.ann, cor_sigs,
#                                  min.tissue.samples = 20,
#                                  sample.rate = 0.9, sample.N = 100,
#                                  N = 1, seed = 1, p.adjust = TRUE,
#                                  method = "spearman")
# 
# saveRDS(TCGA.spearman,
#         file = file.path(data.out.dir, "TCGA.spearman.RDS") )
# 
# #Plotting individual experiments
# # pp = plot_all_experiments(TCGA.spearman$Biliary_AdenoCA, "title")
# 
# TCGA.spearman.summary.nets = lapply(TCGA.spearman, summarize_calcs )
# 
# saveRDS(TCGA.spearman.summary.nets,
#         file = file.path(data.out.dir, "TCGA.spearman.summary.nets.RDS") )
# 
# 
# ### CoDa correlation
# 
# TCGA.coda = calculate_bootstrap_metric(TCGA.full.subset.ann, cor_coda,
#                                   min.tissue.samples = 20,
#                                   sample.rate = 0.9, sample.N = 100,
#                                   N = 1, seed = 1, p.adjust = TRUE,
#                                   rand.add = FALSE)
# 
# saveRDS(TCGA.coda,
#         file = file.path(data.out.dir, "TCGA.coda.RDS") )
# 
# #Plotting individual experiments
# # pp = plot_all_experiments(TCGA.coda$Biliary_AdenoCA, "title")
# 
# TCGA.coda.summary.nets = lapply(TCGA.coda, summarize_calcs )
# 
# saveRDS(TCGA.coda.summary.nets,
#         file = file.path(data.out.dir, "TCGA.coda.summary.nets.RDS") )


### Writing to supplements

TCGA.cooccurrence.summary.nets = readRDS(file = file.path(data.out.dir, "TCGA.cooccurrence.summary.nets.RDS") )

TCGA.bcmi.summary.nets = readRDS(file = file.path(data.out.dir, "TCGA.bcmi.summary.nets.RDS") )

TCGA.spearman.summary.nets = readRDS(file = file.path(data.out.dir, "TCGA.spearman.summary.nets.RDS") )

TCGA.coda.summary.nets = readRDS(file = file.path(data.out.dir, "TCGA.coda.summary.nets.RDS") )

out.file = here("supp_data/tcga_sig_sig_all_metrics_report.xlsx")

if(file.exists(out.fie)) {
    file.remove(out.file)    
}


all.interactions = list(MI = TCGA.bcmi.summary.nets,
                        cooccurrence = TCGA.cooccurrence.summary.nets,
                        CoDa = TCGA.coda.summary.nets,
                        Spearman = TCGA.spearman.summary.nets)

tissues = unique(TCGA.full.subset.ann$Cancer.Types)

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
            
            p.values[sig1, sig2] = min(percentile(TCGA.null.dists[[tissue]][[metric.name]][sig1, sig2, ],
                                                  abs(metric.values[sig1, sig2])),
                                       1 - percentile(TCGA.null.dists[[tissue]][[metric.name]][sig1, sig2, ],
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
            
            p.values[sig1, sig2] = min(percentile(TCGA.null.dists[[tissue]][[metric.name]][sig1, sig2, ],
                                                  abs(metric.values[sig1, sig2])),
                                       1 - percentile(TCGA.null.dists[[tissue]][[metric.name]][sig1, sig2, ],
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


saveRDS(all.sig.interactions, file = here("data/RDS/TCGA/signatures/TCGA.all.significant.interactions.RDS"))
