library(here)

source(here("R/load_packages.R"))

fig.dir = here("figures/sig_sig_analysis")
data.out.dir = here("data/RDS/TCGA", "sig_sig_ints")

for (dir.path in c(fig.dir, data.out.dir) ) {
    if (! file.exists(dir.path)) {
        dir.create(dir.path)
    }
}

TCGA.full.subset.ann = readRDS(here("data", "RDS", "TCGA", "signatures",
                                     "TCGA.full.subset.ann.RDS"))

### Cooccurrence

TCGA.cooccurrence = get_metrics_list(TCGA.full.subset.ann, cooccurrence,
                 min.tissue.samples = 20,
                 sample.rate = 0.9, sample.N = 100, 
                 N = 1, seed = 1, p.adjust = TRUE)

saveRDS(TCGA.cooccurrence, 
        file = file.path(data.out.dir, "TCGA.cooccurrence.RDS") )

# Plotting individual experiments
# pp = plot_all_experiments(TCGA.cooccurrence$Biliary_AdenoCA, "title")

TCGA.cooccurrence.summary.nets = lapply(TCGA.cooccurrence, summarize_calcs )

saveRDS(TCGA.cooccurrence.summary.nets, 
        file = file.path(data.out.dir, "TCGA.cooccurrence.summary.nets.RDS") )

### Mutual information(BCMI)

TCGA.bcmi = get_metrics_list(TCGA.full.subset.ann, bcmi,
                                      min.tissue.samples = 20,
                                      sample.rate = 0.9, sample.N = 100, 
                                      N = 1, seed = 1, p.adjust = TRUE)

saveRDS(TCGA.bcmi, 
        file = file.path(data.out.dir, "TCGA.bcmi.RDS") )

#Plotting individual experiments
# pp = plot_all_experiments(TCGA.bcmi$Biliary_AdenoCA, "title")

TCGA.bcmi.summary.nets = lapply(TCGA.bcmi, summarize_calcs )

saveRDS(TCGA.bcmi.summary.nets, 
        file = file.path(data.out.dir, "TCGA.bcmi.summary.nets.RDS") )

### Pearson correlation

TCGA.pearson = get_metrics_list(TCGA.full.subset.ann, cor_sigs,
                              min.tissue.samples = 20,
                              sample.rate = 0.9, sample.N = 100, 
                              N = 1, seed = 1, p.adjust = TRUE,
                              method = "pearson")

saveRDS(TCGA.pearson, 
        file = file.path(data.out.dir, "TCGA.pearson.RDS") )

#Plotting individual experiments
# pp = plot_all_experiments(TCGA.pearson$Biliary_AdenoCA, "title")

TCGA.pearson.summary.nets = lapply(TCGA.pearson, summarize_calcs )

saveRDS(TCGA.pearson.summary.nets, 
        file = file.path(data.out.dir, "TCGA.pearson.summary.nets.RDS") )


### Spearman correlation

TCGA.spearman = get_metrics_list(TCGA.full.subset.ann, cor_sigs,
                                 min.tissue.samples = 20,
                                 sample.rate = 0.9, sample.N = 100, 
                                 N = 1, seed = 1, p.adjust = TRUE,
                                 method = "spearman")

saveRDS(TCGA.spearman, 
        file = file.path(data.out.dir, "TCGA.spearman.RDS") )

#Plotting individual experiments
# pp = plot_all_experiments(TCGA.spearman$Biliary_AdenoCA, "title")

TCGA.spearman.summary.nets = lapply(TCGA.spearman, summarize_calcs )

saveRDS(TCGA.spearman.summary.nets, 
        file = file.path(data.out.dir, "TCGA.spearman.summary.nets.RDS") )


### CoDa correlation

TCGA.coda = get_metrics_list(TCGA.full.subset.ann, cor_coda,
                                  min.tissue.samples = 20,
                                  sample.rate = 0.9, sample.N = 100, 
                                  N = 1, seed = 1, p.adjust = TRUE,
                                  rand.add = FALSE)

saveRDS(TCGA.coda, 
        file = file.path(data.out.dir, "TCGA.coda.RDS") )

#Plotting individual experiments
# pp = plot_all_experiments(TCGA.coda$Biliary_AdenoCA, "title")

TCGA.coda.summary.nets = lapply(TCGA.coda, summarize_calcs )

saveRDS(TCGA.coda.summary.nets, 
        file = file.path(data.out.dir, "TCGA.coda.summary.nets.RDS") )

