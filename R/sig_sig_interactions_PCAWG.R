library(here)
source(here("R/load_packages.R"))

fig.dir = here("figures/sig_sig_analysis")
data.out.dir = here("data/RDS/TCGA", "sig_sig_ints")

for (dir.path in c(fig.dir, data.out.dir) ) {
    if (! file.exists(dir.path)) {
        dir.create(dir.path)
    }
}

PCAWG.full.subset.ann = readRDS(here("data", "RDS", "PCAWG", "signatures",
                                     "PCAWG.full.subset.ann.RDS"))

### Cooccurrence

PCAWG.cooccurrence = get_metrics_list(PCAWG.full.subset.ann, cooccurrence,
                 min.tissue.samples = 20,
                 sample.rate = 0.9, sample.N = 100, 
                 N = 1, seed = 1, p.adjust = TRUE)

saveRDS(PCAWG.cooccurrence, 
        file = file.path(data.out.dir, "PCAWG.cooccurrence.RDS") )

# Plotting individual experiments
# pp = plot_all_experiments(PCAWG.cooccurrence$Biliary_AdenoCA, "title")

PCAWG.cooccurrence.summary.nets = lapply(PCAWG.cooccurrence, summarize_calcs )

saveRDS(PCAWG.cooccurrence.summary.nets, 
        file = file.path(data.out.dir, "PCAWG.cooccurrence.summary.nets.RDS") )

### Mutual information(BCMI)

PCAWG.bcmi = get_metrics_list(PCAWG.full.subset.ann, bcmi,
                                      min.tissue.samples = 20,
                                      sample.rate = 0.9, sample.N = 100, 
                                      N = 1, seed = 1, p.adjust = TRUE)

saveRDS(PCAWG.bcmi, 
        file = file.path(data.out.dir, "PCAWG.bcmi.RDS") )

#Plotting individual experiments
# pp = plot_all_experiments(PCAWG.bcmi$Biliary_AdenoCA, "title")

PCAWG.bcmi.summary.nets = lapply(PCAWG.bcmi, summarize_calcs )

saveRDS(PCAWG.bcmi.summary.nets, 
        file = file.path(data.out.dir, "PCAWG.bcmi.summary.nets.RDS") )

### Pearson correlation

PCAWG.pearson = get_metrics_list(PCAWG.full.subset.ann, cor_sigs,
                              min.tissue.samples = 20,
                              sample.rate = 0.9, sample.N = 100, 
                              N = 1, seed = 1, p.adjust = TRUE,
                              method = "pearson")

saveRDS(PCAWG.pearson, 
        file = file.path(data.out.dir, "PCAWG.pearson.RDS") )

#Plotting individual experiments
# pp = plot_all_experiments(PCAWG.pearson$Biliary_AdenoCA, "title")

PCAWG.pearson.summary.nets = lapply(PCAWG.pearson, summarize_calcs )

saveRDS(PCAWG.pearson.summary.nets, 
        file = file.path(data.out.dir, "PCAWG.pearson.summary.nets.RDS") )


### Spearman correlation

PCAWG.spearman = get_metrics_list(PCAWG.full.subset.ann, cor_sigs,
                                 min.tissue.samples = 20,
                                 sample.rate = 0.9, sample.N = 100, 
                                 N = 1, seed = 1, p.adjust = TRUE,
                                 method = "spearman")

saveRDS(PCAWG.spearman, 
        file = file.path(data.out.dir, "PCAWG.spearman.RDS") )

#Plotting individual experiments
# pp = plot_all_experiments(PCAWG.spearman$Biliary_AdenoCA, "title")

PCAWG.spearman.summary.nets = lapply(PCAWG.spearman, summarize_calcs )

saveRDS(PCAWG.spearman.summary.nets, 
        file = file.path(data.out.dir, "PCAWG.spearman.summary.nets.RDS") )


### CoDa correlation

PCAWG.coda = get_metrics_list(PCAWG.full.subset.ann, cor_coda,
                                  min.tissue.samples = 20,
                                  sample.rate = 0.9, sample.N = 100, 
                                  N = 1, seed = 1, p.adjust = TRUE,
                                  rand.add = FALSE)

saveRDS(PCAWG.coda, 
        file = file.path(data.out.dir, "PCAWG.coda.RDS") )

#Plotting individual experiments
# pp = plot_all_experiments(PCAWG.coda$Biliary_AdenoCA, "title")

PCAWG.coda.summary.nets = lapply(PCAWG.coda, summarize_calcs )

saveRDS(PCAWG.coda.summary.nets, 
        file = file.path(data.out.dir, "PCAWG.coda.summary.nets.RDS") )

