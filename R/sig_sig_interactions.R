library(here)
fig.dir = here("figures/sig_sig_analysis")

if(! file.exists(fig.dir)){
    dir.create(fig.dir)
}

source(here("R/load_packages.R"))

PCAWG.full.subset.ann = readRDS(here("data", "RDS", "PCAWG", "signatures",
                                     "PCAWG.full.subset.ann.RDS"))

PCAWG.cooccurrence = get_metrics_list(PCAWG.full.subset.ann, cooccurrence,
                 min.tissue.samples = 20,
                 sample.rate = 0.9, sample.N = 100, 
                 N = 1, seed = 1, p.adjust = TRUE)
