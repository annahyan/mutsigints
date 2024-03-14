# library(chrwiseSignatures)
library(ggraph)
library(lsa)
library(gridExtra)
library(cowplot)
library(ggstar)
library(viridis)
library(pheatmap)

library(here)
source(here("R/load_packages.R"))
library(ggrepel)
library(ggpmisc)

fig.dir = here("figures/supp_figures")

if (!file.exists(fig.dir)) {
    dir.create(fig.dir)
}

source(here("R/process_sigs_TCGA.R"))
source(here("R/process_sigs_PCAWG.R"))

tissues = unique(PCAWG.full.subset$Cancer.Types ) 

sample_counts = sapply(tissues, function(x) 
    subset_tissue(PCAWG.full.subset, x) %>% nrow())

sig_counts = sapply(tissues, function(x) 
    subset_tissue(PCAWG.full.subset, x) %>% ncol())

sig_sample_count_data = data.frame(
    cancer.types = tissues,
    sample_counts = sample_counts,
    sig_counts = sig_counts) 

set.seed(5)
sig_sample_cor = sig_sample_count_data %>% 
    ggplot(aes(x = sample_counts, y = sig_counts)) + 
    geom_point() + 
    geom_smooth(method='lm', formula= y~x, color = rgb(0,0,0, alpha = .4), fill = "darkolivegreen3") +
    geom_text_repel(data = subset(sig_sample_count_data, sig_counts > 10 | sample_counts > 100) %>% 
                        sample_n(size = 10),
                    aes(x = sample_counts, y = sig_counts, label = cancer.types),
                    size = 3) +
    stat_fit_glance(method = 'lm',
                    method.args = list(formula = y~x),
                    geom = 'text',
                    aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                    label.x = 280, label.y = 22, size = 4, color = "darkred") +
    theme_classic(base_size = 13) +
    xlab("Sample counts") + ylab("Signature counts")
sig_sample_cor

ggsave(plot = sig_sample_cor, 
       filename = file.path(fig.dir, "sig_sample_correlations.pdf"),
       width = 4, height = 3)
