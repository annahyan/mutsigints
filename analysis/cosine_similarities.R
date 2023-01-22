library(MutationalPatterns)
library(ggplot2)
library(tidyverse)

fig.dir = here("figures/supp_figures")

if (!file.exists(fig.dir)) {
    dir.create(fig.dir)
}

cosmic_v3 = get_known_signatures(muttype = "snv", source = "COSMIC")

cosmic_v3_similarities <- cos_sim_matrix(cosmic_v3, cosmic_v3)

cosmic_v3_simplot = plot_cosine_heatmap(cosmic_v3_similarities, 
                                        cluster_rows = FALSE, cluster_cols = FALSE, plot_values = TRUE)

ggsave(plot = cosmic_v3_simplot, 
       filename = file.path(fig.dir, "cosmic_v3_sim_heatmap.pdf"),
       width = 12, height = 12)

similarities_vector = cosmic_v3_similarities[ upper.tri(cosmic_v3_similarities, diag = FALSE)]

max_y = 6
similarity_histogram = similarities_vector %>% enframe(value = "cosine_similarity") %>% 
    ggplot(aes(x = cosine_similarity)) + 
    geom_histogram(aes(y = ..density..), 
                   color = "darkblue", fill = "white", size = 0.8) + 
    xlab("cosine similarity") + 
    stat_ecdf(aes_(y=bquote(..y.. * .(max_y))), size  =1.) + 
    scale_y_continuous(sec.axis=sec_axis(trans = ~100 * ./max_y, name="percentage")) +
    theme_bw(base_size = 18)


ggsave(plot = similarity_histogram, 
       filename = file.path(fig.dir, "cosmic_v3_sim_histogram.pdf"),
       width = 4, height = 3)
