library(here)
library(tidygraph)
library(igraph)
library(ggraph)

source(here("R/functions.R"))

fig.dir = here("figures/pathways_analysis")

PCAWG.int.tissues = list()

PCAWG.int.tissues$pos.tissues = read.table(
        file = here("supp_data", "PCAWG_sig_pathway_positive_interaction_tissue_summaries.tsv"),
        row.names = 1, h = TRUE, sep = "\t", na.strings = NA)

PCAWG.int.tissues$neg.tissues = read.delim(
    file = here("supp_data", "PCAWG_sig_pathway_negative_interaction_tissue_summaries.tsv"),
    row.names = 1, h = TRUE, sep = "\t", na.strings = NA)


TCGA.int.tissues = list()

TCGA.int.tissues$pos.tissues = read.delim(
    file = here("supp_data", "TCGA_sig_pathway_positive_interaction_tissue_summaries.tsv"),
    row.names = 1, h = TRUE, sep = "\t", na.strings = NA)

TCGA.int.tissues$neg.tissues = read.delim(
    file = here("supp_data", "TCGA_sig_pathway_negative_interaction_tissue_summaries.tsv"),
    row.names = 1, h = TRUE, sep = "\t", na.strings = NA)


compare_PCAWG_TCGA_sig_pathways = function(pcawg.tissues, tcga.tissues) {
    
    total.ints = 0
    common.sigs = intersect(rownames(pcawg.tissues), rownames(tcga.tissues))
    
    common.paths = intersect(colnames(pcawg.tissues), colnames(tcga.tissues))
    
    intersect.ints = matrix(rep("", length(common.sigs) * length(common.paths)),
                            nrow = length(common.sigs),
                            ncol = length(common.paths),
                            dimnames = list(common.sigs, common.paths) )
    
    for (i in 1:length(common.sigs) ) {
        for (j in 1: length(common.paths)) {
            i.sig = common.sigs[i]
            j.path = common.paths[j]
            PCAWG.tiss = tolower(strsplit(pcawg.tissues[i.sig, j.path], ", ")[[1]])
            TCGA.tiss = gsub("breast_cancer", "breast_adenoca", 
                             tolower(strsplit(tcga.tissues[i.sig, j.path], ", ")[[1]]))
            common.tiss = intersect(PCAWG.tiss, TCGA.tiss)
            intersect.ints[i.sig, j.path] = paste0(common.tiss,
                                                              collapse = ", ")
            total.ints = total.ints + length(common.tiss)
        }
    }
    non.empty.rows = which(rowSums(intersect.ints > "" ) > 0)
    non.empty.cols = which(colSums(intersect.ints > "" ) > 0)
    
    intersect.ints = intersect.ints[non.empty.rows, non.empty.cols]
    
    return(list(common.tissues = intersect.ints, total.counts = total.ints))
}

pos.ints = compare_PCAWG_TCGA_sig_pathways(PCAWG.int.tissues$pos.tissues,
                                      TCGA.int.tissues$pos.tissues)

write.table(as.data.frame(order_matrix_rc(pos.ints$common.tissues)), 
            file = here("supp_data", "sig_pathway_common_positive_interaction_tissue_summaries.tsv"),
            row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

neg.ints = compare_PCAWG_TCGA_sig_pathways(PCAWG.int.tissues$neg.tissues,
                                      TCGA.int.tissues$neg.tissues)

write.table(as.data.frame(order_matrix_rc(neg.ints$common.tissues)), 
            file = here("supp_data", "sig_pathway_common_negative_interaction_tissue_summaries.tsv"),
            row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

#### 

make_matrix_square = function(mat) {
    
    rows = rownames(mat)
    cols = colnames(mat)
    
    all.elems = union(rows, cols)
    
    non.zeros = as.data.frame(which(mat != 0, arr.ind = TRUE))
    non.zeros[, "row"] = rows[non.zeros[, "row"]]
    non.zeros[, "col"] = cols[non.zeros[, "col"]]
    
    non.zeros = as.matrix(non.zeros)
    
    out.mat = matrix(rep(0, length(all.elems) ** 2 ), 
                    nrow = length(all.elems),
                    ncol = length(all.elems),
                    dimnames = list(all.elems, all.elems))
    
    out.mat[non.zeros] = mat[which(mat != 0)]
    
    out.mat = out.mat + t(out.mat)
    return(out.mat)
}

### PCAWG and TCGA common interactions network

all.sigs = union( 
    union(rownames(PCAWG.int.tissues$pos.tissues), rownames(TCGA.int.tissues$pos.tissues)),
    union(rownames(PCAWG.int.tissues$neg.tissues), rownames(TCGA.int.tissues$neg)))

all.paths = union( 
    union(colnames(PCAWG.int.tissues$pos.tissues), colnames(TCGA.int.tissues$pos.tissues)),
    union(colnames(PCAWG.int.tissues$neg.tissues), colnames(TCGA.int.tissues$neg)))




c.pos.ints = pos.ints$common.tissues
c.neg.ints = neg.ints$common.tissues

c.pos.ints[ c.pos.ints != "" ] = 1
c.pos.ints[ c.pos.ints == "" ] = 0
c.pos.ints = apply(c.pos.ints, 2, as.numeric)
dimnames(c.pos.ints) = dimnames(pos.ints$common.tissues)
c.pos.ints = make_matrix_square(c.pos.ints)


# c.neg.ints[ c.neg.ints != "" ] = 1
# c.neg.ints[ c.neg.ints == "" ] = 0
# c.neg.ints = apply(c.neg.ints, 2, as.numeric)
# dimnames(c.neg.ints) = dimnames(neg.ints$common.tissues)
# c.neg.ints = make_matrix_square(c.neg.ints)

pos.graph = as_tbl_graph(graph_from_adjacency_matrix(as.matrix(c.pos.ints)))
pos.graph = pos.graph %>%
    activate(edges) %>%
    mutate(type = "positive") %>%
    activate(nodes) %>%
    mutate(ifelse(name %in% all.sigs, "Signature", "Pathway"))

# neg.graph = as_tbl_graph(graph_from_adjacency_matrix(as.matrix(c.neg.ints)))
# neg.graph = neg.graph %>% 
#     activate(edges) %>% 
#     mutate(type = "negative") %>%
#     activate(nodes) %>% 
#     mutate(ifelse(name %in% all.sigs, "Signature", "Pathway"))

# joined.graph = graph_join(pos.graph, neg.graph)

# pp.graph = joined.graph %>% 
pp.graph = pos.graph %>% 
    ggraph(layout = "sugiyama") +
    geom_edge_link(aes(color = type), width = 0.7) + 
    geom_node_point() + 
    geom_node_label(aes(label = name)) + 
    scale_edge_color_manual(# limits = c("negative", "positive"),
                               values= c( rgb(210, 50, 60, maxColorValue = 255), 
                                          rgb(0, 140, 160, maxColorValue = 255))) + 
    theme_graph(base_family = "Helvetica") # + 
    # guides(edge_color = "none")

pp.graph$data[, c("x", "y")] = pp.graph$data[, c("y", "x")]

ggsave(file = file.path(fig.dir, "common_sig_pathway_interactions_PCAWG_TCGA.pdf"), 
       plot = minor_plot(pp.graph, 0.15), 
       width = 5, height = 4)
