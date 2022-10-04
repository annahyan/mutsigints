library(here)

PCAWG.int.tissues = list()

PCAWG.int.tissues$pos.tissues = read.table(
        file = here("supp_data", "PCAWG_positive_interaction_tissue_summaries.tsv"),
        row.names = 1, h = TRUE, sep = "\t", na.strings = NA)

PCAWG.int.tissues$neg.tissues = read.delim(
    file = here("supp_data", "PCAWG_negative_interaction_tissue_summaries.tsv"),
    row.names = 1, h = TRUE, sep = "\t", na.strings = NA)


TCGA.int.tissues = list()

TCGA.int.tissues$pos.tissues = read.delim(
    file = here("supp_data", "TCGA_positive_interaction_tissue_summaries.tsv"),
    row.names = 1, h = TRUE, sep = "\t", na.strings = "")

TCGA.int.tissues$neg.tissues = read.delim(
    file = here("supp_data", "TCGA_negative_interaction_tissue_summaries.tsv"),
    row.names = 1, h = TRUE, sep = "\t", na.strings = "")


compare_PCAWG_TCGA_tissues = function(pcawg.tissues, tcga.tissues) {
    
    total.ints = 0
    common.sigs = intersect(colnames(pcawg.tissues), colnames(tcga.tissues))
    
    intersect.ints = matrix(rep("", length(common.sigs) ** 2),nrow = length(common.sigs),
                            dimnames = list(common.sigs, common.sigs) )
    
    for (i in 1:(length(common.sigs) - 1 ) ) {
        for (j in (i + 1) : length(common.sigs)) {
            i.sig = common.sigs[i]
            j.sig = common.sigs[j]
            PCAWG.tiss = tolower(strsplit(pcawg.tissues[i.sig, j.sig], ", ")[[1]])
            TCGA.tiss = gsub("breast_cancer", "breast_adenoca", 
                             tolower(strsplit(tcga.tissues[i.sig, j.sig], ", ")[[1]]))
            common.tiss = intersect(PCAWG.tiss, TCGA.tiss)
            intersect.ints[i.sig, j.sig] = paste0(common.tiss,
                                                              collapse = ", ")
            total.ints = total.ints + length(common.tiss)
        }
    }
    non.empty.rows = which(rowSums(intersect.ints > "" ) > 0)
    non.empty.cols = which(colSums(intersect.ints > "" ) > 0)
    
    intersect.ints = intersect.ints[non.empty.rows, non.empty.cols]
    
    return(list(common.tissues = intersect.ints, total.counts = total.ints))
}

pos.ints = compare_PCAWG_TCGA_tissues(PCAWG.int.tissues$pos.tissues,
                                      TCGA.int.tissues$pos.tissues)

write.table(as.data.frame(pos.ints), 
            file = here("supp_data", "common_positive_interaction_tissue_summaries.tsv"),
            row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

neg.ints = compare_PCAWG_TCGA_tissues(PCAWG.int.tissues$neg.tissues,
                                      TCGA.int.tissues$neg.tissues)

write.table(as.data.frame(neg.ints), 
            file = here("supp_data", "common_negative_interaction_tissue_summaries.tsv"),
            row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)