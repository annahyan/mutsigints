library(here)

source(here("R/load_packages.R"))

library(readxl)

pathway_dir = here("data/raw/PCAWG/10_onco_pathways/")

pathway_genes_file = file.path(pathway_dir, 
                               "genes-pathways.xlsx")

pathway_names = excel_sheets(pathway_genes_file)


pathway_genes = list()

for (pathname in pathway_names) {
    pathway_genes[[pathname]] = read_excel(pathway_genes_file, sheet = pathname)
}

driver_genes = read.delim(file.path(pathway_dir, 
                                    "TableS3_panorama_driver_mutations_ICGC_samples.controlled.tsv"))

sample_mut_pathways = function(sample_df) {
    
    pathway_status = setNames(rep(0, 10), pathway_names )
    
    pathway_mutations = sapply(pathway_genes, function(x) intersect(x$Gene, sample_df$gene) %>% length)
    pathway_mutations["Telomeric"] = sum(grepl("telomere", sample_df$gene))
    
    return(pathway_mutations)
}

mutated_pathways = sapply(
    split(driver_genes, f = driver_genes$sample_id), 
    FUN = sample_mut_pathways) %>% t() %>% as.data.frame() %>% 
    rownames_to_column("sample_id")


# Mapping sample_id to donor_id -------------------------------------------


PCAWG.sample.sheet = read.delim(
    here("data/raw/PCAWG/metadata",
         "pcawg_sample_sheet.tsv" ) )


pcawg_specimen_hist = read_excel(here("data/raw/PCAWG/metadata", 
                                      "pcawg_specimen_histology_August2016_v9.xlsx"))


pcawg_pathway_donorids = pcawg_specimen_hist [match(PCAWG.sample.sheet[ 
    match(mutated_pathways$sample_id, PCAWG.sample.sheet$aliquot_id),] %>%  
        pull(icgc_specimen_id), pcawg_specimen_hist$`# icgc_specimen_id`),] %>% pull(icgc_donor_id)


pcawg_path_mapping = data.frame(sample_id = mutated_pathways$sample_id)
pcawg_path_mapping = cbind(pcawg_path_mapping, 
                           icgc_specimen_id = PCAWG.sample.sheet[ 
                               match(pcawg_path_mapping$sample_id, PCAWG.sample.sheet$aliquot_id),] %>%  
                               pull(icgc_specimen_id))

pcawg_path_mapping = cbind( pcawg_path_mapping , 
                            pcawg_specimen_hist [
                                match(pcawg_path_mapping$icgc_specimen_id, 
                                      pcawg_specimen_hist$`# icgc_specimen_id`),] %>% 
                                select(icgc_donor_id, histology_abbreviation) )

mutated_pathways$donor_id = pcawg_pathway_donorids

mutated_pathways_tissues = mutated_pathways
mutated_pathways_tissues$Cancer.Types = pcawg_path_mapping[ 
    match(mutated_pathways$donor_id, pcawg_path_mapping$icgc_donor_id), ]$histology_abbreviation

new_col_order = c("sample_id", "donor_id", "Cancer.Types") # bringing these columns first
new_col_order = c(new_col_order, setdiff(colnames(mutated_pathways_tissues), new_col_order))
mutated_pathways_tissues = mutated_pathways_tissues[, new_col_order]



# merging signatures and pathways -----------------------------------------

PCAWG.full.subset.ann = readRDS(here("data/RDS/PCAWG/signatures",
                                     "PCAWG.full.subset.ann.RDS"))

PCAWG.full.subset.ann.pathways = PCAWG.full.subset.ann[ 
    match(mutated_pathways_tissues$donor_id, PCAWG.full.subset.ann$Sample.Names), ]

non_na_indeces = which(!is.na(PCAWG.full.subset.ann.pathways$Sample.Names ) )

PCAWG.full.subset.ann.pathways = PCAWG.full.subset.ann.pathways[non_na_indeces, ]
mutated_pathways_tissues = mutated_pathways_tissues[ non_na_indeces, ]

mutated_pathways_tissues = as.data.frame(mutated_pathways_tissues)

rownames(PCAWG.full.subset.ann.pathways) = PCAWG.full.subset.ann.pathways$Sample.Names
rownames(mutated_pathways_tissues) = PCAWG.full.subset.ann.pathways$Sample.Names

saveRDS(mutated_pathways_tissues, file = here("data/RDS/PCAWG/10_onco_pathways",
                                              "pcawg_pathways.RDS"))

saveRDS(PCAWG.full.subset.ann.pathways, file = here("data/RDS/PCAWG/10_onco_pathways",
                                              "PCAWG.full.subset.ann.pathways.RDS"))
