library(here)

source(here("R/load_packages.R"))

library(readxl)

pathway_dir = here("data/raw/TCGA/10_onco_pathways/")

pathways_file = file.path(pathway_dir, 
                          "1-s2.0-S0092867418303593-mmc4.xlsx")

pathway.mutations = as.data.frame(readxl::read_excel(pathways_file, sheet = "Pathway level", na = "NA"))

rownames(pathway.mutations) = pathway.mutations$SAMPLE_BARCODE

TCGA.full.subset.ann = readRDS(here("data", "RDS", "TCGA", "signatures",
                                    "TCGA.full.subset.ann.RDS"))

signatures.sample.names = TCGA.full.subset.ann$Sample.Names

signatures.sample.names = sapply(signatures.sample.names,
                                 function(x) substr(x, 0, 15))

rownames(TCGA.full.subset.ann) = signatures.sample.names

common.samples = intersect(pathway.mutations$SAMPLE_BARCODE, signatures.sample.names)

pathway.mutations = pathway.mutations[common.samples, ]
TCGA.full.subset.ann.pathways = TCGA.full.subset.ann[common.samples, ]


pathway.mutations = pathway.mutations %>% 
    tibble::add_column(donor_id = rownames(pathway.mutations), 
                                         .after = "SAMPLE_BARCODE") %>% 
    tibble::add_column(Cancer.Types = TCGA.full.subset.ann.pathways$Cancer.Types,
                       .after = "donor_id")


saveRDS(pathway.mutations, file = here("data/RDS/TCGA/10_onco_pathways",
                                              "tcga_pathways.RDS"))

saveRDS(TCGA.full.subset.ann.pathways, file = here("data/RDS/TCGA/10_onco_pathways",
                                                    "TCGA.full.subset.ann.pathways.RDS"))

