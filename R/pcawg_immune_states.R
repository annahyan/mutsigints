library(here)

source(here("R/load_packages.R"))

library(readxl)

out.dir = here("data/RDS/PCAWG/immune_states")

if (!file.exists(out.dir)) {
    dir.create(out.dir)
}

immune.states = read.table(file = here("data/immune_classes/pcawg_predicted_classes.dat"))

con = file(here("data/expression/tophat_star_fpkm_uq.v2_aliquot_gl_filtered.tsv"), "r")
pcawg.RNAseq.sample.names = strsplit(readLines(con, n = 1), "\t")[[1]]
close(con)
pcawg.RNAseq.sample.names = pcawg.RNAseq.sample.names[2:length(pcawg.RNAseq.sample.names)]

immune.states = cbind(pcawg.RNAseq.sample.names, immune.states)

####

PCAWG.sample.sheet = get_sample_sheet(here("data/raw/PCAWG/metadata",
                                           "pcawg_sample_sheet.tsv" ) )

pcawg.specimen.hist = get_specimen_hist(here("data/raw/PCAWG/metadata", 
                                             "pcawg_specimen_histology_August2016_v9.xlsx"))

pcawg.immune.donorids = pcawg.specimen.hist [match(PCAWG.sample.sheet[ 
    match(pcawg.RNAseq.sample.names, PCAWG.sample.sheet$aliquot_id),] %>%  
        pull(icgc_specimen_id), pcawg.specimen.hist$`# icgc_specimen_id`),] %>% pull(icgc_donor_id)

immune.mapped.ids = cbind(pcawg.RNAseq.sample.names, pcawg.immune.donorids) %>% 
    na.omit() %>% as.data.frame()

#### Out

immune.states.mapping = cbind(immune.mapped.ids, immune.states[
    match(immune.mapped.ids$pcawg.RNAseq.sample.names, 
                    immune.states$pcawg.RNAseq.sample.names), c("RF", "DN")])

PCAWG.full.subset.ann = readRDS(here("data/RDS/PCAWG/signatures",
                                     "PCAWG.full.subset.ann.RDS"))

PCAWG.full.subset.ann.immune = PCAWG.full.subset.ann[ 
    match(immune.states.mapping$pcawg.immune.donorids, PCAWG.full.subset.ann$Sample.Names), ]

common.samples = intersect(immune.states.mapping$pcawg.immune.donorids,
                           PCAWG.full.subset.ann.immune$Sample.Names)

PCAWG.full.subset.ann.immune = PCAWG.full.subset.ann.immune[
    match(common.samples, PCAWG.full.subset.ann.immune$Sample.Names), ]
immune.states.mapping = immune.states.mapping[ 
    match(common.samples, immune.states.mapping$pcawg.immune.donorids), ]

rownames(PCAWG.full.subset.ann.immune) = common.samples
rownames(immune.states.mapping) = common.samples

immune.states.mapping$Cancer.Types = PCAWG.full.subset.ann.immune$Cancer.Types

colnames(immune.states.mapping)[1:2] = c("sample_id", "donor_id")
immune.states.mapping = immune.states.mapping[, c("sample_id", "donor_id", 
                                                  "Cancer.Types", "RF", "DN")]

## Adding xgboost classification output

xgboost.pcawg = readRDS(file = file.path(out.dir, "pcawg.xgboost.out.RDS"))

immune.states.mapping$XGboost = paste0("C", xgboost.pcawg[
    match(immune.states.mapping$sample_id, xgboost.pcawg$SampleIDs), "BestCall"])


saveRDS(immune.states.mapping, 
        file = file.path(out.dir, "pcawg_immune_states.RDS"))

saveRDS(PCAWG.full.subset.ann.immune, 
        file = file.path(out.dir, "PCAWG.full.subset.ann.immune.RDS"))
