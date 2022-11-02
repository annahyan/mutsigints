library(openxlsx)
library(rms)
library(survival)
library(survminer)
library(dplyr)
library(ranger)
library(ggfortify)
library(tidyverse)
library(here)

out.dir = here("data/RDS/PCAWG/metadata")

if (!file.exists(out.dir)) {
    dir.create(out.dir)
}

####
# Reading in donor-sample info --------------------------------------------
####

# load(here( "RDatas/PCAWG.sig.bundle.RData") )
PCAWG.full.subset.ann = readRDS(file = here("data/RDS", 
                                            "PCAWG/signatures/PCAWG.full.subset.ann.RDS") ) 

clin.data.dir = here('data/raw/PCAWG/metadata/')

PCAWG.donor.table = read.xlsx(file.path(clin.data.dir, "pcawg_donor_clinical_August2016_v9.xlsx") )
colnames(PCAWG.donor.table)[1] = gsub("#.", "", colnames(PCAWG.donor.table)[1])
PCAWG.sample.table = read.table(file.path(clin.data.dir, "pcawg_sample_sheet.tsv"), h = T, sep = "\t")

PCAWG.histology.table = read.xlsx(file.path(clin.data.dir, "pcawg_specimen_histology_August2016_v9.xlsx"))
colnames(PCAWG.histology.table)[1] = gsub("#.", "", colnames(PCAWG.histology.table)[1])

sig.donorids = as.character(PCAWG.full.subset.ann$Sample.Names)
# sig.donorids = as.character(PCAWG.sample.table[ 
#     match(sig.spids, PCAWG.sample.table$icgc_specimen_id), "icgc_donor_id"])


#### donor_survival_time and donor_interval_of_last_followup ####

PCAWG.sig.survivals = PCAWG.donor.table[ 
    match(sig.donorids, PCAWG.donor.table$icgc_donor_id), 
    c("donor_survival_time", "donor_interval_of_last_followup", 
      "donor_vital_status", "donor_age_at_diagnosis", "donor_sex")]

survivals_time = PCAWG.sig.survivals$donor_survival_time

# intervals_followup = PCAWG.sig.survivals$donor_interval_of_last_followup
# 
# survival.values = pairwise_nonNA(survivals_time, intervals_followup)
# 
# PCAWG.sig.survivals$survival_out = survival.values

colnames(PCAWG.sig.survivals) = gsub('donor_', '', colnames(PCAWG.sig.survivals))

PCAWG.ordered.histologies = PCAWG.histology.table[ 
    match(sig.donorids, PCAWG.histology.table$icgc_donor_id), ]

PCAWG.clin.df = cbind(PCAWG.sig.survivals, 
                      PCAWG.ordered.histologies)
PCAWG.clin.df$vital_status = as.numeric(PCAWG.clin.df$vital_status == "deceased")

saveRDS(PCAWG.clin.df, file = file.path(out.dir, "PCAWG.clin.df.RDS") )

