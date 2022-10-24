library(RTCGA)
library(RTCGA.clinical)
library(survival)
library(survminer)
library(tidyverse)

# "survival_time", "vital_status", "age_at_diagnosis"

out.dir = here("data/RDS/TCGA/metadata")


TCGA.all.survivals = survivalTCGA(ACC.clinical, BLCA.clinical, BRCA.clinical,
                                  CESC.clinical, CHOL.clinical, COAD.clinical,
                                  COADREAD.clinical, DLBC.clinical, ESCA.clinical,
                                  FPPP.clinical, GBM.clinical, GBMLGG.clinical,
                                  HNSC.clinical, KICH.clinical, KIPAN.clinical,
                                  KIRC.clinical, KIRP.clinical, LAML.clinical,
                                  LGG.clinical, LIHC.clinical, LUAD.clinical,
                                  LUSC.clinical, MESO.clinical, OV.clinical, PAAD.clinical,
                                  PCPG.clinical, PRAD.clinical, READ.clinical,
                                  SARC.clinical, SKCM.clinical, STAD.clinical,
                                  STES.clinical, TGCT.clinical, THCA.clinical,
                                  THYM.clinical, UCEC.clinical, UCS.clinical,
                                  UVM.clinical,
                                  extract.cols=c("admin.disease_code", "patient.days_to_birth"))


TCGA.all.survivals$age_at_diagnosis = - as.numeric(TCGA.all.survivals$patient.days_to_birth) / 365

saveRDS(TCGA.all.survivals, file = file.path(out.dir, "TCGA.clin.df.RDS") )
