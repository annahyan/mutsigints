sig1 = "PolE"
sig2 = "WNT"
tissue = "Uterus_AdenoCa"
tissues = tissue
signatures = c(sig1, sig2)

dataset = TCGA.path.sig.aggregate
clin.df = TCGA.clin.df
param.values = list(
    "age.at.diagnosis" = c(TRUE, FALSE),
    "with.total.muts" = c(TRUE, FALSE), 
    "tmb.logged" = c(TRUE),
    "binary.status" = c(FALSE),
    "epistatic" = c(TRUE, FALSE))
filename = here("supp_data/tcga_positive_sig_path_survivals.xlsx")
min.sample.fraction = 0
rm.non.sig.sheets = TRUE
return.only.sig = TRUE



oo = pick_survival_model_int(dataset,
                                   tissue,
                                   signatures,
                                   clin.df,
                                   param.values,
                                   min.sample.fraction = min.sample.fraction,
                                   filename = filename,
                                   rm.non.sig.sheets = FALSE,
                                   return.only.sig = FALSE)

age.at.diagnosis = TRUE
with.total.muts = TRUE
tmb.logged = TRUE
binary.status = FALSE
epistatic = TRUE

