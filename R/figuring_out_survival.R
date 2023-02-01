dataset = PCAWG.full.subset.ann
signatures = c("SBS40", "ROS")
tissues = "Eso_AdenoCA"
clin.df = PCAWG.clin.df
legend_pos = c(0.8, 0.8)
with.total.muts = TRUE
tmb.logged = TRUE
binary.status = FALSE
conf.int = TRUE

#### checking if TCGA 
TCGA = FALSE
if ( substr(dataset$Sample.Names[1], 1, 4) == "TCGA") {
    TCGA = TRUE
}

if (length(signatures) != 2) {
    stop(paste(c("2 signatures are required for this function:", signatures) ),
         collapse = " ")
    return(NULL)
}

tissues.subset = dataset %>% filter(Cancer.Types %in% tissues)

tissues.subset = tissues.subset[ ! duplicated(tissues.subset$Sample.Names), ]

if ( ! TCGA ) { ### for PCAWG
    adjusted.clin.df = clin.df[ match(tissues.subset$Sample.Names, 
                                      clin.df$icgc_donor_id), ]
    
    non_na.indeces = which(!is.na(adjusted.clin.df$icgc_specimen_id))
    
    tissues.subset = tissues.subset[non_na.indeces, ]
    adjusted.clin.df = adjusted.clin.df[ non_na.indeces, ]
    
} else { ### for TCGA
    
    tissues.subset$Sample.Names = substr(tissues.subset$Sample.Names, 1, 12)
    valid.indeces = which(tissues.subset$Sample.Names %in% clin.df$bcr_patient_barcode)
    tissues.subset = tissues.subset[valid.indeces, ]
    
    adjusted.clin.df = clin.df[ match(tissues.subset$Sample.Names, 
                                      clin.df$bcr_patient_barcode), ]
    
    adjusted.clin.df = as_tibble(adjusted.clin.df)
    
    adjusted.clin.df = rename(adjusted.clin.df, survival_time = times)
    adjusted.clin.df = rename(adjusted.clin.df, vital_status = patient.vital_status)
    
}

survival.df = cbind(tissues.subset[, 1:3], tissues.subset[, signatures], 
                    adjusted.clin.df[, c("survival_time", "vital_status", 
                                         "age_at_diagnosis"#, "sex"
                    ) ] )

survival.df = cbind(survival.df, total_muts = 
                        rowSums(tissues.subset[, 4:ncol(tissues.subset)]))





sbs1 = signatures[1]
sbs2 = signatures[2]

colnames(survival.df)[ which(colnames(survival.df) == sbs1)] = "SBS__1"
colnames(survival.df)[ which(colnames(survival.df) == sbs2)] = "SBS__2"

survival.df$exists__1 = as.numeric(survival.df$SBS__1 > 0)
survival.df$exists__2 = as.numeric(survival.df$SBS__2 > 0)
survival.df$exists__12 = as.numeric(survival.df$SBS__1 > 0 & survival.df$SBS__2 > 0)
survival.df$exists__None = as.numeric(survival.df$SBS__1 == 0 & survival.df$SBS__2 == 0)

sig.comb.status = apply(survival.df[, c("exists__1", "exists__2")], MARGIN = 1, 
                        function(x) {
                            if(x[1] == 1 & x[2] == 1) return(paste0(sbs1, "+", sbs2))
                            if(x[1] == 0 & x[2] == 1) return(paste0( sbs2))
                            if(x[1] == 1 & x[2] == 0) return(paste0( sbs1))
                            if(x[1] == 0 & x[2] == 0) return(paste0("None"))
                        } )

survival.df$status = sig.comb.status
survival.df$status = factor(survival.df$status, levels = c("None", sbs1,
                                                           sbs2, paste0(sbs1, "+", sbs2)))


if (binary.status) {
    survival.df$status = ifelse(survival.df$status == paste0(sbs1, "+", sbs2),
                                paste0(sbs1, "+", sbs2), "Others")
    survival.df$status = factor(survival.df$status, 
                                levels = c("Others", paste0(sbs1, "+", sbs2)))
    
}

survival.df$logTMB = log(survival.df$total_muts)


cox.m1 = coxph(Surv(survival_time, vital_status) ~ status, 
               data = survival.df, na.action = na.omit)

cox.m2 = coxph(Surv(survival_time, vital_status) ~ 
                   age_at_diagnosis + status, 
               data = survival.df, na.action = na.omit)

cox.m3 = coxph(Surv(survival_time, vital_status) ~ status + logTMB, 
               data = survival.df, na.action = na.omit)


cox.m4 = coxph(Surv(survival_time, vital_status) ~ age_at_diagnosis + 
                   status + logTMB, 
               data = survival.df, na.action = na.omit)


y = Surv(survival.df$survival_time, survival.df$vital_status)

x = survival.df[, c("age_at_diagnosis", "status", "logTMB")]

glmnet_fit <- glmnet(model.matrix(cox.m4), y, family = "cox")
coef(cox.m4)
