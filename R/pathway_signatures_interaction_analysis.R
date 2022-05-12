library(here)
fig.dir = here("figures/pathways_analysis")

if(! file.exists(fig.dir)){
    dir.create(fig.dir)
}

source(here("R/load_packages.R"))

mutated.pathways.tissues = readRDS(file = here("data/RDS/PCAWG/10_onco_pathways",
                                               "pcawg_pathways.RDS"))

PCAWG.full.subset.ann.pathways = readRDS(file = here("data/RDS/PCAWG/10_onco_pathways",
                                                     "PCAWG.full.subset.ann.pathways.RDS"))


PATH_MIN_TISSUES = 30


skin.data = get_tissue_pathway_activities("Skin_Melanoma", 
                                          sigs.input = PCAWG.full.subset.ann.pathways,
                                          pathways.input = mutated.pathways.tissues)

skin.data$sigs.logged = skin.data$sigs %>% 
    mutate(across(.cols = everything(), ~ log(.x + 1 )))

skin.concat = merge(skin.data$sigs.logged, skin.data$paths, by = "row.names")

skin.rlm = MASS::rlm(skin.concat[, "Ageing"] ~ skin.concat[, "TGF-Beta"])

skin.rlm2 = MASS::rlm(skin.concat[, "TGF-Beta"] ~ skin.concat[, "Ageing"])


skin.lmrob = robustbase::lmrob(skin.concat[, "Ageing"] ~ skin.concat[, "TGF-Beta"])

skin.lmrob2 = robustbase::lmrob(skin.concat[, "TGF-Beta"] ~ skin.concat[, "Ageing"])




skin.rlm = MASS::rlm(skin.concat[, "Ageing"] ~ skin.concat[, "HIPPO"])
sfsmisc::f.robftest(skin.rlm, var = -1)

skin.rlm2 = MASS::rlm(skin.concat[, "HIPPO"] ~ skin.concat[, "Ageing"])
sfsmisc::f.robftest(skin.rlm2, var = -1)

skin.lmrob = robustbase::lmrob(skin.concat[, "Ageing"] ~ skin.concat[, "HIPPO"])

skin.lmrob0 = robustbase::glmrob(skin.concat[, "HIPPO"] ~ 1, family = binomial)
skin.lmrob2 = robustbase::glmrob(skin.concat[, "HIPPO"] ~ 1 + skin.concat[, "Ageing"], family = binomial)

anova(skin.lmrob0, skin.lmrob2)

colSums(skin.concat[, c("HIPPO", "Ageing")] > 0)


skin.concat[, c("HIPPO", "Ageing")] %>% rstatix::wilcox_test(Ageing ~ HIPPO)

skin.concat[, c("HIPPO", "Ageing")] %>% 
    mutate(HIPPO = factor(HIPPO)) %>% 
    ggplot(aes(x = HIPPO, y = Ageing, color = HIPPO, group = HIPPO)) + 
    geom_boxplot() + geom_jitter(color = "black")



rob.lin.mod = MASS::rlm(tissue.concat[, sig] ~ tissue.concat[, pathway])
int.mat[sig, pathway] = summary(rob.lin.mod)$coefficients[, "Value"][2]
p.values[sig, pathway] = tryCatch({
    sfsmisc::f.robftest(rob.lin.mod, var = -1)$p.value},
    error = function(e) {return(1)})


lin.mod = lm(tissue.concat[, sig] ~ tissue.concat[, pathway])
int.mat[sig, pathway] = summary(lin.mod)$coefficients[, "Estimate"][2]
p.values[sig, pathway] = summary(lin.mod)$coefficients[,"Pr(>|t|)"][2]



# Thinking of skin situation ----------------------------------------------

sig = "UV"
path = "RTK RAS"

int.cols = skin.concat[, c(sig, path)]

cont.table = table(as.data.frame(int.cols > 0))

ft = fisher.test(cont.table)


