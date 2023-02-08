library(here)

### Loading commonly used packages

source(here("R/load_packages.R"))


TCGA.signatures.full = read.delim(
    here("data/raw/TCGA/signatures/SigProfiler",
              "TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv"),
    sep = ",")


# filtering raw TCGA data ------------------------------------------------

TCGA.summaries = data.frame(
    counts = rowSums(TCGA.signatures.full %>% select(4:ncol(.))), 
                      accuracy = TCGA.signatures.full$Accuracy)


TCGA.filter_rec_acc = TCGA.signatures.full %>% filter(Accuracy > ACCURACY, 
                                 TCGA.summaries$counts > MUT.MIN.COUNTS) 

selected_tissues = TCGA.filter_rec_acc %>% 
    count(Cancer.Types, sort = TRUE) %>% 
    filter(n > MIN.TISSUE.SIZE) %>% 
    pull(Cancer.Types)

selected_signatures = TCGA.filter_rec_acc %>%
    select(4:ncol(.)) %>% 
    select_if(colSums(. > 0) > SIG.MIN.COUNTS) %>% 
    colnames()

TCGA.full.subset = TCGA.filter_rec_acc %>% filter(Cancer.Types %in% selected_tissues) %>% 
    select(c("Cancer.Types", "Sample.Names", "Accuracy", all_of(selected_signatures))) %>% 
    mutate(Cancer.Types = gsub("-", "_", Cancer.Types))

saveRDS(TCGA.full.subset, 
        file = here("data/RDS", "TCGA/signatures/TCGA.full.subset.RDS") )

# Signatures annotated ----------------------------------------------------

signature_annotations = read.delim(here("data/raw", "signature_annotations.tsv"))

TCGA.sigs.annotated = set_signature_annotation(
    annotations = signature_annotations,
    signature.matrix.df = 
        TCGA.full.subset %>% select(4:ncol(.))  %>% as.data.frame) 

TCGA.full.subset.ann = cbind(TCGA.full.subset[,1:3], TCGA.sigs.annotated)

saveRDS(TCGA.full.subset.ann, 
        file = here("data/RDS", "TCGA/signatures/TCGA.full.subset.ann.RDS"))   

