library(here)

### Loading commonly used packages

source(here("R/load_packages.R"))


PCAWG.sample.sheet = read.delim(
    here("data/raw/PCAWG/metadata",
    "pcawg_sample_sheet.tsv" ) )

PCAWG.signatures.full = read.delim(
    here("data/raw/PCAWG/signatures/SigProfiler",
              "PCAWG_sigProfiler_SBS_signatures_in_samples.csv"),
    sep = ",")


# filtering raw PCAWG data ------------------------------------------------

PCAWG.summaries = data.frame(
    counts = rowSums(PCAWG.signatures.full %>% select(4:ncol(.))), 
                      accuracy = PCAWG.signatures.full$Accuracy)


PCAWG.filter_rec_acc = PCAWG.signatures.full %>% filter(Accuracy > ACCURACY, 
                                 PCAWG.summaries$counts > MUT.MIN.COUNTS) 

selected_tissues = PCAWG.filter_rec_acc %>% 
    count(Cancer.Types, sort = TRUE) %>% 
    filter(n > MIN.TISSUE.SIZE) %>% 
    pull(Cancer.Types)

selected_signatures = PCAWG.filter_rec_acc %>%
    select(4:ncol(.)) %>% 
    select_if(colSums(. > 0) > SIG.MIN.COUNTS) %>% 
    colnames()

PCAWG.full.subset = PCAWG.filter_rec_acc %>% filter(Cancer.Types %in% selected_tissues) %>% 
    select(c("Cancer.Types", "Sample.Names", "Accuracy", all_of(selected_signatures))) %>% 
    mutate( 
        Sample.Names = PCAWG.sample.sheet[ ## using donor IDs as sample names
            match(Sample.Names,
              PCAWG.sample.sheet$icgc_specimen_id) , "icgc_donor_id"]) %>% 
    mutate(Cancer.Types = gsub("-", "_", Cancer.Types))

saveRDS(PCAWG.full.subset, 
        file = here("data/RDS", "PCAWG/signatures/PCAWG.full.subset.RDS") )

# Signatures annotated ----------------------------------------------------

signature.annotations = read.delim(here("data/raw", "signature_annotations.tsv"))

PCAWG.sigs.annotated = set_signature_annotation(annotations = signature.annotations, 
                                                signature.matrix.df =
                                                    PCAWG.full.subset %>% select(4:ncol(.))  %>% as.data.frame) 

PCAWG.full.subset.ann = cbind(PCAWG.full.subset[,1:3], PCAWG.sigs.annotated)

saveRDS(PCAWG.full.subset.ann, 
        file = here("data/RDS", "PCAWG/signatures/PCAWG.full.subset.ann.RDS"))   

