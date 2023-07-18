### Common packages

library(openxlsx)
library(data.table)
# library(chrwiseSignatures)
library(tidyverse)
library(ggpubr)
library(tidygraph)
library(ggraph)
library(myggheatmap)
library(grid)
library(pheatmap)
library(viridis)

source(here("R/functions.R"))

# signature.annotations = read.delim(here("data/raw", "signature_annotations.tsv"))

### Some project invariants

MIN.TISSUE.SIZE = 25
ACCURACY = 0.85
MUT.MIN.COUNTS = 90
SIG.MIN.COUNTS = 0
