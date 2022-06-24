### Common packages

library(openxlsx)
library(data.table)
library(chrwiseSignatures)
library(tidyverse)
library(ggpubr)
library(ggraph)
library(myggheatmap)

source(here("R/functions.R"))

### Some project invariants

MIN.TISSUE.SIZE = 25
ACCURACY = 0.85
MUT.MIN.COUNTS = 90
SIG.MIN.COUNTS = 0
