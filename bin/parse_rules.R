#!/usr/bin/env Rscript

library("dplyr")
library("tibble")

set.seed(0451)

# 1: phyloseq object
# 2: rule-support file
args <- commandArgs(trailingOnly=TRUE)

gut <- readRDS(args[[1]])
rules <- read.table(args[[2]], sep = ",", skip = 1,
                    col.names = c("Rule", "Support"), stringsAsFactors = FALSE)

# split string to get separate results column 
decs <-
  data.frame(do.call(rbind, stringr::str_split(rules$Rule, pattern = "->")),
             rules$Support,
             stringsAsFactors = FALSE)
colnames(decs) <- c("Rule", "Decision", "Support")

interesting_rules <- decs %>%
  group_by(Decision) %>%
  arrange(Decision) %>%
  filter(Support > 5) 

# extract DNA seqs from rules and match them to our taxonomy
seqs <- stringr::str_extract_all(interesting_rules$Rule, pattern = "[ACTG]+")
bugs <- sapply(seqs, function(x) phyloseq::prune_taxa(x, gut))
  
# extract taxonomic names from taxonomy table of phyloseq object into single string
# | Kingdom  | Phylum        | ... | Species  
# | Bacteria | Bacteroidetes | ... | NA      -> Bacteria;Bacteroidetes;...;NA
# collapse list of taxonomic names into single character vector separated by ---
# for writing out to a text file
get_taxonomy <- function(x) {
  paste(as.vector(unlist(
    apply(phyloseq::tax_table(x), 1, paste, collapse = ";")
  )), collapse = "---")
}

interesting_rules$bugs <- sapply(bugs, get_taxonomy)

write.table(interesting_rules %>% select(bugs, Decision, Rule, Support),
            file = "annotated_rules.txt",
            quote = FALSE,
            sep = "\t", 
            row.names = FALSE)
