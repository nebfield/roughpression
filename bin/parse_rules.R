#!/usr/bin/env Rscript

library("dplyr")
library("tibble")

set.seed(0451)

# 1: phyloseq object
# 2: rule-support file
# 3: discretised data

args <- commandArgs(trailingOnly=TRUE)

ps <- readRDS(args[[1]])
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

write.table(interesting_rules %>% select(Decision, Rule, Support),
            file = "annotated_rules.txt",
            quote = FALSE,
            sep = "\t", 
            row.names = FALSE)

bugs <-
  stringr::str_extract_all(interesting_rules$Rule,
                           "\\w+;\\w+;\\w+;\\w+;\\w+;\\w+;\\w+")
bugs <- sapply(bugs, function(x) stringr::str_replace_all(x, ";", "."))
disc_trim <- disc[, c(unlist(bugs), "Class")]

