library("dplyr")
undisc_arff <- foreign::read.arff("~/projects/roughpression/pipeline/results/gut/gut.arff")
disc_arff <- foreign::read.arff("~/projects/roughpression/pipeline/results/gut/discretised.arff")
disc_table <- as.data.frame(disc_arff)

rules <- read.csv("~/projects/roughpression/pipeline/results/gut/rule-support.csv", header = FALSE, stringsAsFactors = FALSE)

# extract names of bacteria from generated rules 
extracted_bugs <- stringr::str_extract_all(rules$V1, "\\w+;\\w+;\\w+;\\w+\\w+;\\w+;\\w+;\\w+")
reduct <- unique(unlist(extracted_bugs))
disc_table_reduct <- disc_table[, reduct]

# add class labels (they get lost in the discretisation)
disc_table_reduct$class <- undisc_arff$Class
disc_table_reduct <- disc_table_reduct %>%
  filter(class != "Remission")
dt <- SF.asDecisionTable(disc_table_reduct, decision.attr = ncol(disc_table_reduct))

# create a rough set using the feature retained by the superreduct 
IND.A <- BC.IND.relation.RST(dt)
roughset <- BC.LU.approximation.RST(dt, IND.A)

accuracy <- function(roughset, class) {
  # the ratio of the size of the family of lower-approximation sets to the size
  # of the family of upper-approximation sets.
  
  lower <- getElement(roughset, "lower.approximation")
  lower_class <- getElement(lower, class)
  
  upper <- getElement(roughset, "upper.approximation")
  upper_class <- getElement(upper, class) # workers of the world unite!
  
  return(length(lower_class) / length(upper_class))
}

quality <- function(roughset, class, sum_objects) {
  # This represents the ratio of all objects classified with certainty to the
  # total number of objects in X.
  lower <- getElement(roughset, "lower.approximation")
  lower_class <- getElement(lower, class)
  
  return(length(lower_class) / sum_objects)
}

accuracy(roughset, "MDD")
quality(roughset, "MDD", sum(disc_table_reduct$class == "MDD"))
accuracy(roughset, "Control")
quality(roughset, "Control", sum(disc_table_reduct$class == "Control"))

sum(disc_table_reduct$class == "Control")
