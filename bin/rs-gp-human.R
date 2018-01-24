data("GlobalPatterns")
set.seed(0451)

GPr  <- transform_sample_counts(GlobalPatterns, function(x) x / sum(x) )
GPfr <- filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)
GP_human <-
  phyloseq::subset_samples(GPfr, SampleType == "Feces" |
                             SampleType == "Skin" |
                             SampleType == "Tongue")
GP_human <- t(GP_human)
GP_y <- phyloseq::sample_data(GP_human)$SampleType
otu <- data.frame(phyloseq::otu_table(GP_human))
ra_annot <- cbind(data.frame(phyloseq::otu_table(GP_human), Class = GP_y, check.names = FALSE))

# make a decision table
dt <-
  RoughSets::SF.asDecisionTable(
    ra_annot,
    decision.attr = ncol(ra_annot),
    indx.nominal = ncol(ra_annot)
  )

# discretise 
cuts <-
  RoughSets::D.discretization.RST(dt, type.method = "unsupervised.quantiles", nOfIntervals = 3)
# warning about some features that can't be discretised (these will be removed
# by the quickreduct)
discretised <- suppressWarnings(RoughSets::SF.applyDecTable(dt, cuts))

# remove superfluous bugs by generating a single superreduct
fs <- RoughSets::FS.quickreduct.RST(discretised)
discretised_reduced <- SF.applyDecTable(discretised, fs)
# only one feature retained in this superreduct!

# create a rough set using the feature retained by the superreduct 
IND.A <- BC.IND.relation.RST(discretised_reduced, feature.set = 1)
roughset <- BC.LU.approximation.RST(discretised_reduced, IND.A)

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

accuracy(roughset, "Feces")
quality(roughset, "Feces", sum(discretised$Class == "Feces"))
accuracy(roughset, "Skin")
quality(roughset, "Skin", sum(discretised$Class == "Skin"))
accuracy(roughset, "Tongue")
quality(roughset, "Tongue", sum(discretised$Class == "Tongue"))

# generate some rules
?RoughSets::RI.indiscernibilityBasedRules.RST()
rules <- RoughSets::RI.indiscernibilityBasedRules.RST(discretised, fs)
# look at the rules
as.character(rules)

# what's the bug we're talking about?
interesting_bug <- stringr::str_sub(names(fs$reduct), start = 2)
print(phyloseq::tax_table(GP_human)[interesting_bug, ])
# Methanobrevibacter smithii lives in the human gut and shouldn't be found on the skin
# https://en.wikipedia.org/wiki/Methanobrevibacter_smithii


