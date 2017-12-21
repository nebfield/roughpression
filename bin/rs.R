set.seed(0451)

library("RoughSets")
library("phyloseq")

oral <- readRDS("../cache/phyloseq_oral/oral.rds")

# samples to remove ------------------------------------------------------------
oral <-
  phyloseq::prune_samples(!(
    phyloseq::sample_names(oral) %in% c(
      "sample4452nd",
      "sample19742nd",
      "sampleOG2nd",
      "sampleDRY2nd"
    )
  ),
  oral)

oral <- phyloseq::transform_sample_counts(oral, function(x) x / sum(x))
oral_trimmed <- phyloseq::subset_samples(oral, SMOKING != "Occasionally")

# instances : rows, features: columns
oraldf <-
  data.frame(
    phyloseq::otu_table(oral_trimmed),
    SMOKING = phyloseq::sample_data(oral_trimmed)$SMOKING,
    gender = phyloseq::sample_data(oral_trimmed)$sex,
    cohort = phyloseq::sample_data(oral_trimmed)$cohort
  ) 

dt <-
  RoughSets::SF.asDecisionTable(
    dataset = oraldf,
    decision.attr = ncol(oraldf)
  ) # last column: decision attribute

reduct <- RoughSets::FS.feature.subset.computation(dt, method = "quickreduct.frst", randomize = TRUE)
dt_fs <- RoughSets::SF.applyDecTable(dt, reduct)

## evaluate index of objects
res.1 <-
  RoughSets::IS.FRIS.FRST(
    decision.table = dt_fs,
    list(
      threshold.tau = 0.75,
      alpha = 1,
      type.aggregation = c("t.tnorm", "lukasiewicz"),
      t.implicator = "lukasiewicz"
    )
  )
dt_is <- RoughSets::SF.applyDecTable(dt_fs, res.1)

# rule induction
control.ri <- list(type.aggregation = c("t.tnorm", "lukasiewicz"),
                   type.relation = c("tolerance", "eq.3"), 
                   t.implicator = "kleene_dienes")
decRules.hybrid <- RoughSets::RI.hybridFS.FRST(dt_is, control.ri)
