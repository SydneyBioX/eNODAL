## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(eNODAL)

## -----------------------------------------------------------------------------
data("Proteomics")

## -----------------------------------------------------------------------------
Z <- Proteomics$Z
Meta1 <- Proteomics$Meta1
Meta2 <- Proteomics$Meta2

## -----------------------------------------------------------------------------
eNODAL_obj = eNODAL_build(Z = Z, Meta1 = Meta1, Meta2 = Meta2, sig_test = "lm", test_func = "lm")

## -----------------------------------------------------------------------------
eNODAL_obj = runeNODAL(eNODAL_obj)

## -----------------------------------------------------------------------------
Cl_df <- get_clusters(eNODAL_obj)

## -----------------------------------------------------------------------------
show_clusters(eNODAL_obj)

## -----------------------------------------------------------------------------
ClusterPlot(eNODAL_obj)

## -----------------------------------------------------------------------------
NGFPlot(eNODAL_obj, "P19096", c("intake.P", "intake.C", "intake.F"))

## -----------------------------------------------------------------------------
NGFPlot(eNODAL_obj, "P19096", c("intake.P", "intake.C", "intake.F"), "Meta2")

## -----------------------------------------------------------------------------
NGFPlot(eNODAL_obj, "Meta1_C1", c("intake.P", "intake.C", "intake.F"))

## -----------------------------------------------------------------------------
Boxplot(eNODAL_obj, "Meta2", "P19096")

## -----------------------------------------------------------------------------
sessionInfo()

