% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Stage1.R
\name{.hybrid_test}
\alias{.hybrid_test}
\title{Hybrid test for eNODAL object.}
\usage{
.hybrid_test(
  eNODAL_obj,
  index,
  type = "interaction",
  test_func = "lm",
  Boot = 50
)
}
\arguments{
\item{eNODAL_obj}{An eNODAL object as output.}

\item{index}{Number of which omics feature to be tested.}

\item{type}{Testing type. Can be chosen from
"linear"(linear vs gam using npntest), "interaction"(interaction effect)
"meta1"(marginal effect of Meta1) and "meta2"(marginal effect of Meta2)}

\item{test_func}{method for fitting the model, can be chosen from "gam" and "lm".}

\item{Boot}{Number of bootstrapping. By default is 50.}
}
\value{
pvalue number.
}
\description{
Hybrid test for eNODAL object.
}
