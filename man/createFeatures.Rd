% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{createFeatures}
\alias{createFeatures}
\title{Calculate interpretable features for eNODAL object}
\usage{
createFeatures(eNODAL_obj, baseline = NULL, cor_method = "spearman")
}
\arguments{
\item{eNODAL_obj}{An eNODAL object as input
Interaction effect only works for one numeric + one categorical variable}

\item{baseline}{Baseline for comparison(like "control" in Meta2). By default is NULL.}

\item{cor_method}{Method for calculate correlation matrix. By default is "spearman".}
}
\value{
eNODAL_obj object with created interpretable features
}
\description{
Calculate interpretable features for eNODAL object
}
