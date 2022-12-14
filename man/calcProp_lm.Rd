% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{calcProp_lm}
\alias{calcProp_lm}
\title{Calculate proportion of each component in linear model}
\usage{
calcProp_lm(eNODAL_obj, index, Cl_sel)
}
\arguments{
\item{eNODAL_obj}{An eNODAL_obj object}

\item{index}{A number specify which column to calculate}

\item{Cl_sel}{A text specify which cluster this variable belongs to}
}
\value{
A vector contain proportion of variance each component explained.
}
\description{
Calculate proportion of each component in linear model
}
