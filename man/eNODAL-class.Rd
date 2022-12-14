% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eNODAL_class.R
\docType{class}
\name{eNODAL-class}
\alias{eNODAL-class}
\alias{eNODAL}
\title{eNODAL class}
\description{
eNODAL class
}
\section{Slots}{

\describe{
\item{\code{Z}}{A dataframe omics data, rows are samples.}

\item{\code{Meta1}}{A dataframe/matrix/vector of nutrition intake, rows are samples.}

\item{\code{Meta2}}{A dataframe/matrix/vector of second meta data, rows are samples.}

\item{\code{Phenotype}}{A dataframe/matrix/vector of phenotype data, rows are samples.}

\item{\code{params}}{A list of parameters for hypothesis testing in eNODAL.}

\item{\code{eNODAL_middle}}{A list contain middle result of eNODAL.
Including possible formula used in GAM/Linear model($formula).}

\item{\code{eNODAL_output}}{A list eNODAL fitting result.}
}}

