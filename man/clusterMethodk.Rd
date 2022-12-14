% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Stage2.R
\name{clusterMethodk}
\alias{clusterMethodk}
\title{Run clustering method (Need to specify number of clusters)}
\usage{
clusterMethodk(Dist, clmethod, k = 5, ...)
}
\arguments{
\item{Dist}{distance matrix to calculate}

\item{clmethod}{Method to run for clustering, can be chosen from:
hclust: Hierachical clustering,
apcluster: Affinity propagation clustering with k clusters,
kmeans: kmeans clustering,}

\item{k}{Number of clusters}

\item{...}{Other parameters can be passed to clusterMethodk.}
}
\value{
A list contain clustering result.
}
\description{
Run clustering method (Need to specify number of clusters)
}
