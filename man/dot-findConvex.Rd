% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{.findConvex}
\alias{.findConvex}
\title{Finding convex hull of given data}
\usage{
.findConvex(x, y, rgnames, res = 101, x.limits = NA, y.limits = NA)
}
\arguments{
\item{x}{A vector of x-axis coordinates.}

\item{y}{A vector of y-axis coordinates.}

\item{rgnames}{Names of the grid.}

\item{res}{Resolution of the convex hull. By default is 100.}

\item{x.limits}{Limits of x-axis. By default is NA, will determined by data.}

\item{y.limits}{Limits of y-axis. By default is NA, will determined by data.}
}
\value{
A dataframe of data in convex hull.
}
\description{
Finding convex hull of given data
}
