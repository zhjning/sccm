% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_relationships.R
\name{seekGroup}
\alias{seekGroup}
\title{seekGroup}
\usage{
seekGroup(var, res, group)
}
\arguments{
\item{var}{a character of a component name in result table}

\item{res}{a data frame contains the treated summary table of the mediation model}

\item{group}{a vector contains names of the components used in the mediation model}
}
\value{
the belonged group of the varaible: "X", "Y", "M1", "M2", "id,cell type"
}
\description{
wrap function for seek group
}
