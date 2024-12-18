% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_settings.R
\name{get_functions_to_view}
\alias{get_functions_to_view}
\title{get_functions_to_view}
\usage{
get_functions_to_view(
  keywords,
  mod = "inspect",
  overlapKeywords = FALSE,
  species = "mouse",
  cat = NULL
)
}
\arguments{
\item{keywords}{keywords to inspect or full names to extract}

\item{mod}{two values accepted. There are two modes. 'inspect' to see available terms to use. 'catch' to obtain gene sets with given terms.}

\item{overlapKeywords}{a boolean value. Default is FALSE. If TRUE, only terms with all keywords occurred in the term name will be used.}

\item{species}{the species used for the analysis.}

\item{cat}{MSigDB category to search. If want to search all, set cat to NULL. Default is C5.}
}
\value{
a list of given functions
}
\description{
get the gene sets of the query functions with given names.
}
