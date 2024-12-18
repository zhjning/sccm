% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_causal_mediation.R
\name{generate_triple_ct_model}
\alias{generate_triple_ct_model}
\title{generate_triple_ct_model}
\usage{
generate_triple_ct_model(alist_of_ct, cells_to_fix, keep_or_remove = "keep")
}
\arguments{
\item{alist_of_ct}{a character list of the name of cells}

\item{cells_to_fix}{a named list or character to filter the the result. position 1,2,3 represent sender, mediator, receiver, respectively. Names of the vector corresponding to the name of the cell types used in alist_of_ct.}

\item{keep_or_remove}{two values can be identified. "keep" to keep the combinations satisfying cells_to_fix's rules. "remove" to remove the combinations satifying cells_to_fix's rules.}
}
\value{
a combination of the sender, mediator, and receiver to use for the causal mediation analysis.
}
\description{
generate triple cell type models

cells to fix: named list, name for cell position, value for cell name, use Null or "" as the name if don't need to filter the fixed position;
}
