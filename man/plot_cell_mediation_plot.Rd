% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_relationships.R
\name{plot_cell_mediation_plot}
\alias{plot_cell_mediation_plot}
\title{plot_cell_mediation_plot}
\usage{
plot_cell_mediation_plot(
  fit,
  maxx = 80,
  maxy = 30,
  height = 5,
  width = 5,
  base_size = 5,
  base_family = "Arial",
  title_text = "",
  title_size = 10,
  usecolor = TRUE,
  cols = NULL,
  clean = TRUE,
  whatLabels = "std",
  useLabel = FALSE,
  regression = TRUE,
  mode = 1
)
}
\arguments{
\item{fit}{fitted mediation model generated by sem from lavaan}

\item{maxx}{the max x value, default 80}

\item{maxy}{the max y value, default 30}

\item{height}{the height of the figure, default 5}

\item{width}{the width of the figure, default 5}

\item{base_size}{the font size used for text}

\item{base_family}{the font used for text}

\item{usecolor}{default TRUE}

\item{cols}{the color list used for square objects, default NULL}

\item{clean}{the ggplot theme}

\item{whatLabels}{default value 'std'}

\item{regression}{a boolean value}

\item{mode}{default 1}
}
\value{
a ggplot object
}
\description{
plot cell mediation plot

The funciton generates the cell mediation relationships.
Each square represents a kind of cell type.
The solid line indicates a significant regulation. The arrow towards the object being regulated.
The dashed line indicates a non-significant regulation. The value on the line indicate the estimated coefficient.
}
