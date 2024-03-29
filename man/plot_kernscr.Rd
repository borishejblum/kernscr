% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_kernscr.R
\name{plot_kernscr_methodsplit}
\alias{plot_kernscr_methodsplit}
\alias{plot_kernscr_kernelsplit}
\title{Plotting functions used in the manuscript}
\usage{
plot_kernscr_methodsplit(
  raw_melted,
  adj_melted,
  kernel,
  method,
  pathway_names = TRUE,
  title = NULL,
  raw_lower_threshold = round(log10(1/10000), 1),
  adj_lower_threshold = round(log10(1/10000 * 70), 1)
)

plot_kernscr_kernelsplit(
  raw_melted,
  adj_melted,
  kernel,
  method,
  pathway_names = TRUE,
  title = NULL,
  raw_lower_threshold = round(log10(1/10000), 1),
  adj_lower_threshold = NULL
)
}
\description{
Plotting functions used in the manuscript
}
\keyword{internal}
