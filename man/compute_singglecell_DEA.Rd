% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_supercells_DEA.R
\name{compute_singglecell_DEA}
\alias{compute_singglecell_DEA}
\title{Computes DEA for single-cell data and for GT cell type annotation}
\usage{
compute_singglecell_DEA(
  sc.ge,
  clusters,
  ident.1 = NULL,
  ident.2 = NULL,
  seed = 12345,
  pval.thresh = 0.05,
  ...
)
}
\arguments{
\item{sc.ge}{single-cell gene expression}

\item{clusters}{single-cell clustering result or GT cell type annotation}

\item{...}{rest of the parameters of \link{compute_supercells_DEA}}
}
\description{
Computes DEA for single-cell data and for GT cell type annotation
}
