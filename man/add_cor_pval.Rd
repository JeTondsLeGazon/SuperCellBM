% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rank_gene_gene_cor.R
\name{add_cor_pval}
\alias{add_cor_pval}
\title{computes p value of gene gene correlations}
\usage{
add_cor_pval(gg.rank, ge, cell.size = NULL, mean1 = FALSE)
}
\arguments{
\item{gg.rank}{output of \link{rank_gene_gene_cor}}

\item{ge}{super-cell gene expression matrix}

\item{membership}{super-cell memberrship vector}
}
\value{
add p value and adjusted p value columns to \code{gg.rank}
}
\description{
computes p value of gene gene correlations
}
