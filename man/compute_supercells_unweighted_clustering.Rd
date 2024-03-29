% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_supercells_clustering.R
\name{compute_supercells_unweighted_clustering}
\alias{compute_supercells_unweighted_clustering}
\title{Coumputes unweighted clustering for super-cell-like list}
\usage{
compute_supercells_unweighted_clustering(
  SC.list,
  SC.GE.list = NULL,
  algs = c("hclust", "kmeans", "seurat"),
  pca_name = "SC_PCA",
  N.clusters = NULL,
  N.clusters.seq = NULL,
  N.comp = NULL,
  DO_silhouette = TRUE,
  ignore.gammas = c(),
  seurat_fields = c(),
  def.resolution = 0.8,
  step.resolution = 0.1,
  max.counter = 20,
  seed = 12345,
  verbose = FALSE
)
}
\arguments{
\item{SC.list}{super-cell-like list}
}
\description{
Coumputes unweighted clustering for super-cell-like list
}
