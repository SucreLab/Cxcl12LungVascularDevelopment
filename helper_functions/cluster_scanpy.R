# Title     : TODO
# Objective : TODO
# Created by: negretn
# Created on: 9/29/20

# defaults are set to match those from Seurat
#' Run clustering in Scanpy with leiden algorithm using retibulate.
#' Defaults are set to match those from Seurat. By default, this writes the seurat object to disk
#' for faster computation later.
#'
#' @param obj Seurat object
#' @param filename where to save the transitional h5ad
#' @param del_temp_files bool, delete the temporary files after computation?
#' @param resolution Resolution for clustering
#' @param n_neighbors How many neighbors for clustering
#' @param n_pcs number of dimensions for calculating clusters
#' @return An ordered factor, put in to Ident(obj)
cluster_with_scanpy <- function(obj, filename, del_temp_files = FALSE, resolution = 0.3, n_neighbors = 20, n_pcs = 15){
  # Save the seurat object to disk so scanpy can read it
  h5seurat <- paste0(filename, ".h5Seurat")
  h5ad <- paste0(filename, ".h5ad")

  if (!file.exists(h5seurat)) { SaveH5Seurat(obj, filename = h5seurat, overwrite = T, verbose = F) }
  if (!file.exists(h5ad)) { Convert(h5seurat, dest = "h5ad", overwrite = T, verbose = F) }

  # Run some python code in R with reticulate package (req for Seurat)
  sc <- import("scanpy")

  # load the data in to scanpy, find neighbors, and cluster
  sc_data <- sc$read(h5ad)
  sc$pp$neighbors(sc_data, n_neighbors=as.integer(n_neighbors), n_pcs=as.integer(n_pcs))
  sc$tl$leiden(sc_data, resolution=resolution)

  # clean up the extra files created during processing
  if (del_temp_files){
    if (file.exists(h5seurat)) { file.remove(h5seurat) }
    if (file.exists(h5ad)) { file.remove(h5ad) }
  }

  # extract the clusters from the adata object
  scanpy_clusters <- sc_data$obs['leiden']
  # return an ordered factor, to put in to Idents(seurat_obj)
  return(ordered(as.factor(scanpy_clusters[,'leiden'])))
}


#' Run the common Seurat clustring, PCA, and UMAP
#'
#' @param dims_umap dimensions to run UMAP
#' @param dims_neighbors How many dims for neighbor finding
#' @param k_param k for finding neighbors
#' @param cluster_res Resolution for clustering
#' @return Seurat object
cluster_pca_umap <- function(obj, dims_umap = 1:15, dims_neighbors = 1:15, k_param = 10, cluster_res = 0.3){
  if ("integrated" %in% names(obj)){
    DefaultAssay(obj) <- "integrated"
  }

  obj <- RunPCA(obj, verbose = F)
  obj <- RunUMAP(obj, dims = dims_umap, verbose = F)
  obj <- FindNeighbors(obj, dims = dims_neighbors, k.param = k_param)
  obj <- FindClusters(obj, resolution = cluster_res)
  if ("integrated" %in% names(obj)){
    DefaultAssay(obj) <- "SCT"
  }
  return(obj)
}

#' This can use a ton of memory, but it will save a lot of time
#' Important: The output is the same order as levels(Idents(obj))
#'
#' This needs to export the Seurat object to all workers - this takes a lot of RAM.
#'
#' @param obj Seurat object
#' @param n_cor number of CPU cores
#' @return List of data tables, one for each numeric cluster
parallelFindAllMarkers <- function(obj){

  all_markers <- future_lapply(levels(Idents(obj)), function(x){ # Five expression patterns
    FindMarkers(obj, ident.1 = x, ident.2 = NULL, test.use = "MAST")
  })

  return(value(all_markers))
}