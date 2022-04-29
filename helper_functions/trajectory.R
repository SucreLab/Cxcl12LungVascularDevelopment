# Title     : Expression trajectories
# Objective : Functions to rapidly categorize expression trajectories
# Created by: Nick Negretti
# Created on: 9/29/20


#' Get a wide dataframe of gene expression means by timepoint
#'
#' @param seurat_obj seurat object to run on, needs to have a $timepoint variable
#' @param only_vargenes A bool, only include variable genes (Default : `TRUE`)
#' @param expression_cutoff exclude genes with less than this much expression base on a sum of all timepoints (Default : `0.5`)
#' @return Returns Wide dataframe
#' @examples
#' Entire workflow:
#' expression_means <- getExpressionMeans(seurat_obj)
#'
#' group_labels <- labelExpressionGroups(expression_groups)
#'
#'
getExpressionMeans <- function(seurat_obj, only_vargenes = TRUE, expression_cutoff = 0.5, gene_list = NULL){

  if (only_vargenes){
    # Trying to be memory efficient with the pipes
    return(
      lapply(unique(seurat_obj$timepoint), function(x){
        rowMeans(as.matrix(GetAssayData(seurat_obj)[VariableFeatures(seurat_obj),
                                                    seurat_obj$timepoint == x]), na.rm = TRUE)
      }) %>%
        do.call(rbind, .) %>%
        .[,colSums(., na.rm = TRUE) > expression_cutoff]
    )
  } else {
    if (!is.null(gene_list)){
      return(
        lapply(unique(seurat_obj$timepoint), function(x){
          rowMeans(as.matrix(GetAssayData(seurat_obj)[gene_list,
                                                      seurat_obj$timepoint == x]), na.rm = TRUE)
        }) %>%
          do.call(rbind, .) %>%
          .[,colSums(., na.rm = TRUE) > expression_cutoff]
      )
    } else {
      return(
        lapply(unique(seurat_obj$timepoint), function(x){
          rowMeans(as.matrix(GetAssayData(seurat_obj)[,seurat_obj$timepoint == x]), na.rm = TRUE)
        }) %>%
          do.call(rbind, .) %>%
          .[,colSums(., na.rm = TRUE) > expression_cutoff]
      )
    }



  }
}

#' Determine what pattern each expression group represents
#'
#' @param expression_groups output from createExpressionGroups
#' @return Returns a vector with each gene and a label
labelExpressionGroups <- function(expression_means){
  if (is.vector(expression_means)) {
    warning("Cluster that only exists at one timepoint, unable to correlate")
    return(rep(NA, length(expression_means)))
  } else if(nrow(expression_means) == 2){
    warning("Cluster that only exists at two timepoints, unable to correlate")
    return(rep(NA, ncol(expression_means)))
  }

  ntimepoint <- nrow(expression_means)
  ngene <- ncol(expression_means)
  #### Expression pattern detector
  ## Model the five expected patterns.
  correlations <- data.frame(increasing = seq(0, 1, length.out = ntimepoint),
                             decreasing = seq(1, 0, length.out = ntimepoint),
                             mid_high = (sin(seq(-1, 1, length.out = ntimepoint) + pi/2) - 0.5) * 2,
                             mid_low = (sin(seq(-1, 1, length.out = ntimepoint) + pi/2 + pi) + 1) * 2,
                             no_change = rep(0.5, ntimepoint)
  )

  out <- future_sapply(c(1:4), function(x){ # Five expression patterns
    future_sapply(c(1:ngene), function(y){
      cor(scale(expression_means)[,y], scale(correlations)[,x])
    })
  })

  rownames(out) <- colnames(expression_means)

  apply(out, 1, function(x){
    if (max(x) > 0.7){
      return(colnames(correlations)[which(x == max(x))])
    } else {
      return(colnames(correlations)[5])
    }
  })
}



get_lower_ci <- function(data){
  b_ci_out <- suppressWarnings(boot.ci(boot(data = data,statistic = function(x,i) median(x[i]),R = 1000)))
  return(b_ci_out$percent[length(b_ci_out$percent) - 1])
}
get_upper_ci <- function(data){
  b_ci_out <- suppressWarnings(boot.ci(boot(data = data,statistic = function(x,i) median(x[i]),R = 1000)))
  return(b_ci_out$percent[length(b_ci_out$percent)])
}