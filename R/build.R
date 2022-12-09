#' Build eNODAL object
#' @param Z A matrix/dataframe of omics data
#' @param Meta1 A dataframe of first set of meta variables.
#' @param Meta2 A dataframe of second set of meta variables.
#' @param Phenotype A dataframe of phenotype output.
#' @param Meta1_name Name of experiment conditions in Meta1. If it is NULL will name "Meta1".
#' @param Meta2_name Name of experiment conditions in Meta2. If it is NULL will name "Meta2".
#' @param ... Other parameters can pass to eNODAL_build
#' See \code{\link{createParams}}
#' @return An eNODAL object
#' @examples
#' data(Proteomics)
#' Z = Proteomics$Z
#' Meta1 = Proteomics$Meta1
#' Meta2 = Proteomics$Meta2
#' eNODAL_obj = eNODAL_build(Z, Meta1, Meta2)
#' @export

eNODAL_build <- function(Z, Meta1, Meta2, Phenotype = NULL,
                         Meta1_name = NULL, Meta2_name = NULL,
                         ...
                         ){

  if(is.data.frame(Z) || is.matrix(Z)){
    Z <- as.matrix(Z)
  }
  if(is.null(Meta1_name)){
    Meta1_name0 = "Meta1"
    Meta1_name = "Meta1"
  }else{
    Meta1_name0 = Meta1_name
  }

  if(is.null(Meta2_name)){
    Meta2_name0 = "Meta2"
    Meta2_name = "Meta2"
  }else{
    Meta2_name0 = Meta2_name
  }

  if(is.data.frame(Meta1)){

    Meta1 <- as.data.frame(Meta1)
    check_meta1 <- check_fac_num(Meta1)
    meta1_type <- check_meta1$type
    Meta1 <- check_meta1$X

  }else if(is.factor(Meta1)){
    Meta1 <- data.frame(Meta1)
#    Meta1 <- check_fac_num(Meta1)$X
    colnames(Meta1) <- Meta1_name
    meta1_type <- "fac"
  }else if(is.vector(Meta1)){
    Meta1 <- data.frame(Meta1)
    colnames(Meta1) <- Meta1_name
    meta1_type <- "num"
  }

  Meta1_name <- colnames(Meta1)

  if(is.data.frame(Meta2)){

    Meta2 <- as.data.frame(Meta2)
    check_meta2 <- check_fac_num(Meta2)
    meta2_type <- check_meta2$type
    Meta2 <- check_meta2$X
  }else if(is.factor(Meta2)){
    Meta2 <- data.frame(Meta2)
#    Meta2 <- check_fac_num(Meta2)$X
    colnames(Meta2) <- Meta2_name
    meta2_type <- "fac"
  }else if(is.vector(Meta2)){
    Meta2 <- data.frame(Meta2)
    colnames(Meta2) <- Meta2_name
    meta2_type <- "num"
  }
  Meta2_name <- colnames(Meta2)

  if(length(Meta1_name) > 1){
    Meta1_tmp <- paste(Meta1_name, collapse = "+")
  }else{
    Meta1_tmp <- Meta1_name
  }

  if(length(Meta2_name) > 1){
    Meta2_tmp <- paste(Meta2_name, collapse = "+")
  }else{
    Meta2_tmp <- Meta2_name
  }

  if(!is.null(Phenotype)){
    if(is.matrix(Phenotype)){
      Phenotype <- as.data.frame(Phenotype)
      check_pheno <- check_fac_num(Phenotype)
      Phenotype_type <- check_pheno$type
    }else if(is.vector(Phenotype)||is.factor(Phenotype)){
      Phenotype <- data.frame(Phenotype)
    }
  }

  params <- createParams(...)

  lm_full <- paste("Y~(",  Meta1_tmp, ")*(", Meta2_tmp, ")", sep = "")
  lm_int_alt <- paste("Y~(",  Meta1_tmp, "):(", Meta2_tmp, ")", sep = "")
  lm_int_null <- paste("Y~(",  Meta1_tmp, ")+(", Meta2_tmp, ")", sep = "")
  #  lm_int_alt <- paste("Y~(",  Meta1_tmp, "):(", Meta2_tmp, ")", sep = "")
  lm_meta1 <- paste("Y~(",  Meta1_tmp, ")", sep = "")
  lm_meta2 <- paste("Y~(", Meta2_tmp, ")", sep = "")

  if((params$gam) == TRUE){
    gam_meta1_tmp <- paste(colnames(Meta1)[meta1_type == "num"], collapse = ",")
    gam_meta2_tmp <- paste(colnames(Meta2)[meta2_type == "num"], collapse = ",")
    gam_k <- params$gam_k
    if(is.null(gam_k)){
      if(sum(meta1_type == "num")){
        gam_tmp1 <- paste0("s(",gam_meta1_tmp)
      }else{
        gam_tmp1 <- ""
      }

      if(sum(meta2_type== "num")){
        gam_tmp2 <- paste0("s(",gam_meta2_tmp)
      }else{
        gam_tmp2 <- ""
      }
    }else{
      if(sum(meta1_type== "num")){
        gam_tmp1 <- paste0("s(",gam_meta1_tmp, ",k=", gam_k)
      }else{
        gam_tmp1 <- ""
      }

      if(sum(meta2_type== "num")){
        gam_tmp2 <- paste0("s(",gam_meta2_tmp, ",k=", gam_k)
      }else{
        gam_tmp2 <- ""
      }
    }
    gam1_fac <- paste0(colnames(Meta1)[meta1_type == "fac"], collapse = "+")
    gam2_fac <- paste0(colnames(Meta2)[meta2_type == "fac"], collapse = "+")

    gam_tmp_full <- ""
    for(i in 1:sum(meta2_type == "fac")){
      if(nchar(gam_tmp1) > 0){
        gam_tmp_full <- paste0(gam_tmp_full, gam_tmp1, ",by=",
                               colnames(Meta2)[meta2_type == "fac"][i], ")+")
      }
    }

    for(i in 1:sum(meta2_type == "fac")){
      if(nchar(gam_tmp2) > 0){
        gam_tmp_full <- paste0(gam_tmp_full, gam_tmp2, ",by=",
                               colnames(Meta1)[meta2_type == "fac"][i], ")+")
      }
    }

    if(sum(meta1_type == "fac")){
      gam_tmp_full <- paste0("Y~",gam_tmp_full,
                         paste0(colnames(Meta1)[meta1_type == "fac"],
                                collapse = "+"))
    }else{
      gam_tmp_full <- paste0("Y~",gam_tmp_full)
    }

    # Remove last `+`
    if(sum(meta2_type == "fac")){
      gam_full <- paste0(gam_tmp_full,
                         paste0(colnames(Meta2)[meta2_type == "fac"],
                                collapse = "+"))
    }else{
      gam_full <- substr(gam_tmp_full, 1, nchar(gam_tmp_full) - 1)
    }
    gam_int_null <- "Y~"
    gam_meta1 <- "Y~"
    gam_meta2 <- "Y~"
    if(nchar(gam_tmp1)){
      gam_int_null <- paste0(gam_int_null,gam_tmp1,")+")
      gam_meta1 <- paste0(gam_meta1, gam_tmp1,")+")
    }
    if(nchar(gam_tmp2)){
      gam_int_null <- paste0(gam_int_null,gam_tmp2,")+")
      gam_meta2 <- paste0(gam_meta2, gam_tmp2,")+")
    }
    if(nchar(gam1_fac)){
      gam_int_null <- paste0(gam_int_null,gam1_fac,"+")
      gam_meta1 <- paste0(gam_meta1, gam1_fac, "+")
    }
    if(nchar(gam2_fac)){
      gam_int_null <- paste0(gam_int_null,gam2_fac,"+")
      gam_meta2 <- paste0(gam_meta2, gam2_fac, "+")
    }
    gam_int_null <- substr(gam_int_null, 1, nchar(gam_int_null) - 1)
    gam_meta1 <- substr(gam_meta1, 1, nchar(gam_meta1) - 1)
    gam_meta2 <- substr(gam_meta2, 1, nchar(gam_meta2) - 1)

  }else{
    gam_full = gam_int_null = gam_meta1 = gam_meta2 = ""
  }


#  Meta_tmp <- as.data.frame(cbind(Meta1, Meta2))
#  Design_all <- model.matrix(formula(model_all), data = Meta_tmp)[,-1]
#  Design_int_null <- model.matrix(formula(model_int_null), data = Meta_tmp)[,-1]
#  Design_int_alt <- model.matrix(formula(model_int_alt), data = Meta_tmp)[,-1]
#  Design_meta1 <- model.matrix(formula(model_meta1), data = Meta_tmp)[,-1]
#  Design_meta2 <- model.matrix(formula(model_meta2), data = Meta_tmp)[,-1]

#  grp <- rep("Interaction_eff", ncol(Design_all))
#  grp[1:(ncol(Design_int_null))] <- c(rep("Meta1", ncol(Design_meta1)),
#                                          rep("Meta2", ncol(Design_meta2)))
  eNODAL_middle <- list()
  formula_list <- list(linear = c("full" = lm_full, "int_null" = lm_int_null,
                                  "int_alt" = lm_int_alt,
                                  "meta1" = lm_meta1, "meta2" = lm_meta2),
                       gam = c("full" = gam_full, "int_null" = gam_int_null,
                               "meta1" = gam_meta1, "meta2" = gam_meta2))
  eNODAL_middle$formula = formula_list
  eNODAL_middle$type = list(Meta1 = meta1_type, Meta2 = meta2_type)
  eNODAL_middle$name = list(Meta1 = Meta1_name0, Meta2 = Meta2_name0)
  eNODAL_obj <- eNODAL(Z = Z,
                       Meta1 = Meta1,
                       Meta2 = Meta2,
                       eNODAL_middle = eNODAL_middle,
                       params = params,
                       Phenotype = Phenotype)
  return(eNODAL_obj)
}

#' Create parameters list for eNODAL object.
#' @param gam_k Parameter k in GAM model. By default is 29.
#' @param gam Whether use create GAM formula. By default is TRUE.
#' If no continuous variable, set FALSE.
#' @param test_method Testing method.
#' Can be chosen from "F", "globaltest", "Tmax", "Chisq", "Cp". By default is "F"
#' @param h_adj Indicator of whether use hierarchical adjustment. By default is TRUE
#' @param sig_level Significance level of different test. By default is set below:
#' LC = 0.01, Linear = 0.01, Interaction = 0.05, Meta1 = 0.05, Meta2 = 0.05.
#' @param adapt Indicator of whether to use adaptive way to decide testing type.
#' @param test_func Testing function used. By default is "gam".
#' @param save_dist Save distance result. By default is FALSE.
#' @param q0 Parameter q0 in apcluster. Byt default is 0. If it is NULL, will not use apcluster.
#' @param knn.k Number of k-nearest neighbour to create knn graph in Louvian algorithm.
#' By default is 40, If it is NULL, will not use Louvian.
#' @param dist_thresh Threshold of distance when creating knn graph in Louvian algorithm.
#' By default is 0.4, If it is NULL, will not use Louvian.
#' @param eps Parameter eps in dbscan. If it is NULL, will not use dbscan.
#' @param minPts Parameter minPts in dbscan. If it is NULL, will not use dbscan.
#' @param minClusterSize Parameter minClusterSize in dbscan. If it is NULL, will not use dbscan.
#' @param kcluster Number of clusters if using fixed number of clustering method.
#' @param consensus_param Parameter list for consensus clustering.
#' Format: list(clmethod = ..., params = ...), available method and parameters listed above.
#' By default is Louvian, with knn.k = 15, dist_thresh = 0.3.
#' @param adj_method Pvalue adjust method. See p.adjust. By default is "BH".
#' @param sig_test Method for test between sig vs. non-sig. Can be choosen from "LC-test", "lm" or "gam.
#' @param ... Other parameters can be passed to eNODAL.
#' @return A list contain parameters in eNODAL object.
#' @export

createParams <- function(gam_k = 29, gam = TRUE, test_method = "F",
                         h_adj = TRUE, adj_method = "BH",
                         sig_level = c("LC" = 0.01, "Linear" = 0.01,
                                       "Sig" = 0.01, "Interaction" = 0.05,
                                       "Meta1" = 0.05, "Meta2" = 0.05),
                         adapt = FALSE, test_func = "gam", sig_test = "LC-test",
                         save_dist = FALSE, q0 = 0, knn.k = 40, dist_thresh = 0.4,
                         eps = NULL, minPts = NULL, minClusterSize = NULL,
                         kcluster = NULL,
                         consensus_param = list(clmethod = "louvian",knn.k = 15,
                                                dist_thresh = 0.75), ...){

  Clparam <- createClParam(q0, knn.k, dist_thresh, eps, minPts, minClusterSize)

  return(list(gam_k = gam_k, gam = gam, test_method = test_method, save_dist = save_dist,
              h_adj = h_adj, adj_method = adj_method, sig_level = sig_level,
              adapt = adapt, test_func = test_func, sig_test = sig_test,
              Clparam = Clparam, consensus_param = consensus_param))
}

#' Create parameters list for stage II clustering.
#' @param q0 Parameter q0 in apcluster. Byt default is 0. If it is NULL, will not use apcluster.
#' @param knn.k Number of k-nearest neighbour to create knn graph in Louvian algorithm.
#' By default is 40, If it is NULL, will not use Louvian.
#' @param dist_thresh Threshold of distance when creating knn graph in Louvian algorithm.
#' By default is 0.4, If it is NULL, will not use Louvian.
#' @param eps Parameter eps in dbscan. If it is NULL, will not use dbscan.
#' @param minPts Parameter minPts in dbscan. If it is NULL, will not use dbscan.
#' @param minClusterSize Parameter minClusterSize in dbscan. If it is NULL, will not use dbscan.
#' @return A list of methods and parameters.
createClParam <- function(q0 = 0, knn.k = 40, dist_thresh = 0.4,
                          eps = NULL, minPts = NULL, minClusterSize = NULL){
  Param0 <- list()
  if(!is.null(q0)){
    tmp_param <- list(clmethod = "apcluster", q = q0)
    Param0 <- append(Param0, list(tmp_param))
  }

  if(!is.null(knn.k) || !is.null(dist_thresh)){
    tmp_param <- list(clmethod = "louvian")
    if(!is.null(knn.k)){
      tmp_param$knn.k <- knn.k
    }
    if(!is.null(dist_thresh)){
      tmp_param$dist_thresh <- dist_thresh
    }
    Param0 <- append(Param0, list(tmp_param))
  }

  if(!is.null(eps) || !is.null(minPts)){
    tmp_param <- list(clmethod = "dbscan")
    if(!is.null(eps)){
      tmp_param$eps <- eps
    }
    if(!is.null(minPts)){
      tmp_param$minPts <- minPts
    }
    Param0 <- append(Param0, list(tmp_param))
  }

  if(!is.null(minClusterSize)){
    tmp_param <- list(clmethod = "Dcutree",
                      minClusterSize = minClusterSize)
    Param0 <- append(Param0, list(tmp_param))
  }
  return(Param0)
}

#' Print all the clustering of eNODAL result
#' @param eNODAL, an eNODAL object.
#' @return Clustering result of eNODAL object.
show_clusters <- function(eNODAL){

  Cl_df <- eNODAL@eNODAL_output$Cluster_res
  if(is.null(Cl_df)){
    cat("Currently no clustering result.")
  }else{
    sig0 <- Cl_df$Sig0
    if(!is.null(sig0)){
      sig0_tab <- table(sig0)
      text_tmp <- c()
      for(i in 1:length(sig0_tab)){
        text_tmp <- paste0(text_tmp, names(sig0_tab)[i], "(", sig0_tab[i], "), ")
      }
      text_tmp <- substr(text_tmp, 1, nchar(text_tmp) - 2)
      cat("Level0, significant vs. non-significant under experimental condition: \n",
                         text_tmp, "\n")
    }
    sig1 <- Cl_df$Sig1
    if(!is.null(sig1)){
      sig1_tab <- table(sig1)
      text_tmp <- c()
      for(i in 1:length(sig1_tab)){
        text_tmp <- paste0(text_tmp, names(sig1_tab)[i], "(", sig1_tab[i], "), ")
      }
      text_tmp <- substr(text_tmp, 1, nchar(text_tmp) - 2)
      cat("Level1, marginal effect vs. interaction effect: \n",
                      text_tmp,
                      "\n")
    }
    sig2 <- Cl_df$Sig2
    if(!is.null(sig2)){
      sig2_tab <- table(sig2)
      text_tmp <- c()
      for(i in 1:length(sig2_tab)){
        text_tmp <- paste0(text_tmp, names(sig2_tab)[i], "(", sig2_tab[i], "), ")
      }
      text_tmp <- substr(text_tmp, 1, nchar(text_tmp) - 2)
      cat("Level2, unsupervised clustering within each cluster: \n",
                     text_tmp)
    }
  }
}

#' Get clustering result of eNODAL object
#' @param eNODAL, an eNODAL object.
#' @return Clustering result of eNODAL object.
get_clusters <- function(eNODAL){
  Cl_df <- eNODAL@eNODAL_output$Cluster_res
  return(Cl_df)
}
