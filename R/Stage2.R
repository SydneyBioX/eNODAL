#' Run subclustering using with eNODAL object.
#' @param eNODAL Input eNODAL object.
#' @return eNODAL object with clustering result(Stage II).
#' @export
runSubclust <- function(eNODAL){

  ClParam <- eNODAL@params$Clparam
  Cl_df <- eNODAL@eNODAL_output$Cluster_res
  sig_levels <- na.omit(unique(Cl_df$Sig1))
  X <- eNODAL@Z
  V_exp <- eNODAL@eNODAL_middle$VarExplained
  if(is.null(V_exp)){
    eNODAL <- calcProp(eNODAL)
    V_exp <- eNODAL@eNODAL_middle$VarExplained
  }

  consensus_param <- eNODAL@params$consensus_param
  save_dist <- eNODAL@params$save_dist
#  if(!is.null(V_exp)){
#    V_prop = V_exp[,1:3]/rowSums(V_exp[,1:3])
#  }else{
#    V_prop = NULL
#  }

  Cl_df$Sig2 <- NA
  for(sig1 in sig_levels){
    idx_tmp <- which(Cl_df$Sig1 == sig1)
    X_sel <- X[,idx_tmp]
    V_sel <- V_exp[Cl_df$varname[idx_tmp],]
    Dist_list <- createDistList(X_sel, V_sel, dist_method = "all")
    if(save_dist){
      eNODAL@eNODAL_middle$Dist_list = Dist_list
    }
    cluster <- createConsensus(Dist_list, ClParam)
    Param_cons <- consensus_param
    Param_cons$Dist = 1 - cluster$cons
    cluster <- do.call(clusterMethod, Param_cons)$cluster
    cluster_name <- paste(sig1, "_C", cluster, sep = "")
    names(cluster_name) <- names(cluster)
    Cl_df$Sig2[match(names(cluster_name),Cl_df$varname)] <- cluster_name
  }
  eNODAL@eNODAL_output$Cluster_res <- Cl_df
  return(eNODAL)
}


#' Run clustering method (do not need to specify number of clusters)
#' @param Dist distance matrix to calculate
#' @param clmethod Method to run for clustering, can be chosen from:
#' apcluster: Affinity propagation clustering,
#' louvian: Louvian algorithm use knn graph,
#' dbscan: Density-based spatial clustering of applications with noise,
#' Dcutree: Adaptive branch pruning of hierarchical clustering.
#' @param ... Other parameters can be passed to clusterMethod.
#' @importFrom stats hclust as.dist
#' @return A list contain clustering result.

clusterMethod <- function(Dist, clmethod, ...){

  if(clmethod == "apcluster"){
    apres <- apcluster::apcluster(s = 1 - Dist, ...)
    exemplars <- apres@exemplars
    cluster_tmp <- apres@clusters
    clusters <- c()
    for(i in 1:length(cluster_tmp)){
      clusteri <- c(cluster_tmp[[i]])
      tmp <- rep(i, length(clusteri))
      names(tmp) <- names(clusteri)
      clusters <- c(clusters, tmp)
    }
    clusters <- clusters[colnames(Dist)]
    return(list(cluster = clusters, exemplars = names(exemplars)))
  }else if(clmethod == "louvian"){
    knn_graph <- make.knn.graph(Dist, ...)
    g <- knn_graph
    louvianres <- igraph::cluster_louvain(g)
    clusters <- louvianres$membership
    names(clusters) <- colnames(Dist)
    return(list(cluster = clusters, graph = g, layout = knn_graph$layout))
  }else if(clmethod == "dbscan"){
    dbcluster <- dbscan::dbscan(as.dist(Dist),...)
    clusters <- dbcluster$cluster + 1
    names(clusters) <- colnames(Dist)
    return(list(cluster = clusters))
  }else if(clmethod == "Dcutree"){
    hcluster <- hclust(as.dist(Dist))
    clusters <- dynamicTreeCut::cutreeDynamic(hcluster, distM = Dist, verbose = 0, ...)
    names(clusters) <- colnames(Dist)
    return(list(cluster = clusters))
  }

}

#' Run clustering method (Need to specify number of clusters)
#' @param Dist distance matrix to calculate
#' @param k Number of clusters
#' @param clmethod Method to run for clustering, can be chosen from:
#' hclust: Hierachical clustering,
#' apcluster: Affinity propagation clustering with k clusters,
#' kmeans: kmeans clustering,
#' @param ... Other parameters can be passed to clusterMethodk.
#' @importFrom stats hclust cutree as.dist
#' @return A list contain clustering result.

clusterMethodk <- function(Dist, clmethod, k = 5, ...){

  if(clmethod == "hclust"){
    hcluster <- hclust(as.dist(Dist),...)
    clusters = cutree(hcluster, k)
    return(list(cluster = clusters))
  }else if(clmethod == "apcluster"){
    apres <- apcluster::apclusterK(s = 1 - Dist, K = k,...)
    exemplars <- apres@exemplars
    cluster_tmp <- apres@clusters
    clusters <- c()
    for(i in 1:length(cluster_tmp)){
      clusteri <- c(cluster_tmp[[i]])
      tmp <- rep(i, length(clusteri))
      names(tmp) <- names(clusteri)
      clusters <- c(clusters, tmp)
    }
    clusters <- clusters[colnames(Dist)]
    return(list(cluster = clusters, exemplars = names(exemplars)))}
#   Partition around median
#  }else if(clmethod == "pam"){
#    pamcluster <- pam(as.dist(Dist), k = k, diss = T,...)
#    clusters <- pamcluster$clustering
#    names(clusters) <- colnames(Dist)
#    exemplars <- pamcluster$medoids
#    return(list(cluster = clusters, exemplars = names(exemplars)))
#  }
}

#clusterMethodKx <- function(x, method, k = 5,...){
#  if(method == "kmeans"){
#    kmcluster <- kmeans(x, centers = k,...)
#    clusters <- kmcluster$cluster
#    names(clusters) <- colnames(Dist)
#    return(list(cluster = clusters))

#  }else if(method == "pam"){
#    pamcluster <- pam(x, k = k,...)
#    clusters <- pamcluster$clustering
#    names(clusters) <- colnames(Dist)
#    exemplars <- pamcluster$medoids
#    return(list(cluster = clusters, exemplars = names(exemplars)))
#  }
#}

#' Create distance matrix with respect to the variance explained.
#' @param V_exp Variance explained matrix (see calcProp).
#' @return An eNODAL object with calculated Atchison distance.
createwD2 <- function(V_exp){

#  Cl_df <- eNODAL@eNODAL_output$Cluster_res
#  V_exp <- eNODAL@eNODAL_middle$VarExplained
#  if(is.null(V_exp)){
#    eNODAL <- calcProp(eNODAL)
#    V_exp <- eNODAL@eNODAL_middle$VarExplained
#  }
  V_prop <- V_exp[,1:3]
  V_prop <- V_prop/rowSums(V_prop)

#  V_prop0 <- V_prop[Cl_df$varname[which(Cl_df$Sig1 == Cl_sel)],,drop = F]
  V_props <- colSums(V_prop)
  V_prop1 <- V_prop[,V_props != 0,drop = F]
  if(ncol(V_prop1) == 1){
    D2 <- matrix(0, nrow = nrow(V_prop1), ncol = nrow(V_prop1))
    rownames(D2) <- rownames(V_prop1)
    colnames(D2) <- rownames(V_prop1)
  }else{
    V_prop1 <- V_prop1/rowSums(V_prop1)
    D2 <- calcdist(V_prop1, "Atchison")
  }

  return(D2)
}

#createwD <- function(eNODAL, manual_Feature, manual_grp, V_prop, D2 = NULL,
#                     sig_level = c("Meta1" = 0.01, "Meta2" = 0.05, "Int" = 0.1),
#                     trans = "none", weighted_type = 1){

#  manual_Feature[abs(manual_Feature) < qnorm(1 - sig/2)] <- 0
#  if(trans == "binary"){
#    manual_Feature <- sign(manual_Feature)
#  }

#  if(is.null(D2)){
#    D2 <- createwD2(V_prop)
#    if(max(D2)!=0){
#      D2_norm <- D2/max(D2)
#    }else{
#      D2_norm <- D2
#    }

#  }

#  feat_num <- table(manual_grp)
#  if(weighted_type == 1){
#    wmanual_Feature <- manual_Feature
#    wmanual_Feature[,manual_grp == "N"] <- manual_Feature[,manual_grp == "N"] * V_prop[,"N"] / sqrt(feat_num["N"])
#    wmanual_Feature[,manual_grp == "T"] <- manual_Feature[,manual_grp == "T"] * V_prop[,"T"] / sqrt(feat_num["T"])
#    wmanual_Feature[,manual_grp == "NxT"] <- manual_Feature[,manual_grp == "NxT"] * V_prop[,"NxT"] / sqrt(feat_num["NxT"])
#    D1 <- as.matrix(dist(wmanual_Feature))
#  }else if(weighted_type == 2){
#    wMat_N <- outer(V_prop[,"N"], V_prop[,"N"], `+`)/2
#    wMat_T <- outer(V_prop[,"T"], V_prop[,"T"], `+`)/2
#    wMat_I <- outer(V_prop[,"NxT"], V_prop[,"NxT"], `+`)/2
#
#    D_N <- as.matrix(dist(manual_Feature[,manual_grp == "N"]))
#    D_T <- as.matrix(dist(manual_Feature[,manual_grp == "T"]))
#    D_I <- as.matrix(dist(manual_Feature[,manual_grp == "NxT"]))
#    D1 <- wMat_N*D_N/feat_num["N"] + wMat_T*D_T/feat_num["T"] +
#      wMat_I*D_I/(feat_num["NxT"])

#  }else if(weighted_type == 3){
#    D1 <- as.matrix(dist(manual_Feature))
#  }

#  if(max(D2) != 0){
#    D2_norm <- D2/max(D2)
#  }else{
#    D2_norm <- D2
#  }

#  D1_norm <- D1/max(D1)
#  wD <- (1 - D2_norm)*D1_norm + D2_norm*D2_norm
#  return(list(D1 = D1, D2 = D2, wD = wD,
#              Feature = manual_Feature,
#              Feature_grp = manual_grp))
#}

#' Run correlation distance
#' @param X_sel Matrix used to calculate correlation distance()
#' @param dist_type Type of correlation to calculate distance, can be chosen from
#' "pearson" and "spearman"
#' @return Distance matrix

createDist_base <- function(X_sel, dist_type){

  if(dist_type == "pcc"){
    Dist <- calcdist(X_sel, "pearson")
  }else if(dist_type == "spc"){
    Dist <- calcdist(X_sel, "spearman")
  }
  return(Dist)
}

#' Calculate weighted distance
#' @param X_sel Matrix used to calculate correlation distance.
#' @param V_exp Matrix of explained variance (see calcProp).
#' @param Dist_base Base distance, create from createDist_base.
#' @param D2 Distance matrix for weight.
#' @return Weighted distance matrix
createwDist <- function(X_sel, V_exp, Dist_base = NULL, D2 = NULL){

  if(is.null(D2)){
    D2 <- createwD2(V_exp)
  }
  if(max(D2) > 0){
    D2_norm <- D2/max(D2)
  }else{
    D2_norm <- D2
  }

  Dist_norm <- Dist_base/max(Dist_base)
  wDist <- (1 - D2_norm) *Dist_norm + D2_norm^2

  return(wDist)
}

#' Calculate list of distance matrix
#' @param X_sel Matrix used to calculate correlation distance.
#' @param V_exp Matrix of explained variance (see calcProp).
#' @param dist_method Types of distance matrix to be calculated.
#' Can be chosen from "base" or "all".
#' @return Weighted distance matrix
#'
createDistList <- function(X_sel, V_exp = NULL,
                           dist_method = "all"){

  if(dist_method == "base"){
    D_pcc <- createDist_base(X_sel, dist_type = "pcc")
    D_pcc <- D_pcc/max(D_pcc)
    D_spc <- createDist_base(X_sel, dist_type = "spc")
    D_spc <- D_spc/max(D_spc)
    Dist_list <- list(D_pcc = D_pcc, D_spc = D_spc)
  }else if(dist_method == "all"){
    D_pcc <- createDist_base(X_sel, dist_type = "pcc")
    D_pcc <- D_pcc/max(D_pcc)
    D_spc <- createDist_base(X_sel, dist_type = "spc")
    D_spc <- D_spc/max(D_spc)
    D2 <- createwD2(V_exp)
    wD_pcc <- createwDist(X_sel, V_exp, D_pcc, D2)
    wD_spc <- createwDist(X_sel, V_exp, D_spc, D2)
#    wD1 <- createwDist(X_sel, "wD", V_prop, manual_Feature, manual_grp, Dist_base = NULL, D2, weighted_type = 1,...)
#    wD2 <- createwDist(X_sel, "wD", V_prop, manual_Feature, manual_grp, Dist_base = NULL, D2, weighted_type = 2,...)
#    wD3 <- createwDist(X_sel, "wD", V_prop, manual_Feature, manual_grp, Dist_base = NULL, D2, weighted_type = 3,...)
    Dist_list <- list(D_pcc = D_pcc, D_spc = D_spc, wD_pcc = wD_pcc, wD_spc = wD_spc)
#    ,wD1 = wD1, wD2 = wD2, wD3 = wD3)
  }
  return(Dist_list)
}

#' Consensus clustering
#' @param Dist_list A list of distance matrix(see createDistList)
#' @param Param0 Parameters for consensus clustering(see setParam)
#' @return Consensus matrix.
createConsensus <- function(Dist_list, Param0){

  #  Param0 <- createClParam(...)
  num_var <- ncol(Dist_list[[1]])
  cons <- matrix(0, nrow = num_var, ncol = num_var)
  k <- 0
  for(i in 1:length(Dist_list)){
    for(j in 1:length(Param0)){
      Dist_tmp <- Dist_list[[i]]
      param <- Param0[[j]]
      param$Dist <- Dist_tmp/max(Dist_tmp)
      cluster <- do.call(clusterMethod, param)
      if(length(unique(cluster$cluster)) > 1){
        Adj_mat <- cvtCooccur(cluster$cluster)
        cons <- cons + Adj_mat
        k <- k+1
      }
    }
  }
  colnames(cons) <- colnames(Dist_list[[1]])
  rownames(cons) <- rownames(Dist_list[[1]])
  cons0 <- cons
  cons <- cons/k

  return(list(cons = cons, cons0 = cons0))
}
