#' Add annotation for variable name
#' @param eNODAL_obj An fitted eNODAL object.
#' @param from Initial id type in eNODAL object.
#' @param to To id type in eNODAL_obj object.
#' @param db Annotation database.
#' @returns An eNODAL_obj object.
#' @export
runAnno <- function(eNODAL_obj, from = "UNIPROT",
                    to = c("SYMBOL", "ENTREZID"),
                    db = "org.Mm.eg.db"){

  Cl_df <- eNODAL_obj@eNODAL_output$Cluster_res
  var_name <- Cl_df$varname
  transfer_tmp <- clusterProfiler::bitr(var_name, fromType = from, toType = to,
                                        OrgDb = db, drop = FALSE)
  for(i in 1:length(to)){
    Cl_df[[to[i]]] = transfer_tmp[,i + 1]
  }
  eNODAL_obj@eNODAL_output$Cluster_res <- as.data.frame(Cl_df)

  return(eNODAL_obj)
}

#' Check each column of a dataframe is a factor or numerical
#' @param X A dataframe
#' @param max_level maximum number of levels,
#' i.e. if the unique value of a column less than max_level,
#' will set to be factor. By default is 7.
#' @return a list contain checked matrix X and type of values for each output.
#' @export
check_fac_num <- function(X, max_level = 7){

  if(is.null(ncol(X))){
    X <- matrix(X, ncol = 1)
  }
  p <- ncol(X)
  res <- c()
  for(i in 1:p){

    tmp <- X[,i]
    tmp <- as.factor(as.character(tmp))
    if(length(levels(tmp)) > max_level){
      res[i] <- "num"
      X[,i] <- as.numeric(as.character(tmp))
    }else{
      res[i] <- "fac"
      X[,i] <- tmp
    }
  }
  return(list(X= X, type = res))
}

#' Convert a dataframe with factor and numerical to numerical with dummy code.
#' @param X A dataframe
#' @param max_level maximum number of levels,
#' i.e. if the unique value of a column less than max_level,
#' will set to be factor. By default is 7.
#' @importFrom stats model.matrix
#' @return a list contain checked matrix X and type of values for each output.
#' @export
cvt_fac_num <- function(X, max_level = 7){

  check_res <- check_fac_num(X, max_level)
  type = check_res$type
  X <- check_res$X
  Y <- matrix(0, nrow = nrow(X), ncol = 1)
  type_Y <- c()
  name_Y <- c()
  if(is.null(colnames(X))){
    colnames(X) <- paste("V", 1:ncol(X), sep = "")
  }
  k = 1
  for(i in 1:ncol(X)){
    tmp <- X[,i]
    if(type[i] == "num"){

      Y <- cbind(Y, tmp)
      type_Y[k] <- "num"
      name_Y[k] <- colnames(X)[i]
      k <- k + 1
    }else{
      tmp <- psych::dummy.code(tmp)
      tmp <- tmp[,-1]
      if(is.vector(tmp)){
        Y <- cbind(Y, tmp)
        type_Y[k] <- "fac"
        name_Y[k] <- colnames(X)[i]
        k <- k + 1
      }else{
        Y <- cbind(Y, tmp)
        type_Y[k:(k + ncol(tmp) - 1)] <- "fac"
        name_Y[k:(k + ncol(tmp) - 1)] <- paste(colnames(X)[i], colnames(tmp), sep = "_")
        k <- k + ncol(tmp)
      }
    }
  }

  Y <- Y[,-1]
  colnames(Y) <- name_Y
  return(list(X = Y, type = type_Y))

}

#' Unify all the type of X into matrix
#' @param X Input data.
#' @return Matrix form of data.
.unify_class <- function(X){

  cls_X <- class(X)

  if(cls_X == "data.frame"){
    X <- as.matrix(X)
  }else if(!("matrix" %in% cls_X)){
    X <- matrix(X, ncol = 1)
  }
  return(X)
}

#' Make knn graph
#' @param Dist distance matrix
#' @param knn.k k-nearest neighbor. By default is 15
#' @param dist_thresh minimal distance threshold (k-nearest neighbor should smaller than dist_thresh)
#' @return The corresponding knn graph
#' @export
make.knn.graph<-function(Dist, knn.k = 15, dist_thresh = 0.3){
  knn.k <- min(knn.k, nrow(Dist) - 1)
  dist_thresh <- max(dist_thresh, min(Dist[upper.tri(Dist)])*1.2)
  edges <- mat.or.vec(0,2)
  for (i in 1:nrow(Dist)){

    matches <- setdiff(order(Dist[i,],decreasing = F)[1:(knn.k+1)],i)
    matches_dist <- Dist[i, matches]
    idx_int <- which(matches_dist < dist_thresh)
    if(length(idx_int) < 2){
      matches <- matches[1]
    }else{
      matches <- matches[matches_dist < dist_thresh]
    }

    edges <- rbind(edges,cbind(rep(i,length(matches)),matches))
    edges <- rbind(edges,cbind(matches,rep(i,length(matches))))
  }

  graph <- igraph::graph_from_edgelist(edges,directed=F)
#  igraph::V(graph)$frame.color <- NA

  return(graph)
}

#' Fisher Transformation
#' @param r correlation coefficient
#' @returns Fisher transformed coefficient

fisher_trans <- function(r){

  r[r > 1 - 1e-5] = 1-1e-5
  r[r < -1 + 1e-5] = -1 + 1e-5
  fr <- 0.5*log((1+r)/(1-r))
  return(fr)
}


#' Calculate Gower distance (among cols)
#' @param B A matrix for Gower distance calculation
#' @importFrom stats var
#' @returns A distance matrix.
#' @export
#'
calcGower <- function(B){
  G_B <- matrix(0, nrow = nrow(B), ncol = nrow(B))
  for(i in 1:ncol(B)){
    if(var(B[,i]) !=0 ){
      B_i <- B[,i]
      G_B_i <- outer(B_i, B_i,"-")
      G_B_i <- abs(G_B_i) / max(abs(G_B_i))
      G_B <- G_B + G_B_i
    }
  }
  return(G_B)
}


#' Calculate centered log ratio(clr) transformation of a matrix
#' @param B A matrix for clr transformation
#' @returns Transformed matrix.
#' @export
#'
trans_clr <- function(B){

  B[B < 1e-8] = 1e-8
  B <- t(apply(B, 1, function(x){x / sum(x)}))
  B_clr <- t(scale(log(t(B)), center = T, scale = F))

  return(B_clr)
}


#' Calculate Atchison distance of a matrix (among cols)
#' @param B A matrix
#' @importFrom stats dist
#' @returns Atchison distance matrix
#' @export
#'
calcAtchison <- function(B){

  B_clr <- trans_clr(B)
  aDist <- as.matrix(dist(B_clr))
  return(aDist)

}

#' Calculate different types of distance (among cols)
#' @param B A matrix
#' @param method Method for distance matrix calculation, can be chosen from
#' "gower", "euclid", "pearson", "spearman" and "Atchison"
#' @importFrom stats dist cor
#' @returns Transformed matrix.
#' @export
#'
calcdist <- function(B, method = "pearson"){

  if(method == "gower"){
    Dist <- calcGower(B)
    colnames(Dist) <- rownames(B)
    rownames(Dist) <- rownames(B)
  }else if(method == "euclid"){
    Dist <- as.matrix(dist(scale(B)))
    colnames(Dist) <- rownames(B)
    rownames(Dist) <- rownames(B)
  }else if(method == "pearson"){
    Dist <- as.matrix(1 - cor(B, method = "pearson"))
    colnames(Dist) <- colnames(B)
    rownames(Dist) <- colnames(B)
  }else if(method == "spearman"){
    Dist <- as.matrix(1 - cor(B, method = "spearman"))
    colnames(Dist) <- colnames(B)
    rownames(Dist) <- colnames(B)
  }else if(method == "Atchison"){
    Dist <- calcAtchison(B)
    colnames(Dist) <- rownames(B)
    rownames(Dist) <- rownames(B)
  }
  return(Dist)
}


#' Convert clustering result to co-occurrence matrix
#' @param cl_res A vector of clustering membership
#' @importFrom stats model.matrix
#' @return A matrix of co-occurrence.
#' @export
cvtCooccur <- function(cl_res){

  cl_res <- as.factor(cl_res)
  cl_onehot <- model.matrix(~ 0 + cl_res)

  res <- cl_onehot %*% t(cl_onehot)
  return(res)
}

#' Transformation of distance matrix
#' @param dists A matrix of distance
#' @param method transformation method, can be chosen from"pca" and "laplacian"
#' "pca" calculate the first several principal components and "laplacian" calculate last several eigen values.
#' @param d_max Maximum number of eigen decomposition. By default is 4.
#' @importFrom stats prcomp
#' @return A matrix of co-occurrence.
#' @export
calctrans <- function(dists, method, d_max = 4){
  d_max <- min(d_max, ncol(dists))
  #  dist_row <- rowMeans(dists)
  #  dist_col <- colMeans(dists)
  #  dist_mean <- mean(dists)
  #  dist_mds <- (-0.5*(dists - matrix(rep(dist_row, nrow(dists)), nrow = nrow(dists), byrow = F) -
  #                       matrix(rep(dist_col, ncol(dists)), nrow = nrow(dists), byrow = T) + dist_mean))
  if (method == "pca") {
    t <- prcomp(dists, scale. = TRUE, center = TRUE)
    return(t$x[,1:d_max])}
#  } else if (method == "laplacian") {
#    L <- norm_laplacian(dists)
#    l <- eigen(L)
    #    l <- prcomp(L, scale. = TRUE, center = TRUE)
    # sort eigenvectors by their eigenvalues
    #    return((l$vectors[, order(l$values)])[,1:d_max])
#    return((l$vectors[, order(l$values)])[,1:d_max])
#  }
}


#' Calculate proportion of each component in linear model
#' @param eNODAL_obj An eNODAL_obj object
#' @param index A number specify which column to calculate
#' @param Cl_sel A text specify which cluster this variable belongs to
#' @importFrom stats var lm formula anova
#' @return A vector contain proportion of variance each component explained.
calcProp_lm <- function(eNODAL_obj, index, Cl_sel){

  X <- eNODAL_obj@Z
  Meta1 <- eNODAL_obj@Meta1
  Meta2 <- eNODAL_obj@Meta2
  Xi <- X[,index]
  V_exp <- rep(0, 5)
  lm_formula <- eNODAL_obj@eNODAL_middle$formula$linear

  var_tot <- var(Xi) * (nrow(X) - 1)
  data_g <- data.frame(Y = Xi, Meta1, Meta2)

  if(Cl_sel == "Add"){

    fita <- lm(formula(lm_formula[2]), data = data_g)
    fitm1 <- lm(formula(lm_formula[3]), data = data_g)
    fitm2 <- lm(formula(lm_formula[4]), data = data_g)

    var_m1 <- anova(fitm2, fita)$`Sum of Sq`[2]
    var_m2 <- anova(fitm1, fita)$`Sum of Sq`[2]
    var_E <- anova(fita)$`Sum Sq`[2]

    V_exp[c(1,2,4,5)] <- c(var_m1, var_m2, var_E, var_tot)

  }else if(Cl_sel == "Int"){

    fitall <- lm(formula(lm_formula[1]), data = data_g)
    fita <- lm(formula(lm_formula[2]), data = data_g)
    fitm1 <- lm(formula(lm_formula[3]), data = data_g)
    fitm2 <- lm(formula(lm_formula[4]), data = data_g)

    var_Int <- anova(fita, fitall)$`Sum of Sq`[2]
    var_m1 <- anova(fitm2, fitall)$`Sum of Sq`[2]
    var_m2 <- anova(fitm1, fitall)$`Sum of Sq`[2]
    var_E <- anova(fitall)$`Sum Sq`[2]

    V_exp <- c(var_m1, var_m2, var_Int, var_E, var_tot)

  }else if(Cl_sel == "Meta1"){
    fitm1 <- lm(formula(lm_formula[3]), data = data_g)

    var_E <- anova(fitm1)$`Sum Sq`[2]
    var_m1 <- var_tot - var_E
    V_exp <- c(var_m1, 0, 0, var_E, var_tot)
  }else if(Cl_sel == "Meta2"){
    fitm2 <- lm(formula(lm_formula[4]), data = data_g)

    var_E <- anova(fitm2)$`Sum Sq`[2]
    var_m2 <- var_tot - var_E
    V_exp <- c(0, var_m2, 0, var_E, var_tot)
  }
  names(V_exp) <- c("Meta1","Meta2", "Int", "E", "Tot")
  return(V_exp)

}

#' Calculate proportion of each component in gam model
#' @param eNODAL_obj An eNODAL_obj object
#' @param index A number specify which column to calculate
#' @param Cl_sel A text specify which cluster this variable belongs to
#' @importFrom stats var formula
#' @return A vector contain proportion of variance each component explained.
calcProp_gam <- function(eNODAL_obj, index, Cl_sel){

  X <- eNODAL_obj@Z
  Meta1 <- eNODAL_obj@Meta1
  Meta2 <- eNODAL_obj@Meta2
  Xi <- X[,index]
  V_exp <- rep(0, 5)

  gam_formula <- eNODAL_obj@eNODAL_middle$formula$gam
  data_g <- data.frame(Y = Xi, Meta1, Meta2)

  var_tot <- var(Xi)*(nrow(X) - 1)

  if(Cl_sel == "Add"){

    fita <- mgcv::gam(formula(gam_formula[2]), data = data_g)
    fitm1 <- mgcv::gam(formula(gam_formula[3]), data = data_g)
    fitm2 <- mgcv::gam(formula(gam_formula[4]), data = data_g)

    var_E <- var_tot * (1 - mgcv::anova.gam(fita)$dev.expl)
    var_m1 <- mgcv::anova.gam(fitm2, fita)$Deviance[2]
    var_m2 <- mgcv::anova.gam(fitm1, fita)$Deviance[2]

    V_exp[c(1,2,4,5)] <- c(var_m1, var_m2, var_E, var_tot)

  }else if(Cl_sel == "Int"){

    fitall <- mgcv::gam(formula(gam_formula[1]), data = data_g)
    fita <- mgcv::gam(formula(gam_formula[2]), data = data_g)
    fitm1 <- mgcv::gam(formula(gam_formula[3]), data = data_g)
    fitm2 <- mgcv::gam(formula(gam_formula[4]), data = data_g)

    var_E <- var_tot * (1 - mgcv::anova.gam(fitall)$dev.expl)
    var_Int <- mgcv::anova.gam(fita, fitall)$Deviance[2]
    var_m1 <- mgcv::anova.gam(fitm2, fitall)$Deviance[2]
    var_m2 <- mgcv::anova.gam(fitm1, fitall)$Deviance[2]

    V_exp <- c(var_m1, var_m2, var_Int, var_E, var_tot)

  }else if(Cl_sel == "Meta1"){
    fitm1 <- mgcv::gam(formula(gam_formula[3]), data = data_g)

    var_E <- var_tot * (1 - mgcv::anova.gam(fitm1)$dev.expl)
    var_m1 <- var_tot - var_E
    V_exp <- c(var_m1, 0, 0, var_E, var_tot)
  }else if(Cl_sel == "Meta2"){
    fitm2 <- mgcv::gam(formula(gam_formula[4]), data = data_g)

    var_E <- var_tot * (1 - mgcv::anova.gam(fitm2)$dev.expl)
    var_m2 <- var_tot - var_E
    V_exp <- c(0, var_m2, 0, var_E, var_tot)
  }
  names(V_exp) <- c("Meta1","Meta2", "Int", "E", "Tot")
  return(V_exp)
}

#' Calculate proportion of each component in gam model
#' @param eNODAL_obj An eNODAL_obj object
#' @return An eNODAL_obj object with calculated proportion.
#' @export
calcProp <- function(eNODAL_obj){

  Cl_df <- eNODAL_obj@eNODAL_output$Cluster_res
  if(is.null(Cl_df)){
    message("Please run hypothesis test first! i.e. eNODAL_obj <- runHT(eNODAL_obj)")
    return(eNODAL_obj)
  }
  X <- eNODAL_obj@Z
  Meta1 <- eNODAL_obj@Meta1
  Meta2 <- eNODAL_obj@Meta2
  idx_sig <- which(!is.na(Cl_df$Sig1))
  X_tmp <- X[,idx_sig]

  V_exp <- matrix(0, nrow = ncol(X_tmp), ncol = 5)
  for(i in 1:ncol(X_tmp)){
    Cl_sel <- Cl_df[idx_sig[i],"Sig1"]
    if(Cl_df$model[idx_sig[i]] == "lm"){
      V_tmp <- calcProp_lm(eNODAL_obj, idx_sig[i], Cl_sel)
    }else if(Cl_df$model[idx_sig[i]] == "gam"){
      V_tmp <- calcProp_gam(eNODAL_obj, idx_sig[i], Cl_sel)
    }
    V_exp[i,] <- V_tmp
  }
  colnames(V_exp) <- c("Meta1","Meta2", "Int", "E", "Tot")
  rownames(V_exp) <- colnames(X_tmp)
  eNODAL_obj@eNODAL_middle$VarExplained <- V_exp
  return(eNODAL_obj)
}

#' Calculate interpretable features for eNODAL object
#' @param eNODAL_obj An eNODAL object as input
#' Interaction effect only works for one numeric + one categorical variable
#' @param baseline Baseline for comparison(like "control" in Meta2). By default is NULL.
#' @param cor_method Method for calculate correlation matrix. By default is "spearman".
#' @returns eNODAL_obj object with created interpretable features
#' @importFrom stats formula model.matrix cor na.omit t.test
#' @export
createFeatures <- function(eNODAL_obj, baseline = NULL, cor_method = "spearman"){

  Meta1 <- eNODAL_obj@Meta1
  Meta2 <- eNODAL_obj@Meta2
  Z <- eNODAL_obj@Z

  data_g0 <- data.frame(Meta1, Meta2)
  meta1_type <- eNODAL_obj@eNODAL_middle$type$Meta1
  meta2_type <- eNODAL_obj@eNODAL_middle$type$Meta2

  Int_formula <- substr(eNODAL_obj@eNODAL_middle$formula$linear["int_alt"],2,100)
  Int_term <- model.matrix(formula(Int_formula), data = data_g0)[,-1]

  all_formula <- substr(eNODAL_obj@eNODAL_middle$formula$linear["full"],2,100)
  all_term <- model.matrix(formula(all_formula), data = data_g0)[,-1]
  all_col <- ncol(all_term)
  F_res <- matrix(0, nrow = ncol(Z), ncol = all_col)
  colnames(F_res) <- paste("V", 1:ncol(F_res), sep = "")
  grp <- c()

  # Marginal effect Meta1
  k <- 1
  for(i in 1:length(meta1_type)){
    if(meta1_type[i] == "num"){
      tmp <- cor(Meta1[,i], Z, method = cor_method)
      tmp <- fisher_trans(tmp)
      tmp <- tmp/(1/sqrt(nrow(Z) - 3))
      F_res[,k] <- tmp
      colnames(F_res)[k] <- colnames(Meta1)[i]
      k <- k + 1
      grp = c(grp, "Meta1")
    }else{

      fac_tab <- na.omit(unique(Meta1[,i]))

      if(is.null(baseline)){
        baseline <- fac_tab[1]
      }else if(!(baseline %in% fac_tab)){
        baseline <- fac_tab[1]
      }
      Z_baseline <- Z[Meta1[,i] == baseline,]
      fac_other <- setdiff(fac_tab, baseline)
      for(j in 1:length(fac_other)){
        Z_tmp <- Z[Meta1[,i] == fac_other[j],]
        tmp <- sapply(1:ncol(Z_tmp), function(x){
          test_tmp <- t.test(Z_tmp[,x], Z_baseline[,x])$statistic
        })
        F_res[,k] <- tmp
        colnames(F_res)[k] <- paste0(colnames(Meta1)[i],fac_other[j])
        k <- k + 1
        grp = c(grp, "Meta1")
      }
    }
  }


  # Marginal effect Meta2
  for(i in 1:length(meta2_type)){
    if(meta2_type[i] == "num"){
      tmp <- cor(Meta2[,i], Z, method = cor_method)
      tmp <- fisher_trans(tmp)
      tmp <- tmp/(1/sqrt(nrow(Z) - 3))
      F_res[,k] <- tmp
      colnames(F_res)[k] <- colnames(Meta2)[i]
      k <- k + 1
      grp = c(grp, "Meta2")
    }else{

      fac_tab <- na.omit(unique(Meta2[,i]))

      if(is.null(baseline)){
        baseline <- fac_tab[1]
      }else if(!(baseline %in% fac_tab)){
        baseline <- fac_tab[1]
      }
      Z_baseline <- Z[Meta2[,i] == baseline,]
      fac_other <- setdiff(fac_tab, baseline)
      for(j in 1:length(fac_other)){
        Z_tmp <- Z[Meta2[,i] == fac_other[j],]
        tmp <- sapply(1:ncol(Z_tmp), function(x){
          test_tmp <- t.test(Z_tmp[,x], Z_baseline[,x])$statistic
        })
        F_res[,k] <- tmp
        colnames(F_res)[k] <- paste0(colnames(Meta2)[i],fac_other[j])
        k <- k + 1
        grp = c(grp, "Meta2")
      }
    }
  }

  # Interaction effect

  if(is.null(baseline)){
    baseline <- strsplit(colnames(Int_term)[1],":")[[1]][[2]]
  }
  idx_baseline <- grep(baseline, colnames(Int_term))
  if(!is.null(idx_baseline)){
    baseline <- strsplit(colnames(Int_term)[1],":")[[1]][[2]]
    idx_baseline <- grep(baseline, colnames(Int_term))
  }
  for(j in 1:length(idx_baseline)){
    tmp <- Int_term[,idx_baseline[j]]
    tmp_other <- strsplit(colnames(Int_term)[idx_baseline[j]], ":")[[1]]
    comp_name <- tmp_other[!grepl(baseline, tmp_other)]

    contrast_tmp <- grep(comp_name, colnames(Int_term))
    idx_tmp0 <- which(tmp != 0)
    cor_tmp0 <- cor(tmp[idx_tmp0], Z[idx_tmp0,], method = cor_method)
    fisher_tmp0 <- fisher_trans(cor_tmp0)
#      fisher_tmp00 <- fisher_tmp0/(1/sqrt(length(idx_tmp) - 3))
    contrast_tmp <- setdiff(contrast_tmp, idx_baseline[j])
    for(jj in 1:length(contrast_tmp)){
      tmp1 <- Int_term[,contrast_tmp[jj]]
      idx_tmp1 <- which(tmp1 != 0)
      cor_tmp1 <- cor(tmp1[idx_tmp1], Z[idx_tmp1,], method = cor_method)
      fisher_tmp1 <- fisher_trans(cor_tmp1)
#        fisher_tmp11 <- fisher_tmp1/(1/sqrt(length(idx_tmp1) - 3))
      sN <- sqrt(1/(length(idx_tmp0) - 3) + 1/(length(idx_tmp1) - 3))
      F_res[,k] <- (fisher_tmp1 - fisher_tmp0)/sN
      colnames(F_res)[k] <- colnames(Int_term)[contrast_tmp[jj]]
      k <- k + 1
      grp = c(grp, "Int")
      if(k > ncol(F_res)){
        break
      }
    }
    if(k > ncol(F_res)){
      break
    }
  }

  eNODAL_obj@eNODAL_middle$Features = F_res
  eNODAL_obj@eNODAL_middle$F_type = grp
  return(eNODAL_obj)
}

#' Change names of eNODAL_obj clustering result
#' @param eNODAL_obj An eNODAL_obj object as input
#' @returns eNODAL_obj object with names specified by Meta_name.
#'
.changename <- function(eNODAL_obj){

  Meta1_name <- eNODAL_obj@eNODAL_middle$name$Meta1
  Meta2_name <- eNODAL_obj@eNODAL_middle$name$Meta2

  Cl_df <- eNODAL_obj@eNODAL_output$Cluster_res
  if(is.null(Cl_df)){
    return(eNODAL_obj)
  }else{
    if(!all(is.na(Cl_df$Sig1))){
      Cl_df$Sig1 <- gsub("Meta1", Meta1_name, Cl_df$Sig1)
      Cl_df$Sig1 <- gsub("Meta2", Meta2_name, Cl_df$Sig1)
    }
    if(!all(is.na(Cl_df$Sig2))){
      Cl_df$Sig2 <- gsub("Meta1", Meta1_name, Cl_df$Sig2)
      Cl_df$Sig2 <- gsub("Meta2", Meta2_name, Cl_df$Sig2)
    }
    eNODAL_obj@eNODAL_output$Cluster_res <- Cl_df
  }
  return(eNODAL_obj)

}

#' Print all the clustering of eNODAL result
#' @param eNODAL_obj, an eNODAL object.
#' @return Clustering result of eNODAL object.
#' @export
show_clusters <- function(eNODAL_obj){

  Cl_df <- eNODAL_obj@eNODAL_output$Cluster_res
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
#' @param eNODAL_obj, an eNODAL object.
#' @return Clustering result of eNODAL object.
#' @export
get_clusters <- function(eNODAL_obj){
  Cl_df <- eNODAL_obj@eNODAL_output$Cluster_res
  return(Cl_df)
}

