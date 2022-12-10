#' Run hypothesis testing using with eNODAL object.
#' @param eNODAL_obj Input eNODAL object.
#' @param Boot Number of Bootstrapping used in test (npntest). By default is 50.
#' @return eNODAL object with clustering result(Stage I).
#' @importFrom stats lm p.adjust formula model.matrix
#' @export
runHT <- function(eNODAL_obj, Boot = 50){

  h_adj <- eNODAL_obj@params$h_adj
  adj_method <- eNODAL_obj@params$adj_method
  sig_level <- eNODAL_obj@params$sig_level
  adapt <- eNODAL_obj@params$adapt
  test_func <- eNODAL_obj@params$test_func
  test_method <- eNODAL_obj@params$test_method
  sig_test <- eNODAL_obj@params$sig_test

  Z <- eNODAL_obj@Z
  Meta1 <- eNODAL_obj@Meta1
  Meta2 <- eNODAL_obj@Meta2
  P_mat <- matrix(NA, nrow = ncol(Z), ncol = 5)

  lm_formula <- eNODAL_obj@eNODAL_middle$formula$linear
  gam_formula <- eNODAL_obj@eNODAL_middle$formula$gam

  p_sig <- c()

  for(i in 1:ncol(Z)){
    if(sig_test == "LC-test"){
      dat <- model.matrix(formula(substr(lm_formula[1],
                                         2, nchar(lm_formula[1]))),
                          data.frame(Meta1, Meta2))[,-1]

      p_sig[i] <- LC_test0(dat, Z[,i], n = Boot)
      sig0 <- sig_level["LC"]
    }else if(sig_test == "gam"){
      dat <- data.frame(Y = Z[,i], Meta1, Meta2)
      calc0 <- "Y~1"
      calc1 <- gam_formula[1]
      p_sig[i] <- .nested_test(dat, calc0, calc1,
                            func = mgcv::gam, method = test_method)$p
      sig0 <- sig_level["Sig"]
    }else if(sig_test == "lm"){

      dat <- data.frame(Y = Z[,i], Meta1, Meta2)
      calc0 <- "Y~1"
      calc1 <- lm_formula[1]
      p_sig[i] <- .nested_test(dat, calc0, calc1,
                            func = lm, method = test_method)$p
      sig0 <- sig_level["Sig"]
    }
  }
  p_sig_adj <- p.adjust(p_sig, method = adj_method)

  if(adj_method == "none"){
    p_sig_adj <- p_sig
  }else{
    p_sig_adj <- p.adjust(p_sig, method = adj_method)
  }
  P_mat[,1] <- p_sig_adj

  Cl <- rep("sig", ncol(Z))
  Cl[p_sig_adj > sig0] <- "non-sig"
  Cl_df <- data.frame(Sig0 = Cl)
  idx_sig <- which(p_sig_adj <= sig0)
  Z_tmp <- Z[,idx_sig]
  Cl_df_model <- rep(NA, ncol(Z))

  if(adapt){
    p_lin <- hybrid_test(eNODAL_obj, "linear", test_func, Boot)
    p_lin <- sapply(p_lin,`[[`,1)
    P_mat[idx_sig,2] <- p_lin

    Cl_df_lin <- rep("lm", ncol(Z_tmp))
    Cl_df_lin[p_lin <= sig_level["Linear"]] <- "gam"
    Cl_df_model[idx_sig] <- Cl_df_lin
    Cl_df$model = Cl_df_model
  }else{
    Cl_df_model[idx_sig] <- test_func
    Cl_df$model<- Cl_df_model
  }
  test_func_list <- Cl_df$model[idx_sig]

  p_int_tmp <- rep(NA, ncol(Z_tmp))
  p_Meta1_tmp <- rep(NA, ncol(Z_tmp))
  p_Meta2_tmp <- rep(NA, ncol(Z_tmp))
  Rsq_res <- matrix(0, nrow = ncol(Z_tmp), ncol = 6)

  for(i in 1:ncol(Z_tmp)){

    # Test interaction
    int_res <- .hybrid_test(eNODAL_obj, index = idx_sig[i],
                            type = "interaction",
                            test_func = test_func_list[i])

    p_int_tmp[i] <- int_res$p
    Rsq_res[i,1:2] <- int_res$Rsq

    # Test N
    Meta1_res <- .hybrid_test(eNODAL_obj, index = idx_sig[i],
                              type = "meta1",
                              test_func = test_func_list[i])

    p_Meta1_tmp[i] <- Meta1_res$p
    Rsq_res[i,3:4] <- Meta1_res$Rsq

    # Test T
    Meta2_res <- .hybrid_test(eNODAL_obj, index = idx_sig[i],
                              type = "meta2",
                              test_func = test_func_list[i])

    p_Meta2_tmp[i] <- Meta2_res$p
    Rsq_res[i,5:6] <- Meta2_res$Rsq

  }

  if(h_adj){
    q <- sapply(1:4, function(i){
      ncol(model.matrix(formula(substr(lm_formula[i],
                                       2, nchar(lm_formula[i]))),
                        data.frame(Meta1, Meta2))) - 1
    })
    p_int_tmp <- p_int_tmp * q[1]/(q[1] - q[2])
    p_Meta2_tmp <- p_Meta2_tmp * q[2]/q[4]
    p_Meta1_tmp <- p_Meta1_tmp * q[2]/q[3]
  }
  P_mat[idx_sig,3] = p_int_tmp
  P_mat[idx_sig,4] = p_Meta1_tmp
  P_mat[idx_sig,5] = p_Meta2_tmp

  Cl_int_tmp <- rep("Add", ncol(Z_tmp))
  Cl_int_tmp[p_int_tmp < sig_level["Interaction"]] = "Int"
  names(Cl_int_tmp) <- colnames(Z_tmp)

  Cl_tmp <- rep(NA, ncol(Z))

  idx_add <- which(Cl_int_tmp == "Add")
  Cl_tmp1 <- Cl_int_tmp[idx_add]
  p_Meta2_tmp1 <- p_Meta2_tmp[idx_add]
  p_Meta1_tmp1 <- p_Meta1_tmp[idx_add]
  Cl_tmp1[(p_Meta2_tmp1 <= sig_level["Meta2"])&(p_Meta1_tmp1 > sig_level["Meta1"])] = "Meta2"
  Cl_tmp1[(p_Meta2_tmp1 > sig_level["Meta2"])&(p_Meta1_tmp1 <= sig_level["Meta1"])] = "Meta1"
  Cl_int_tmp[idx_add] <- Cl_tmp1

  Cl_tmp[idx_sig] <- Cl_int_tmp
  Cl_df$Sig1 <- Cl_tmp

  Cl_df$varname <- colnames(Z)
  colnames(P_mat) <- c("Sig", "Linear", "Interaction", "Meta1", "Meta2")
  Rsq_res <- Rsq_res[,c(2,1,5,3)]
  colnames(Rsq_res) <- c("Interaction", "Add", "Meta1", "Meta2")
  rownames(Rsq_res) <- colnames(Z_tmp)
  rownames(P_mat) <- colnames(Z)
  rownames(Cl_df) <- colnames(Z)
  eNODAL_obj@eNODAL_output$Cluster_res <- Cl_df
  eNODAL_obj@eNODAL_output$Pvalue <- P_mat
  eNODAL_obj@eNODAL_output$Rsq_res <- Rsq_res
  return(eNODAL_obj)
}

#' nonparametric ANOVA test from Zhou's paper using bootstrap.
#' @param data_g Matrix or data frame for model fitting.
#' @param calc0 A string of formula for null hypothesis
#' @param calc1 A string of formula for alternative hypothesis
#' @param func0 Fitting function used for null hypothesis.
#' @param func1 Fitting function used for alternative hypothesis.
#' @param Boot Number of bootstrapping. By default is 50.
#' @importFrom stats formula lm
#' @return A list contain bootstrapping pvalue and fitted R.square
#'
.npntest <- function(data_g, calc0, calc1, func0 = lm,
                     func1 = mgcv::gam, Boot = 50){

  data_g <- as.data.frame(data_g)
#  data_g1 <- as.data.frame(data_g1)

  n <- nrow(data_g)
  M0 <- func0(formula = formula(calc0), data = data_g)
  M1 <- func1(formula = formula(calc1), data = data_g)
  E0 <- M0$residuals
  E1 <- M1$residuals
  Y_pred0 <- M0$fitted.values
  flag = 0
  if(sum(E0^2) < sum(E1^2)){
    flag = 1
    F0 <- (sum(E1^2) - sum(E0^2))/sum(E0^2)
  }else{
    F0 <- (sum(E0^2) - sum(E1^2))/sum(E1^2)
  }

  F_star <- c()
  for(b in 1:Boot){
    E1_star <- sample(E1, n, replace = T)
    Z_star <- Y_pred0 + E1_star
    data0_star <- data_g
    data0_star$Y <- Z_star
    M0_star <- func0(formula = formula(calc0), data = data0_star)

    data1_star <- data_g
    data1_star$Y <-Z_star

    M1_star <- func1(formula = formula(calc1), data = data1_star)
    E0_star <- M0_star$residuals
    E1_star <- M1_star$residuals

    if(flag){
      F_star[b] <- (sum(E1_star^2) - sum(E0_star^2))/sum(E0_star^2)
    }else{
      F_star[b] <- (sum(E0_star^2) - sum(E1_star^2))/sum(E1_star^2)
    }

  }
  p_np <- sum(F_star > F0)/Boot

  if("gam" %in% class(M0)){
    Rsq0 <- summary(M0)$r.sq
  }else{
    Rsq0 <- summary(M0)$adj.r.squared
  }

  if("gam" %in% class(M1)){
    Rsq1 <- summary(M1)$r.sq
  }else{
    Rsq1 <- summary(M1)$adj.r.squared
  }

  return(list(p = p_np, Rsq = c(Rsq0, Rsq1)))
}

#' Nested test.
#' @param data_g Matrix or data frame for model fitting.
#' @param calc0 A string of formula for null hypothesis
#' @param calc1 A string of formula for alternative hypothesis
#' @param func Fitting function used for test. Can be chosen from lm and gam.
#' By default is lm.
#' @param method Testing method.
#' Can be chosen from "F", "globaltest", "Tmax", "Chisq", "Cp". By default is "F"
#' @importFrom stats formula anova na.omit coef
#' @return A list contain bootstrapping pvalue and fitted R.square
#'
.nested_test <- function(data_g, calc0, calc1, func = lm, method = "F"){

  data_g <- as.data.frame(data_g)
  M0 <- func(formula(calc0), data = data_g)
  M1 <- func(formula(calc1), data = data_g)

  if("gam" %in% class(M0)){
    Rsq0 <- summary(M0)$r.sq
  }else{
    Rsq0 <- summary(M0)$adj.r.squared
  }

  if("gam" %in% class(M1)){
    Rsq1 <- summary(M1)$r.sq
  }else{
    Rsq1 <- summary(M1)$adj.r.squared
  }
  if(method %in% c("F", "Chisq", "Cp")){
    anova_res <- anova(M0, M1, test = "F")
    #  print(anova_res)
    p <- na.omit(anova_res$`Pr(>F)`)
  }else if(method == "globaltest"){
    gt_res <- globaltest::gt(formula(calc0), formula(calc1), data = data_g)
    p <- gt_res@result[1]
  }else if(method == "Tmax"){
    beta0 <- coef(M0)
    beta1 <- coef(M1)
    var_list <- setdiff(names(beta1), names(beta0))
    tmp <- summary(M1)$coefficients
    tmp <- tmp[var_list,]
    p <- tmp[,4]
    p <- min(p)
  }

  return(list(p = p, Rsq = c(Rsq0, Rsq1)))
}

#' Nested test for linear model.
#' @param Y An omics feature as output.
#' @param Meta1 A dataframe of first set of meta variables.
#' @param Meta2 A dataframe of second set of meta variables.
#' @param calc0 A string of formula for null hypothesis.
#' @param calc1 A string of formula for alternative hypothesis.
#' @param method Testing method.
#' Can be chosen from "F", "globaltest", "Tmax", "Chisq", "Cp". By default is "F".
#' @importFrom stats lm
#' @return A list contain bootstrapping pvalue and fitted R.square
#' @export

lm_test <- function(Y, Meta1, Meta2, calc0, calc1, method = "F"){

  data_g <- data.frame(Y, Meta1, Meta2)
  p <- .nested_test(data_g, calc0, calc1, lm)
  return(p)
}

#' Hybrid test for eNODAL object.
#' @param eNODAL_obj An eNODAL object as output.
#' @param index Number of which omics feature to be tested.
#' @param type Testing type. Can be chosen from
#' "linear"(linear vs gam using npntest), "interaction"(interaction effect)
#' "meta1"(marginal effect of Meta1) and "meta2"(marginal effect of Meta2)
#' @param test_func method for fitting the model, can be chosen from "gam" and "lm".
#' @param Boot Number of bootstrapping. By default is 50.
#' @importFrom stats lm
#' @return pvalue number.

.hybrid_test <- function(eNODAL_obj, index, type = "interaction",
                         test_func = "lm", Boot = 50){

  X <- eNODAL_obj@Z
  Y <- X[,index]
  Meta1 <- eNODAL_obj@Meta1
  Meta2 <- eNODAL_obj@Meta2
  data_g <- data.frame(Y = Y, Meta1, Meta2)
  method = eNODAL_obj@params$test_method
  if(type == "linear"){
    calc0 = eNODAL_obj@eNODAL_middle$formula$linear[1]
    calc1 = eNODAL_obj@eNODAL_middle$formula$gam[1]
    p <- .npntest(data_g, calc0, calc1, lm, mgcv::gam, Boot)
    return(p)
  }
  if(test_func == "lm"){
    func <- lm
    lm_formula <- eNODAL_obj@eNODAL_middle$formula$linear
  }else if(test_func == "gam"){
    func <- mgcv::gam
    gam_formula <- eNODAL_obj@eNODAL_middle$formula$gam
  }

  if(type == "interaction"){
    if(test_func == "lm"){
      calc0 <- lm_formula["int_null"]
      calc1 <- lm_formula["full"]
    }else if(test_func == "gam"){
      calc0 <- gam_formula["int_null"]
      calc1 <- gam_formula["full"]
    }
  }else if(type == "meta1"){
    if(test_func == "lm"){
      calc0 <- lm_formula["meta2"]
      calc1 <- lm_formula["int_null"]
    }else if(test_func == "gam"){
      calc0 <- gam_formula["meta2"]
      calc1 <- gam_formula["int_null"]
    }
  }else if(type == "meta2"){
    if(test_func == "lm"){
      calc0 <- lm_formula["meta1"]
      calc1 <- lm_formula["int_null"]
    }else if(test_func == "gam"){
      calc0 <- gam_formula["meta1"]
      calc1 <- gam_formula["int_null"]
    }
  }
  p <- .nested_test(data_g, calc0, calc1, func, method)
  if(!length(p$p)){
    p$p <- 1
  }
  return(p)
}

#' Hybrid test for eNODAL_obj object.
#' @param eNODAL_obj An eNODAL object as output.
#' @param type Testing type. Can be chosen from
#' "linear"(linear vs gam using npntest), "interaction"(interaction effect)
#' "meta1"(marginal effect of Meta1) and "meta2"(marginal effect of Meta2)
#' @param test_func method for fitting the model, can be chosen from "gam" and "lm".
#' @param Boot Number of bootstrapping. By default is 50.
#' @return pvalue number.
#' @export
hybrid_test <- function(eNODAL_obj, type = "interaction",
                        test_func = "lm", Boot = 50){
  p <- ncol(eNODAL_obj@Z)
  pvalue <- c()
  for(i in 1:p){
    p_value[i] <- .hybrid_test(eNODAL_obj, i, type, test_func, Boot)
  }
  return(p_value)
}

