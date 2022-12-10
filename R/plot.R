#' determine whether points is in convex hull
#' @param testpts A dataframe contain testing point
#' @param calpts A dataframe contain point to determine convex hull.
#' @param hull Create hull to test whether testpts in hull.
#'  By default will create by geometry::convhulln(calpts).
#' @param tol Tolerance. By defualt is mean(mean(abs(calpts)))*sqrt(.Machine$double.eps)).
#' @import ggplot2
#' @return An indicator weather each point belongs to the
#'  convex hull spanned by calpts.
.inhull <- function(testpts, calpts,
                    hull=geometry::convhulln(calpts),
                    tol=mean(mean(abs(calpts)))*sqrt(.Machine$double.eps)) {

  calpts <- as.matrix(calpts)
  testpts <- as.matrix(testpts)
  p <- dim(calpts)[2]
  cx <- dim(testpts)[1] # rows in testpts
  nt <- dim(hull)[1] # number of simplexes in hull
  nrmls <- matrix(NA, nt, p)

  degenflag <- matrix(TRUE, nt, 1)

  for (i in 1:nt) {

    nullsp<-t(MASS::Null(t(calpts[hull[i,-1],] -
                             matrix(calpts[hull[i,1],],p-1,p, byrow=TRUE))))

    if (dim(nullsp)[1] == 1){
      nrmls[i,]<-nullsp
      degenflag[i]<-FALSE
    }
  }

  if(length(degenflag[degenflag]) > 0) warning(length(degenflag[degenflag])," degenerate faces in convex hull")
  nrmls <- nrmls[!degenflag,]
  nt <- dim(nrmls)[1]

  center = apply(calpts, 2, mean)
  a <- calpts[hull[!degenflag,1],]

  nrmls <- nrmls/matrix(apply(nrmls, 1, function(x) sqrt(sum(x^2))), nt, p)

  dp <- sign(apply((matrix(center, nt, p, byrow=TRUE)-a) * nrmls, 1, sum))
  nrmls <- nrmls*matrix(dp, nt, p)

  aN <- diag(a %*% t(nrmls))
  val <- apply(testpts %*% t(nrmls) - matrix(aN, cx, nt, byrow=TRUE), 1,min)

  val[abs(val) < tol] <- 0
  as.integer(sign(val))
}

#' Finding convex hull of given data
#' @param x A vector of x-axis coordinates.
#' @param y A vector of y-axis coordinates.
#' @param rgnames Names of the grid.
#' @param res Resolution of the convex hull. By default is 100.
#' @param x.limits Limits of x-axis. By default is NA, will determined by data.
#' @param y.limits Limits of y-axis. By default is NA, will determined by data.
#' @import ggplot2
#' @return A dataframe of data in convex hull.
.findConvex<-function(x, y, rgnames, res=101, x.limits=NA, y.limits=NA){
  hull<-cbind(x,y)[grDevices::chull(cbind(x,y)),]

  # Either use the specifiec limits for the grid, or use pretty
  if(length(which(is.na(x.limits) == T)) > 1){
    px<-base::pretty(x)
  }else{
    px<-x.limits
  }
  if((length(which(is.na(y.limits) == T)) > 1)){
    py<- base::pretty(y)
  }else{
    py<-y.limits
  }

  # Get the matrix
  x.new<-seq(min(px, na.rm=T),max(px, na.rm=T),len=res)
  y.new<-seq(min(py, na.rm=T),max(py, na.rm=T),len=res)
  ingrid<-as.data.frame(base::expand.grid(x.new,y.new))
  Fgrid<-ingrid
  Fgrid[(sp::point.in.polygon(ingrid[,1], ingrid[,2], hull[,1],hull[,2])==0),]<-NA
  names(Fgrid)<-rgnames
  return(Fgrid)
}

#' Function to make nice surfaces using ggplot
#' @param GAM A gam function for prediction
#' @param data A dataframe contain all observed data
#' @param XYZ XY axis used for plot, Z axis will be used as median, can be multiple Z.
#' @param surf_min minimal value of the fitted surface.
#' If it is NA, the minimal value is determined by fitted value.
#' @param surf_max maximal value of the fitted surface.
#' If it is NA, the maximal value is determined by fitted value.
#' @param x.limits Limits of x axis.
#' If it is NA, the limits is determined by input data.
#' @param y.limits Limits of y axis.
#' If it is NA, the limits is determined by input data.
#' @param exclude Point exclude from prediction used in predict.gam. By default is NA.
#' @param lab_size label size in contour.
#' @param nlevels Number of levels in contour. By default is 5.
#' @param contour_at A vector indicates where to draw the contour.
#' @param skip Number of contours to skip.
#' @param axis_size axis.text size.
#' @param axis_lab_size axis.title size.
#' @param predict_val set additional predict parameters for gam fitting.
#'  e.g. predict_val = data.frame(Meta2 = "control")
#' @param z.val Value of z to be plotted. Plot will condition on z.val and fit f(x,y,z=z.val).
#' @param subtitle Subtitle of the figure.
#' @import ggplot2
#' @return A ggplot type object of NGF.
#' @export
ggSurface<-function(GAM, data, XYZ = NA,
                    predict_val=NA, surf_min=NA,
                    surf_max=NA, x.limits=NA, y.limits=NA, z.val=NA,
                    exclude = NA, subtitle="", lab_size=3, nlevels=5,
                    contour_at=NA, skip=0, axis_size=15, axis_lab_size=15){


  # This specifies the color scheme for surface
  rgb.palette<- grDevices::colorRampPalette(c("blue","cyan","yellow","red"),
                                            space="Lab", interpolate="linear")
  map<-rgb.palette(256)

  # List to hold the plots
  plots_list<-list()

  # List for the order of plots
  if(is.null(XYZ)){
    nutrient.order <- GAM$smooth[[1]]$term
  }else{
    nutrient.order <- XYZ
  }

  # Values to predict over, if they are unspecified
  if(is.na(x.limits)[1] == T){
    x.limits<-c(floor(min(data[,nutrient.order[1]])),
                ceiling(max(data[,nutrient.order[1]])))
  }
  if(is.na(y.limits)[1] == T){
    y.limits<-c(floor(min(data[,nutrient.order[2]])),
                ceiling(max(data[,nutrient.order[2]])))
  }

  # If we do not specify values to slice at, use the 25, 50, and 75 % quantile
  if(any(is.na(z.val) == T)){
    z.val<- apply(data[,nutrient.order[-c(1,2)], drop = F], 2, function(x){median(x)})
  }

  # Fitted list to hold some results for later
  x.new<-seq(min(x.limits, na.rm=T), max(x.limits, na.rm=T), len=501)
  y.new<-seq(min(y.limits, na.rm=T), max(y.limits, na.rm=T), len=501)
  z.new<-z.val
  predictors<-as.data.frame(base::expand.grid(x.new, y.new))
  z.new1 <- matrix(rep(z.new, nrow(predictors)), nrow = nrow(predictors), byrow = T)
  predictors <- data.frame(predictors, z.new1)
  names(predictors)<-nutrient.order
  in.poly<-as.numeric(.inhull(predictors,
                              data[,nutrient.order]) != -1)

  # Add the predictors for the additional 'confounders'
  predictors<-cbind(predictors, predict_val)
  predictors<-predictors[-which(in.poly == 0),]

  # Do the predictions
  predictions <-mgcv::predict.gam(GAM, newdata=predictors,
                                 type="response", exclude=exclude,
                                 newdata.guaranteed=T)

  # Get the proedictions for the kth trait
  predictions_k<-predictions

  # Find the min and max values across all predictions
  mn<-surf_min
  mx<-surf_max
  if(is.na(mn)==T){
    mn<-min(predictions_k, na.rm=T)
  }
  if(is.na(mx)==T){
    mx<-max(predictions_k, na.rm=T)
  }
  locs<-(range(predictions_k, na.rm=TRUE) - mn) / (mx-mn) * 256

  plot_data<-predictors
  plot_data$fit<-predictions_k
  plot_data$x<-plot_data[,nutrient.order[1]]
  plot_data$y<-plot_data[,nutrient.order[2]]
  plot_data <- plot_data
  # Set the contour
  if(is.na(contour_at)[1] == T){
    contour_use<-signif((max(predictions_k, na.rm=T)-min(predictions_k, na.rm=T))/nlevels, 1)
  }else{
    contour_use<-contour_at
  }

  # Make the plot
  plot<-ggplot(plot_data, aes(x=x, y=y)) +
    geom_raster(aes(fill=fit), show.legend=F, interpolate=F, na.rm=T) +
    scale_fill_gradientn(colors=map[locs[1]:locs[2]]) +
    geom_contour(data=plot_data, aes(x=x, y=y, z=fit), na.rm=T, color="black", binwidth=contour_use) +
    metR::geom_label_contour(data=plot_data, aes(x=x, y=y, z=fit), size=lab_size, binwidth=contour_use, skip=skip) +
    theme_bw() +
    theme(axis.text=element_text(size=axis_size), axis.title=element_text(size=axis_lab_size)) +
    theme(title=element_text(size=15)) +
    xlim(x.limits) +
    ylim(y.limits)

  return(plot)

}

#' Generate NGF using eNODAL_obj object
#' @param eNODAL An eNODAL_obj object
#' @param Zaxis An index or name of vector to select a variable for Zaxis.
#' @param XYaxis Variables of XYaxis,
#' can be chosen from all continuous variable in meta1 and meta2.
#' If more than one vector provided, will automatically do choose(n,2)
#' @param By A factor used for `by` in gam fitting.
#' If it is NULL, will not fit with by.
#' @param scaled Whether scale for Z variable.
#' @param type Can be chosen from "Int", "Add", "Meta1" and "Meta2",
#' for different types of visualization.
#' By default is NULL, will be chosen from eNODAL_obj fitting result.
#' @param formula0 Formula used for gam fitting. By default is NULL.
#'  Will automatically determined by other parameters.
#' @param ... Other parameters can be passed to NGFPlot.
#' @import ggplot2
#' @importFrom stats median prcomp
#' @seealso \code{\link{ggSurface}}
#' @return A list of NGF plot.
#' @examples
#' data(eNODAL_example)
#' p = NGFPlot(eNODAL_example, 1, c("intake.P", "intake.C"))
#' @export
NGFPlot <- function(eNODAL_obj, Zaxis, XYaxis, By = NULL, scaled = T,
                    type = NULL, formula0 = NULL, ...){

  rgb.palette<- grDevices::colorRampPalette(c("blue","cyan","yellow","red"),
                                            space="Lab", interpolate="linear")
  map<-rgb.palette(256)

  Z <- eNODAL_obj@Z
  Meta1 <- eNODAL_obj@Meta1
  Meta2 <- eNODAL_obj@Meta2
  gam_k <- eNODAL_obj@params$gam_k
  Meta1_name <- eNODAL_obj@eNODAL_middle$name$Meta1
  Meta2_name <- eNODAL_obj@eNODAL_middle$name$Meta2

  Cl_df <- eNODAL_obj@eNODAL_output$Cluster_res
  if(any(grepl(paste0("sig|Int|Add|",Meta1_name,"|",Meta2_name), Zaxis))){
    if(Zaxis == "sig"){
      idx_tmp <- which(Cl_df$Sig0 == Zaxis)
    }else if(Zaxis %in% c("Int","Add",Meta1_name, Meta2_name)){
      idx_tmp <- which(Cl_df$Sig1 == Zaxis)
    }else{
      idx_tmp <- which(Cl_df$Sig2 == Zaxis)
    }
    Z_tmp <- Z[,idx_tmp]
    Z_tmp_med <- apply(Z_tmp, 1, median)
    Z_pca <- prcomp(Z_tmp, center = TRUE, scale. = TRUE)
    Z0 <- Z_pca$x[,1]
    Z0 <- Z0*sign(sum(Z0*Z_tmp_med))
  }else{
    Z0 <- Z[,Zaxis]
  }

  if(scaled){
    Z0 <- scale(Z0)
  }
  data_g <- data.frame(outcome = Z0, Meta1, Meta2)

  if(is.null(formula)){
    GAM <- mgcv::gam(formula0, data = data_g)
  }else{
    if(is.null(By)){
      formula0 <- paste0("outcome~s(",paste0(XYaxis, collapse = ","), ",k=",gam_k,
      ")")
    }else{
      formula0 <- paste0("outcome~s(",paste0(XYaxis, collapse = ","), ",k=",gam_k,
      ",by=", By,")+",By)
    }
    GAM <- mgcv::gam(formula(formula0), data = data_g)
  }
#  if(is.null(type)){
#    type = Cl_df[Zaxis,"Sig1"]
#  }
#  if(!(type %in% c("Int", "Add", "Meta1", "Meta2"))){
#    message("Type must be one of 'Int', 'Add', 'Meta1' and 'Meta2'")
#    return(NULL)
#  }

  p_list <- list()
  if(is.null(By)){
    if(length(XYaxis) == 1){
      name_sel <- colnames(Z)
      if(is.numeric(Zaxis)){
        name_sel <- name_sel[Zaxis]
      }else{
        name_sel <- Zaxis
      }
      p <- ggplot(data_g) + geom_point(aes_string(x = XYaxis, y = "outcome")) +
        labs(x = XYaxis, y = name_sel) +
        theme_bw()
      p_list[[1]] <- p
    }else if(length(XYaxis) == 2){
      p <- ggSurface(GAM = GAM, data = data_g, XYZ = XYaxis) +
        labs(x = XYaxis[1], y = XYaxis[2])
      p_list[[1]] <- p
    }else{
      comb <- utils::combn(XYaxis, 2)
      surf_min<-NA
      surf_max<-NA
      ranges_list<-list()
      for(i in 1:ncol(comb)){
        XY_tmp <- comb[,i]
        Z_tmp <- setdiff(XYaxis,XY_tmp)
        XYZ_tmp <- c(XY_tmp, Z_tmp)
        p_tmp <- ggSurface(GAM = GAM, data = data_g, XYZ = XYZ_tmp) +
          labs(x = XY_tmp[1], y = XY_tmp[2])
        surf_min<-c(surf_min, min(p_tmp$data$fit))
        surf_max<-c(surf_max, max(p_tmp$data$fit))
        ranges_list[[i]]<-c(range(pretty(range(p_tmp$data[,1]))),
                            range(pretty(range(p_tmp$data[,2]))))
        p_list[[i]] <- p_tmp
      }

      surf_min<-min(surf_min, na.rm=T)
      if(surf_min > 0){
        surf_min<-surf_min * 0.98
      }
      if(surf_min < 0){
        surf_min<-surf_min * 1.02
      }
      surf_max<-max(surf_max, na.rm=T)
      if(surf_max > 0){
        surf_max<-surf_max * 1.02
      }
      if(surf_max < 0){
        surf_max<-surf_max * 0.98
      }

      for(i in 1:length(p_list)){
        predictions_k <- p_list[[i]]$data$fit
        locs<-(range(predictions_k, na.rm=TRUE) - surf_min) / (surf_max-surf_min) * 256
        p_list[[i]] <-suppressMessages(p_list[[i]] + xlim(ranges_list[[i]][c(1,2)]) +
          ylim(ranges_list[[i]][c(3,4)]) +
          scale_fill_gradientn(colors=map[locs[1]:locs[2]]))

      }
    }
    return(p_list)
  }else{
    groups <- sort(unique(data_g[,By]))
    if(length(XYaxis) == 1){
      name_sel <- colnames(Z)
      if(is.numeric(Zaxis)){
        name_sel <- name_sel[Zaxis]
      }else{
        name_sel <- Zaxis
      }
      p <- ggplot(data_g) + geom_point(aes_string(x = XYaxis, y = "outcome",
                                                  color = By)) +
        labs(x = XYaxis, y = name_sel) + theme_bw()
      p_list[[1]] <- p
    }else if(length(XYaxis) == 2){
      surf_min<-NA
      surf_max<-NA
      ranges_list<-list()
      for(j in 1:length(groups)){
        data_gj <- data_g[data_g[,By] == groups[j],]
        df_tmp <- data.frame(By = groups[j])
        colnames(df_tmp) <- By
        GAM <- GAM
        p_tmp <- ggSurface(GAM = GAM, data = data_gj, XYZ = XYaxis,
                           predict_val=df_tmp) + labs(x = XYaxis[1], y =XYaxis[2], tag = groups[j])
        surf_min<-c(surf_min, min(p_tmp$data$fit))
        surf_max<-c(surf_max, max(p_tmp$data$fit))
        ranges_list[[j]]<-c(range(pretty(range(p_tmp$data[,1]))),
                            range(pretty(range(p_tmp$data[,2]))))
        p_list[[j]] <- p_tmp
      }

      surf_min<-min(surf_min, na.rm=T)
      if(surf_min > 0){
        surf_min<-surf_min * 0.98
      }
      if(surf_min < 0){
        surf_min<-surf_min * 1.02
      }
      surf_max<-max(surf_max, na.rm=T)
      if(surf_max > 0){
        surf_max<-surf_max * 1.02
      }
      if(surf_max < 0){
        surf_max<-surf_max * 0.98
      }

      surf_min <- surf_min
      surf_max <- surf_max
      for(i in 1:length(p_list)){
        predictions_k <- p_list[[i]]$data$fit
        locs<-(range(predictions_k, na.rm=TRUE) - surf_min) / (surf_max-surf_min) * 256
        p_list[[i]] <- suppressMessages(p_list[[i]] + xlim(ranges_list[[i]][c(1,2)]) +
          ylim(ranges_list[[i]][c(3,4)]) + scale_fill_gradientn(colors=map[locs[1]:locs[2]]))
      }

    }else{

      comb <- utils::combn(XYaxis, 2)
      surf_min<-NA
      surf_max<-NA
      ranges_list<-list()

      for(i in 1:ncol(comb)){
        XY_tmp <- comb[,i]
        Z_tmp <- setdiff(XYaxis,XY_tmp)
        XYZ_tmp <- c(XY_tmp, Z_tmp)
        for(j in 1:length(groups)){
          data_gj <- data_g[data_g[,By] == groups[j],]
          df_tmp <- data.frame(By = groups[j])
          colnames(df_tmp) <- By
          p_tmp <- ggSurface(GAM = GAM, data = data_gj, XYZ = XYZ_tmp,
                             predict_val=df_tmp) + labs(x = XYZ_tmp[1], y = XYZ_tmp[2], tag = groups[j])
          surf_min<-c(surf_min, min(p_tmp$data$fit))
          surf_max<-c(surf_max, max(p_tmp$data$fit))
          ranges_list[[(i-1)*length(groups) + j]]<-c(range(pretty(range(p_tmp$data[,1]))),
                                                 range(pretty(range(p_tmp$data[,2]))))
          p_list[[(i-1)*length(groups) + j]] <- p_tmp
        }

      }

      surf_min<-min(surf_min, na.rm=T)
      if(surf_min > 0){
        surf_min<-surf_min * 0.98
      }
      if(surf_min < 0){
        surf_min<-surf_min * 1.02
      }
      surf_max<-max(surf_max, na.rm=T)
      if(surf_max > 0){
        surf_max<-surf_max * 1.02
      }
      if(surf_max < 0){
        surf_max<-surf_max * 0.98
      }

      for(i in 1:length(p_list)){
        predictions_k <- p_list[[i]]$data$fit
        locs<-(range(predictions_k, na.rm=TRUE) - surf_min) / (surf_max-surf_min) * 256
        p_list[[i]] <- suppressMessages(p_list[[i]] + xlim(ranges_list[[i]][c(1,2)]) +
          ylim(ranges_list[[i]][c(3,4)]) + scale_fill_gradientn(colors=map[locs[1]:locs[2]]))
      }

    }
    return(p_list)
  }
}

#' Generate clustering result plot for eNODAL_obj object
#' @param eNODAL An eNODAL_obj object.(Already run two stage clustering)
#' @import ggplot2
#' @return A concentric plot of eNODAL_obj clustering result
#' @export
ClusterPlot <- function(eNODAL_obj){

  Cl_df <- eNODAL_obj@eNODAL_output$Cluster_res

  if(is.null(Cl_df)){
    message("Please run processing first!")
    return(NULL)
  }

  plt_df <- data.frame(name = "", value = 0, level = 0, type = NA)
  if(!is.null(Cl_df$Sig0)){

    summary_sig0 <- table(Cl_df$Sig0)
    tmp_df <- data.frame(name = names(summary_sig0), value = as.vector(summary_sig0),
                         level = 1, type = names(summary_sig0))
    plt_df <- rbind(plt_df, tmp_df)

  }

  if(!is.null(Cl_df$Sig1)){
    summary_sig1 <- table(Cl_df$Sig1)
    tmp_df <- data.frame(name = names(summary_sig1), value = as.vector(summary_sig1),
                         level = 2, type = names(summary_sig1))
    plt_df <- rbind(plt_df, tmp_df)
  }

  if(!is.null(Cl_df$Sig2)){
    summary_sig2 <- table(Cl_df$Sig2)
    type_sig2 <- sapply(strsplit(names(summary_sig2), "_"), `[[`, 1)
    tmp_df <- data.frame(name = names(summary_sig2), value = as.vector(summary_sig2),
                         level = 3, type = type_sig2)
    plt_df <- rbind(plt_df, tmp_df)
  }

  plt_df$type <- as.factor(plt_df$type)
  plt_df$level <- factor(plt_df$level)
  p <- ggplot(plt_df, aes(x = level, y = value, fill = type, alpha = level)) +
    geom_col(width = 1, color = "gray90", size = 0.25, position = position_stack()) +
    coord_polar(theta = "y") +
    geom_text(aes(label = name), size = 2.5, position = position_stack(vjust = 0.5)) +
    scale_x_discrete(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    scale_alpha_manual(values = c("0" = 0, "1" = 1, "2" = 1, "3" = 0.7), guide = "none") +
    scale_fill_discrete(na.translate = F) +
#    scale_fill_manual(values = c("gray80", "#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854"), na.translate = F) +
    labs(x = NULL, y = NULL, fill = "Group") +
    theme_minimal()
  return(p)
}
#' Generate Boxplot using eNODAL_obj object
#' @param eNODAL_obj An eNODAL_obj object
#' @param Zaxis An index or name of vector to select a variable for Zaxis.
#' @param XYaxis Variables of XYaxis,
#' can be chosen from all continuous variable in meta1 and meta2.
#' If more than one vector provided, will automatically do choose(n,2)
#' @param scaled Whether scale for Z variable.
#' @param ... Other parameters can be passed to NGFPlot.
#' @seealso \code{\link{ggSurface}}
#' @return A list of NGF plot.
#' @import ggplot2
#' @importFrom stats median prcomp
#' @export
Boxplot <- function(eNODAL_obj, XYaxis, Zaxis, scaled = F, ...){

  Cl_df <- eNODAL_obj@eNODAL_output$Cluster_res

  if(is.null(Cl_df)){
    message("Please run processing first!")
    return(NULL)
  }
  Z <- eNODAL_obj@Z
  Meta1 <- eNODAL_obj@Meta1
  Meta2 <- eNODAL_obj@Meta2
  gam_k <- eNODAL_obj@params$gam_k
  Meta1_name <- eNODAL_obj@eNODAL_middle$name$Meta1
  Meta2_name <- eNODAL_obj@eNODAL_middle$name$Meta2

  if(any(grepl(paste0("sig|Int|Add|",Meta1_name,"|",Meta2_name), Zaxis))){
    if(Zaxis == "sig"){
      idx_tmp <- which(Cl_df$Sig0 == Zaxis)
    }else if(Zaxis %in% c("Int","Add",Meta1_name, Meta2_name)){
      idx_tmp <- which(Cl_df$Sig1 == Zaxis)
    }else{
      idx_tmp <- which(Cl_df$Sig2 == Zaxis)
    }
    Z_tmp <- Z[,idx_tmp]
    Z_tmp_med <- apply(Z_tmp, 1, median)
    Z_pca <- prcomp(Z_tmp, center = TRUE, scale. = TRUE)
    Z0 <- Z_pca$x[,1]
    Z0 <- Z0*sign(sum(Z0*Z_tmp_med))
  }else{
    Z0 <- Z[,Zaxis]
  }

  if(scaled){
    Z0 <- scale(Z0)
  }
  data_g <- data.frame(outcome = Z0, Meta1, Meta2)


  if(length(XYaxis) == 1){
    p <- ggplot(data_g) +
      geom_boxplot(aes_string(x = XYaxis[1], y = "outcome")) +
      labs(x = XYaxis[1], y = Zaxis) + theme_bw()
  }else if(length(XYaxis) == 1){
    p <- ggplot(data_g) +
      geom_boxplot(aes_string(x = XYaxis[1], y = "outcome", fill = XYaxis[2])) +
      labs(x = XYaxis[1], y = Zaxis) + theme_bw()
  }

  return(p)
}
