#' check validity
#' @param object A eNODAL object.
check_validity <- function(object){

  n_Z = nrow(object@Z)
  n_M1 = nrow(object@Meta1)
  n_M2 = nrow(object@Meta2)

  if(!all(n_M1 == n_M2)){
    errors <- c("Numbers of two meta data are not equal!")
    return(errors)
  }
  if(!all(n_Z == n_M1)){
    errors <- c("Numbers of samples in omics data not equal
                to meta data!")
    return(errors)
  }

  if(!is.null(object@Phenotype)){
    n_P = nrow(object@Phenotype)

    if(!all(n_Z == n_P)){
      errors <- c("Numbers of samples in omics data not equal
                to metabolomic response!")
      return(errors)
    }
  }

  return(TRUE)
}

#' @importClassesFrom S4Vectors DataFrame DFrame
setClassUnion("data.frameORNULL", c("data.frame","matrix", "vector", "NULL"))

setClassUnion("characterORNULL", c("character", "NULL"))

setClassUnion("listORNULL", c("list", "NULL"))

#' eNODAL class
#' @slot Z A dataframe omics data, rows are samples.
#' @slot Meta1 A dataframe/matrix/vector of nutrition intake, rows are samples.
#' @slot Meta2 A dataframe/matrix/vector of second meta data, rows are samples.
#' @slot Phenotype A dataframe/matrix/vector of phenotype data, rows are samples.
#' @slot params A list of parameters for hypothesis testing in eNODAL.
#' @slot eNODAL_middle A list contain middle result of eNODAL.
#' Including possible formula used in GAM/Linear model($formula).
#' @slot eNODAL_output A list eNODAL fitting result.
#' @importFrom methods new

eNODAL <- setClass("eNODAL",
                  slots = c(Z = "data.frameORNULL",
                            Meta1 = "data.frameORNULL",
                            Meta2 = "data.frameORNULL",
                            Phenotype = "data.frameORNULL",
                            params = "listORNULL",
                            eNODAL_middle = "list",
                            eNODAL_output = "list"
                  ),
                  prototype = list(
                    params = list()
                  ))

#' @importFrom S4Vectors coolcat
#'
setMethod("show", "eNODAL", function(object) {

  p = length(object@Z)

  cat("Omics dim: ", dim(object@Z), '\n')
  show_eNODAL(object@Meta1, "Meta data1")
  show_eNODAL(object@Meta2, "Meta data2")
  if(!is.null(object@Phenotype)){
    cat("Phenotype data dim:", dim(object@Phenotype), "\n")
  }

  S4Vectors::coolcat("eNODAL output names (%d): %s\n",
                     names(object@eNODAL_output))

  S4Vectors::coolcat("params (%d): %s\n", names(object@params))

})

show_eNODAL <- function(X, header = "Meta data1"){

  if("numeric" %in% class(X)){
    cat(header, "is a vector, mean = ", mean(X), "\n")
  }else if("factor" %in% class(X)){
    cat(header,"is a factor of",  length(levels(X)), "levels = ", levels(X), "\n")
  }else{
    cat(header, "dim:", dim(X), "\n")
  }

}

#' @importFrom S4Vectors setValidity2
#'
S4Vectors::setValidity2("eNODAL", check_validity)
