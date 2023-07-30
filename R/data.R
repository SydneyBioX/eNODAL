#' @title Part of processed mouse liver nutriomics proteomics dataset
#' @description  A list of three type of data:
#' Proteomics data: A dataframe contain 127 mice and 1000 proteome.
#' Nutrition condition: A dataframe of 127 mice and 4 nutrition variables.
#' Drug treatment: A factor contain drug intake for each mouce.
#' @format A list of three type of data:
#' \describe{
#'   \item{Z}{A dataframe of proteome data}
#'   \item{Meta1}{A dataframe of nutrition intake condition in experiment}
#'   \item{Meta2}{A factor of drug treatment}
#' }
#' @usage data(Proteomics, package = 'eNODAL')
"Proteomics"

#' @title An example of fitted eNODAL object
#' @description An example of fitted eNODAL object
#' Generate from following code
#' data(Proteomics)
#' eNODAL_obj <- eNODAL_build(Z = Z, Meta1 = Meta1, Meta2 = Meta2,
#'  sig_test = "lm", test_func = "lm")
#' eNODAL_obj <- runeNODAL(eNODAL_obj)
#' @usage data(eNODAL_example, package = 'eNODAL')
"eNODAL_example"

#' @title Full proteomics mouse nutrition experimental data
#' @description A list of three type of data:
#' @format A list of three type of data:
#' \describe{
#'   \item{Meta}{All phenotypical data}
#'   \item{Prot_raw}{Raw proteomics data}
#'   \item{Prot_RUV}{RUV processed proteomics data}
#' }
#' @usage data(Proteomics_full, package = 'eNODAL')
"Proteomics_full"
