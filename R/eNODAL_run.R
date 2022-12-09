#' Run eNODAL using with eNODAL object.
#' @param eNODAL Input eNODAL object.
#' @return eNODAL with clustering result.
#' @examples
#' data(Proteomics)
#' Z = Proteomics$Z
#' Meta1 = Proteomics$Meta1
#' Meta2 = Proteomics$Meta2
#' eNODAL_obj = eNODAL_build(Z, Meta1, Meta2, sig_test = "lm", test_func = "lm")
#' eNODAL_obj = runeNODAL(eNODAL_obj)
#' @export
runeNODAL <- function(eNODAL){
  eNODAL <- runHT(eNODAL)
  eNODAL <- runSubclust(eNODAL)
  eNODAL <- .changename(eNODAL)
  return(eNODAL)
}
