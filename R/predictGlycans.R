#' Function to predict glycan masses and m/z values
#'
#' This function will predict all possible glycan molecules within the constraints of the parameters.
#' @param dp1 Minimum degree of polymerisation. Requires an integer. 
#' @param dp2 Maximum degree of polymerisation. Requires an integer.
#' @param ESI_mode ESI mode used. Accepts 'pos' or 'neg'
#' @param scan_range1 Lower end of scan range. Requires an integer.
#' @param scan_range2 Upper end of scan range. Requires an integer.
#' @param pent_option Should pentose monomers be included? 1 = yes, 0 = no. Defaults to no.
#' @param modifications Modifications to be considered. Any combination of 'carboxyl', 'phosphate', 'deoxy', 'nacetyl', 'omethyl', 'anhydrobridge', 'oacetyl', 'unsaturated', 'alditol', 'amino', 'dehydrated', 'sulphate' or 'all' or 'none' (default)
#' @param nmod_max Maximum number of modifications per monomer on average (default 1). Does not take into account unsaturated, alditol or dehydrated.
#' @param double_sulphate Can monomers be double-sulphated: 0 for no (default), 1 for yes. For yes you MUST give a value of at least 2 to nmod_max.
#' @param label Are sugars labelled? Currently only accepts none (default, leave out) or 'procainamide'.
#' @export
#' @examples
#' mz_predicted <- predictGlycans(dp1 = 1, dp2 = 8, ESI_mode = 'neg', scan_range1 = 175, scan_range2 = 1400, pent_option = 1, modifications = c("sulphate", "deoxy", "carboxyl"),nmod_max = 2, double_sulphate = 1)

predictGlycans <- function(dp1,dp2,ESI_mode, scan_range1, scan_range2,
                          pent_option=NULL,modifications=NULL,
                          nmod_max=NULL,double_sulphate=NULL,label=NULL){
  path <- paste(system.file(package="glycanPredict"), "sugarMassesPredict-r.py", sep="/")
  #check if pandas installed
  if(!reticulate::py_module_available("pandas")){
    reticulate::py_install("pandas")
  }
  if(!reticulate::py_module_available("numpy")){
    reticulate::py_install("numpy")
  }
  reticulate::source_python(path)
  dp1 = as.integer(dp1)
  dp2 = as.integer(dp2)
  scan_range1 = as.integer(scan_range1)
  scan_range2 = as.integer(scan_range2)
  if(!is.null(pent_option)){
    pent_option = as.integer(pent_option)
  } else{
    pent_option = as.integer(0)
  }
  if(!is.null(nmod_max)){
    nmod_max = as.integer(nmod_max)
  } else{
    nmod_max = as.integer(1)
  }
  if(!is.null(double_sulphate)){
    double_sulphate = as.integer(double_sulphate)
  } else{
    double_sulphate = as.integer(0)
  }
  if(!is.null(modifications)){
    if(is.vector(modifications)){
      modifications = as.list(modifications)
    }
  } else{
    modifications = "none"
  }
  if(is.null(label)){
    label = "none"
  }
  df <- predict_sugars(dp1 = dp1, dp2 = dp2, ESI_mode = ESI_mode,
                       scan_range1 = scan_range1, scan_range2 = scan_range2,
                       pent_option = pent_option, modifications = modifications,
                       label = label, nmod_max = nmod_max, 
                       double_sulphate = double_sulphate)
  return(df)
}
