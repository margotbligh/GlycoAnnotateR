#' Function to predict glycan masses and m/z values
#'
#' @include setClass.R
#'
#' @description This function will predict all possible glycan molecules within the constraints of the parameters.
#' @param param A \code{predictGlycansParam} object containing all parameters for the
#'     predictGlycans function. See \link[glycanPredict]{predictGlycansParam-class}
#'
#' @export
#' @examples
#' pgp <- predictGlycansParam()
#' pgp@@dp <- c(1,7)
#' pgp@@ESI_mode <- 'neg'
#' pgp@@scan_range <- c(150, 1300)
#' pgp@@modifications <- c('sulphate', 'carboxyl')
#' pgp@@double_sulphate <- TRUE
#' predicted.df <- predictGlycans(param = pgp)
#' 
#' @details 
#' This function is intended for “calculation” of all possible glycans (or sugars) 
#' within a set of constraining parameters (contained within the 
#' \code{predictGlycansParam} object). Specifically, the user indicates which 
#' monomer types (hexose only or hexose and pentose), degree of polymerisation (length) 
#' range and modification types should be included, the desired maximum 
#' for the average number of modifications per monomer and whether 
#' mono-/oligosaccharides are procainamide-labelled or not. There is also the 
#' option for whether or not double sulphation of a single monomer is possible or not. 
#' The function “builds” names, formulas and masses for all sugars possible 
#' within the constraining parameters. The user also provides two parameters 
#' related to mass spectrometry: ionisation mode (\code{ESI_mode}) and 
#' scan range (\code{scan_range}). m/z values of ions are calculated 
#' depending on the ionisation mode and modifications. Sugars which contain 
#' no ions with m/z values within the given scan range are removed. The final 
#' output is returned as a (wide format) dataframe. This package was written 
#' for annotation of mass spec data (especially LC-MS) but if used for 
#' other purposes either ionisation mode can be given and very wide scan ranges. 
#' The function works by sourcing a python file and then using the function 
#' encoded in the python script.
#' 
#' **Word of caution: please note that this tool will predict some sugars that 
#' are not really ‘possible’ as the nature of sugar chemistry means that it 
#' would take a long time to add in all the constraints!**
#' 
#' For more details see the vignette:
#' \code{vignette("glycanPredict", package = "glycanPredict")}
#' 
#' @seealso 
#' glycanPredict::predictGlycansParam()
#' 

predictGlycans <- function(param){
  path <- paste(system.file(package="glycanPredict"), "sugarMassesPredict-r_package-maldi_20230317.py", sep="/")
  #check if pandas installed
  if(!reticulate::py_module_available("pandas")){
    reticulate::py_install("pandas")
  }
  if(!reticulate::py_module_available("numpy")){
    reticulate::py_install("numpy")
  }
  #source python script
  reticulate::source_python(path)
  #get parameters from object
  #and reset class types for python
  dp = as.list(as.integer(param@dp))
  ESI_mode = param@ESI_mode
  scan_range = as.list(as.integer(param@scan_range))
  pent_option = param@pent_option
  nmod_max = as.integer(param@nmod_max)
  modifications = as.list(param@modifications)
  double_sulphate = param@double_sulphate
  label = param@label
  df <- predict_sugars(dp = dp, ESI_mode = ESI_mode,
                       scan_range = scan_range,
                       pent_option = pent_option, modifications = modifications,
                       label = label, nmod_max = nmod_max, 
                       double_sulphate = double_sulphate)
  return(df)
}
