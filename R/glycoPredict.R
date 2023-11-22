#' Function to predict glycan masses and m/z values
#'
#' @include setClass.R
#'
#' @description This function will predict all possible glycan molecules within the constraints of the parameters.
#' @param param A \code{glycoPredictParam} object containing all parameters for the
#'     glycoPredict function. See \link[GlycoAnnotateR]{glycoPredictParam-class}
#'
#' @export
#' @examples
#' pgp <- glycoPredictParam()
#' pgp@@dp <- c(1,7)
#' pgp@@polarity <- 'neg'
#' pgp@@scan_range <- c(150, 1300)
#' pgp@@modifications <- c('sulphate', 'carboxylicacid')
#' pgp@@double_sulphate <- TRUE
#' predicted.df <- glycoPredict(param = pgp)
#' 
#' @details 
#' This function is intended for “calculation” of all possible glycans (or sugars) 
#' within a set of constraining parameters (contained within the 
#' \code{glycoPredictParam} object). Specifically, the user indicates which 
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
#' 
#' @seealso 
#' glycoAnnotateR::glycoPredictParam()
#' 

glycoPredict <- function(param){
  path <- paste(system.file(package="GlycoAnnotateR"), "sugarMassesPredict.py", sep="/")
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
  polarity = param@polarity
  scan_range = as.list(as.integer(param@scan_range))
  pent_option = param@pent_option
  nmod_max = as.integer(param@nmod_max)
  modifications = as.list(param@modifications)
  double_sulphate = param@double_sulphate
  label = param@label
  ion_type = param@ion_type
  adducts = as.list(param@adducts)
  naming = as.list(param@naming)
  glycan_linkage = as.list(param@glycan_linkage)
  df <- predict_sugars(dp = dp, polarity = polarity,
                       scan_range = scan_range,
                       pent_option = pent_option, modifications = modifications,
                       label = label, nmod_max = nmod_max, ion_type = ion_type,
                       double_sulphate = double_sulphate, adducts = adducts,
                       naming = naming, glycan_linkage = glycan_linkage)
  format = param@format
  library(magrittr)
  if(format == "long"){
    as.num = function(x, na.strings = "NA") {
      stopifnot(is.character(x))
      na = x %in% na.strings
      x[na] = "0"
      x = as.numeric(x)
      x[na] = NA_real_
      x
    }
    df.l <- df %>% 
      #make long
      tidyr::pivot_longer(cols = starts_with("[M"),
                          names_to = "ion",
                          values_to = "mz") %>% 
      #remove ions outside scan range
      tidyr::drop_na(mz) %>% 
      #calculate ion formula
      dplyr::mutate(C = stringr::str_split_i(formula, "C", 2) %>% 
                      sub("\\D.*", "", .) %>% 
                      as.num(),
                    H = stringr::str_split_i(formula, "H", 2) %>% 
                      sub("\\D.*", "", .) %>% 
                      as.num(),
                    N = stringr::str_split_i(formula, "N", 2) %>% 
                      sub("\\D.*", "", .) %>% 
                      as.num(),
                    N = dplyr::case_when(grepl("N", formula) & is.na(N) ~ 1,
                                         TRUE ~ N),
                    O = stringr::str_split_i(formula, "O", 2) %>% 
                      sub("\\D.*", "", .) %>% 
                      as.num(),
                    P = stringr::str_split_i(formula, "P", 2) %>% 
                      sub("\\D.*", "", .) %>% 
                      as.num(),
                    P = dplyr::case_when(grepl("P", formula) & is.na(P) ~ 1,
                                         TRUE ~ P),
                    S = stringr::str_split_i(formula, "S", 2) %>% 
                      sub("\\D.*", "", .) %>% 
                      as.num(),
                    S = dplyr::case_when(grepl("S", formula) & is.na(S) ~ 1,
                                         TRUE ~ S),
                    ion_effect = gsub("\\[M|\\].*", "", ion),
                    delta_H = sub(".*([+-]\\d*H).*", "\\1", ion_effect) %>% 
                      sub("[-+]\\d[^H].*|[-+][A-G, I-Z].*", "", .) %>% 
                      sub("H", "", .) %>% 
                      sub("^-$", -1, .) %>% 
                      sub("^\\+$", 1, .) %>% 
                      as.num())
    df.l$delta_H[df.l$ion_effect == "+NH4"] <- 4
    df.l <- df.l %>% 
      dplyr::mutate(delta_N = sub(".*([+-]\\d*N[^a]).*", "\\1", ion_effect) %>% 
                      sub("[+-]Na", "", .) %>% 
                      sub("[-+]\\d[^N].*|[-+][A-M, O-Z].*|[A-M, O-Z]", "", .) %>% 
                      sub("N", "", .) %>% 
                      sub("^-$", -1, .) %>% 
                      sub("^\\+$", 1, .) %>% 
                      as.num(),
                    delta_Cl = sub(".*([+-]\\d*Cl).*", "\\1", ion_effect) %>%  
                      sub("[-+]\\d[^Cl].*|[-+][A-B, D-Z].*", "", .) %>% 
                      sub("Cl", "", .) %>% 
                      sub("^-$", "-1", .) %>% 
                      sub("^\\+$", "1", .) %>% 
                      as.num(na.strings = "+CHOO"),
                    delta_Na = sub(".*([+-]\\d*Na).*", "\\1", ion_effect) %>% 
                      sub("[-+]\\d[^Na].*|[-+][A-M, O-Z].*", "", .) %>% 
                      sub("Na", "", .) %>% 
                      sub("^-$", -1, .) %>% 
                      sub("^\\+$", 1, .) %>% 
                      as.num(na.strings = "+NH4"),
                    delta_K = sub(".*([+-]\\d*K).*", "\\1", ion_effect) %>% 
                      sub("[-+]\\d[^K].*|[-+][A-J, L-Z].*", "", .) %>% 
                      sub("K", "", .) %>% 
                      sub("^-$", -1, .) %>% 
                      sub("^\\+$", 1, .) %>% 
                      as.num())
    df.l[is.na(df.l)] <- 0
    df.l <- df.l %>% 
      dplyr::mutate(ion_formula = paste0("C", C, 
                                         "Cl", delta_Cl,
                                         "H", H + delta_H,
                                         "K", delta_K,
                                         "N", N + delta_N,
                                         "Na", delta_Na,
                                         "O", O,
                                         "S", S, "P", P) %>% 
                      gsub("[A-Z]0|Na0|Cl0", "", .) %>% 
                      gsub("(\\b|\\D)1(\\b|\\D)", "\\1\\2", .))
    df <- df.l %>% 
      dplyr::select(!matches("delta_|^[[:upper:]][a,c]?$|_effect"))
    
  }
  return(df)
}