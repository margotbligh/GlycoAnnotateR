#' Predict masses and m/z values of theoretical glycans
#'
#' @include setClass.R
#'
#' @description \code{glycoPredict()} predicts all possible glycan within the
#' constraints set by the \code{glycoPredictParam} object.
#' @param param A \code{glycoPredictParam} object. See \link[GlycoAnnotateR]{glycoPredictParam}
#'
#' @export
#'
#' @examples
#' gpp <- glycoPredictParam()
#' gpp@@dp <- c(1,7)
#' gpp@@polarity <- 'neg'
#' gpp@@scan_range <- c(150, 1300)
#' gpp@@modifications <- c('sulfate', 'carboxylicacid')
#' predicted.df <- glycoPredict(param = gpp)
#'
#' @details
#' \code{glycoPredict()} is used to predict masses and mass to charge ratios of all theoretically
#' possible glycans within a set of constraining parameters (defined in the
#' \code{glycoPredictParam} object). This package was written
#' for annotation of mass spec data (especially LC-MS) but if used for
#' other purposes either ionisation mode and very wide scan ranges can be given.
#' The function works by sourcing a python file and then using the function
#' encoded in the python script.
#'
#' @seealso
#' \link[GlycoAnnotateR]{glycoPredictParam}
#'

glycoPredict <- function(param){
  path <- paste(system.file(package="GlycoAnnotateR"), "sugarMassesPredict.py", sep="/")
  #check if pandas installed
  if(!reticulate::py_suppress_warnings(reticulate::py_module_available("pandas"))){
    reticulate::py_install("pandas")
  }
  if(!reticulate::py_suppress_warnings(reticulate::py_module_available("numpy"))){
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
  double_sulfate = param@double_sulfate
  label = param@label
  ion_type = param@ion_type
  adducts = as.list(param@adducts)
  naming = as.list(param@naming)
  glycan_linkage = as.list(param@glycan_linkage)
  modification_limits = param@modification_limits

  #message(paste("Glycans will be predicted according to the following glycoPredictParam() object:\n", str(param)))

  df <- predict_sugars(dp = dp, polarity = polarity,
                       scan_range = scan_range,
                       pent_option = pent_option, modifications = modifications,
                       label = label, nmod_max = nmod_max, ion_type = ion_type,
                       double_sulfate = double_sulfate, adducts = adducts,
                       naming = naming, glycan_linkage = glycan_linkage,
                       modification_limits = modification_limits)
  format = param@format
  library(magrittr)
  library(tidyr)
  library(dplyr)
  library(data.table)
  library(magrittr)
  library(stringr)

  #read in formatted composition db
  db <- fread(paste0(system.file(package="GlycoAnnotateR"),
                     '/db/compositions_formatted.tsv'))

  #get shared modifications
  modifications_shared <- names(db)[names(db) %in% param@modifications]
  if(param@pent_option == TRUE){
    modifications_shared <- c(modifications_shared, 'pent')
  }

  #get names to drop after matching
  db_names_to_drop <- names(db)[!grepl('_id', names(db))]
  db_names_to_drop <- c(db_names_to_drop,
                        param@modifications) %>%
    unique()

  #match predictions with db
  df <- df %>%
    dplyr::mutate(mass_rounded = round(mass, 2)) %>%
    dplyr::left_join(db %>%
                       dplyr::mutate(mass_rounded = round(mass, 2)),
                     by = c('mass_rounded', 'hex', modifications_shared)) %>%
    dplyr::select(!any_of(db_names_to_drop)) %>%
    dplyr::select(!any_of(c('mass_rounded', 'mass.y'))) %>%
    dplyr::rename(mass = `mass.x`)

  #change NAs to nones
  df$glyconnect_id[is.na(df$glyconnect_id)] <- 'none'
  df$glytoucan_id[is.na(df$glytoucan_id)] <- 'none'

  if(format == "long"){
    as.num = function(x, na.strings = "NA") {
      stopifnot(is.character(x))
      na = x %in% na.strings
      x[na] = "0"
      x = as.numeric(x)
      x[na] = NA_real_
      x
    }
    suppressWarnings(df.l <- df %>%
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
                                       as.num(),
                                     multiple_ammonium = dplyr::case_when(grepl('\\dNH4', ion) ~
                                                                            stringr::str_split_i(ion_effect,
                                                                                                 '\\+|N', 2) %>%
                                                                            as.numeric(),
                                                                          TRUE ~ NA),
                                     delta_H = dplyr::case_when(grepl('\\dNH4', ion) ~
                                                                  delta_H + (multiple_ammonium*4),
                                                                grepl('\\+NH4', ion) ~ 4,
                                                                grepl('\\+CHOO', ion) ~ 1,
                                                                TRUE ~ delta_H),
                                     delta_N = sub(".*([+-]\\d*N[^a]).*", "\\1", ion_effect) %>%
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
                                       as.num(),
                                     delta_C = dplyr::case_when(grepl('\\+CHOO', ion) ~ 1,
                                                                TRUE ~ 0),
                                     delta_O = dplyr::case_when(grepl('\\+CHOO', ion) ~ 2,
                                                                TRUE ~ 0)))
    df.l[is.na(df.l)] <- 0
    df.l <- df.l %>%
      dplyr::mutate(ion_formula = paste0("C", C + delta_C,
                                         "Cl", delta_Cl,
                                         "H", H + delta_H,
                                         "K", delta_K,
                                         "N", N + delta_N,
                                         "Na", delta_Na,
                                         "O", O + delta_O,
                                         "S", S, "P", P) %>%
                      gsub("[A-Z]0|Na0|Cl0", "", .) %>%
                      gsub("(\\b|\\D)1(\\b|\\D)", "\\1\\2", .))
    df <- df.l %>%
      dplyr::select(!matches("delta_|^[[:upper:]][a,c]?$|_effect|multiple_")) %>%
      dplyr::mutate(charge = stringr::str_split_i(ion, '\\]', 2) %>%
                      sub('\\+$', '+1', .) %>%
                      sub('\\-$', '-1', .))

  }

  if (nrow(df) == 0){
    warning('Output has zero rows! Check your scan range, adducts/polarity and DP range are sensible')
  }
  return(df)
}
