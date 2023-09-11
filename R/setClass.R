#' predictGlycansParam class
#' 
#' @description 
#' An S4 class that contains the parameters for glycan prediction.
#' Once created the predictGlycansParam object should be supplied to the function
#' \link[glycanPredict]{predictGlycans} for glycan prediction constrained by the parameters contained
#' in the \code{predictGlycansParam} object.
#' 
#' @seealso glycanPredict::predictGlycans()
#' 
#' @export predictGlycansParam
#' @exportClass predictGlycansParam
#' 
#' @usage 
#' pgp <- predictGlycansParam()
#' df <- predictGlycans(param = pgp)
#' 
#' @slot dp Degree of polymerisation range (numeric, length 2).
#' @slot polarity ionisation mode used. Accepts 'pos' or 'neg'.
#' @slot scan_range Scan range used during MS. (numeric, length 2).
#' @slot pent_option Logical: Should pentose monomers be included? 
#' @slot modifications Modifications to be considered. Any combination of 'carboxyl', 'phosphate', 'deoxy', 'nacetyl', 'omethyl', 'anhydrobridge', 'oacetyl', 'unsaturated', 'alditol', 'amino', 'dehydrated', 'sulphate' or 'all' or 'none' (default)
#' @slot nmod_max Maximum number of modifications per monomer on average (default 1). Does not take into account unsaturated, alditol or dehydrated.
#' @slot double_sulphate Logical: can monomers be double-sulphated. If \code{TRUE} you MUST give a value of at least 2 to nmod_max.
#' @slot label Are sugars labelled? Currently only accepts 'none' or 'procainamide'.
#' @slot ion_type Ionisation type. Currently accepted ESI and MALDI. Impacts ions.
#' 
#' @inherit predictGlycans details
#' 
#' @inherit predictGlycans examples
#'

predictGlycansParam = setClass("predictGlycansParam",
         slots = c(
           dp = "numeric",
           polarity = "character",
           scan_range = "numeric",
           pent_option = "logical",
           modifications = "character",
           nmod_max = "numeric",
           double_sulphate = "logical",
           label = "character",
           ion_type = "character",
           format = "character",
           adducts = "character"
         ),
         prototype = prototype(
           dp = c(1, 6),
           polarity = "neg",
           scan_range = c(175, 1400),
           pent_option = FALSE,
           modifications = "none",
           nmod_max = 1,
           double_sulphate = FALSE,
           label = "none",
           ion_type = "ESI",
           format = "long",
           adducts = c("H", "Cl", "nH")
         ),
         validity = function(object) {
           msg <- character()
           if (length(object@dp) != 2 | any(object@dp < 0))
             msg <- c(msg, paste0("'dp' has to be a numeric",
                                  " of length 2 with only positive",
                                  " values."))
           if(!any(object@polarity %in% c("neg", "pos")))
             msg <- c(msg, paste0("'polarity' has to be a character",
                                  " containing 'pos' and/or 'neg"))
           if (length(object@scan_range) != 2 | any(object@scan_range < 0))
             msg <- c(msg, paste0("'scan_range' has to be a numeric",
                                  " of length 2 with only positive",
                                  " values."))
           if (length(object@pent_option) != 1 | 
               class(object@pent_option) != "logical")
             msg <- c(msg, paste0("'pent_option' has to be a logical of",
                                  " length 1."))
           possible_modifications <-  c('none', 'all', 'carboxyl', 'phosphate', 
                                        'deoxy', 'nacetyl', 'omethyl',
                                        'anhydrobridge', 'oacetyl', 'unsaturated', 
                                        'alditol', 'amino','dehydrated','sulphate')
           if (!all(object@modifications %in% possible_modifications))
             msg <- c(msg, paste0("valid options for 'modifications' are: ",
                                  paste0("'", possible_modifications, "'",
                                         collapse = ", "), "."))
           if (length(object@nmod_max) != 1 | any(object@nmod_max <= 0))
             msg <- c(msg, paste0("'nmod_max' has to be positive numeric",
                                  " of length 1."))
           if (length(object@double_sulphate) != 1 | 
               class(object@double_sulphate) != "logical")
             msg <- c(msg, paste0("'double_sulphate' has to be a logical of",
                                  " length 1."))
           possible_labels <- c("none", "procainamide", "proca","procA", "ProA",
                                "2-ap","2-AP","pa","PA","2-aminopyridine",
                                "2-aa", "2-AA","aba", "ABA","2-aminobenzoic acid",
                                "2-ab", "2-AB", "ab", "AB", "2-aminobenzamide",
                                "pmp", "PMP", "1-phenyl-3-methyl-5-pyrazolone")
           if(length(object@label) != 1 | !any(object@label %in% possible_labels))
             msg <- c(msg, paste0("'label' has to be a character",
                                  " of length 1",
                                  " containing ",
                                  possible_labels))
           if(length(object@format) != 1 | !any(object@format %in% c("wide", "long")))
             msg <- c(msg, paste0("'format' has to be a character",
                                  " of length 1. the only allowed options",
                                  " are 'wide' and 'long'"))
           if(!object@ion_type %in% c("ESI", "MALDI") | length(object@ion_type) != 1
             msg <- c(msg, paste0("'ion_type' has to be a character",
                                  " of length 1. the options are 'ESI' or 'MALDI'"))
             }
           possible_adducts <- c("H", "Na", "NH4", "K", "Cl", "CHOO", "nH")
           if (!all(object@adducts %in% possible_adducts))
             msg <- c(msg, paste0("valid options for 'adducts' are: ",
                                  paste0("'", possible_adducts, "'",
                                         collapse = ", "), "."))
           if (length(msg) >= 1)
             print(msg)
           else
             TRUE
         })
