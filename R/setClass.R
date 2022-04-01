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
#' @usage 
#' pgp <- predictGlycansParam()
#' df <- predictGlycans(param = pgp)
#' 
#' @slot dp Degree of polymerisation range (numeric, length 2).
#' @slot ESI_mode ESI mode used. Accepts 'pos' or 'neg'.
#' @slot scan_range Scan range used during MS. (numeric, length 2).
#' @slot pent_option Logical: Should pentose monomers be included? 
#' @slot modifications Modifications to be considered. Any combination of 'carboxyl', 'phosphate', 'deoxy', 'nacetyl', 'omethyl', 'anhydrobridge', 'oacetyl', 'unsaturated', 'alditol', 'amino', 'dehydrated', 'sulphate' or 'all' or 'none' (default)
#' @slot nmod_max Maximum number of modifications per monomer on average (default 1). Does not take into account unsaturated, alditol or dehydrated.
#' @slot double_sulphate Logical: can monomers be double-sulphated. If \code{TRUE} you MUST give a value of at least 2 to nmod_max.
#' @slot label Are sugars labelled? Currently only accepts 'none' or 'procainamide'.
#' 
#' @inherit predictGlycans details
#' 
#' @inherit predictGlycans examples
#' 
#' @export

predictGlycansParam = setClass("predictGlycansParam",
         slots = c(
           dp = "numeric",
           ESI_mode = "character",
           scan_range = "numeric",
           pent_option = "logical",
           modifications = "character",
           nmod_max = "numeric",
           double_sulphate = "logical",
           label = "character"
         ),
         prototype = prototype(
           dp = c(1, 6),
           ESI_mode = "neg",
           scan_range = c(175, 1400),
           pent_option = FALSE,
           modifications = "none",
           nmod_max = 1,
           double_sulphate = FALSE,
           label = "none"
         ),
         validity = function(object) {
           msg <- character()
           if (length(object@dp) != 2 | any(object@dp < 0))
             msg <- c(msg, paste0("'dp' has to be a numeric",
                                  " of length 2 with only positive",
                                  " values."))
           if(!any(object@ESI_mode %in% c("neg", "pos")))
             msg <- c(msg, paste0("'ESI_mode' has to be a character",
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
           if (!(object@modifications) %in% possible_modifications)
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
           if(length(object@label) != 1 | !any(object@label %in% c("none", 
                                                                   "procainamide")))
             msg <- c(msg, paste0("'label' has to be a character",
                                  " of length 1",
                                  " containing 'none' or 'procainamide"))
           if (length(msg) >= 1)
             print(msg)
           else
             TRUE
         })
