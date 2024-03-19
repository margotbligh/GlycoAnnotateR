#' glycoPredictParam class
#'
#' @description
#' An S4 class that contains the parameters for \link[GlycoAnnotateR]{glycoPredict()} function.
#'
#' @seealso GlycoAnnotateR::glycoPredict()
#'
#' @export glycoPredictParam
#' @exportClass glycoPredictParam
#'
#' @slot dp Degree of polymerisation range. It must be a numeric vector of length 2. Default: c(1,6).
#' @slot polarity Ionisation mode used. Accepts 'pos' or 'neg' (default).
#' @slot scan_range Scan range used during MS. It must be a numeric vector of length 2. Default: c(175, 1400).
#' @slot pent_option Logical. Should pentose monomers be included? Default: \code{FALSE}.
#' @slot modifications Modifications to be considered. Any combination of 'carboxylicacid', 'sialicacid',
#' 'phosphate', 'deoxy', 'nacetyl', 'omethyl', 'anhydrobridge', 'oacetyl',
#' 'unsaturated', 'alditol', 'amino', 'dehydrated', 'sulfate', 'aminopentyllinker' or 'all' or 'none' (default)
#' @slot nmod_max Maximum number of modifications per monomer, calculated by the number of modifications over the number of monomers (default 1).
#' Does not take into account unsaturated, alditol or dehydrated.
#' @slot double_sulfate Logical. Can monomers be double-sulfated. If \code{TRUE}, nmod_max needs to have a value of at least 2.
#' @slot label Are sugars labelled by reductive amination? Current supported labels are: "none", "procainamide","2-aminobenzoic acid",
#' "2-aminobenzamide", "1-phenyl-3-methyl-5-pyrazolone".
#' @slot ion_type Ionisation type. Currently accepted ESI and MALDI. Impacts ions.
#' @slot naming Notation for molecule names. Uses commonly accepted abbreviations. Possibilities: 'IUPAC' (default), 'Oxford', 'GlycoCT'
#' @slot adducts Adduct types to be included. Options are 'H', 'Na', 'K', 'NH4', 'Cl' 
#' and 'CHOO'. See [here](https://margotbligh.github.io/GlycoAnnotateR/#output-and-other-parameters) 
#' for detailed description of which adducts are generated,
#' @slot glycan_linkage Option to implement filters for O- and N-glycans. Possibilities: 'none' (default), 'nglycan' or 'oglycan'.
#' @slot modification_limits Option to implement user created filters. Must be a named list, with names as modifications and values as limits.
#' @slot format Output format. Options are 'long' (default) or wide.
#' 
#' @inherit glycoPredict details
#' @inherit glycoPredict examples

glycoPredictParam = setClass("glycoPredictParam",
         slots = c(
           dp = "numeric",
           polarity = "character",
           scan_range = "numeric",
           pent_option = "logical",
           modifications = "character",
           nmod_max = "numeric",
           double_sulfate = "logical",
           label = "character",
           ion_type = "character",
           format = "character",
           adducts = "character",
           naming = "character",
           glycan_linkage = "character",
           modification_limits = "ANY"
         ),
         prototype = prototype(
           dp = c(1, 6),
           polarity = "neg",
           scan_range = c(175, 1400),
           pent_option = FALSE,
           modifications = "none",
           nmod_max = 1,
           double_sulfate = FALSE,
           label = "none",
           ion_type = "ESI",
           format = "long",
           adducts = "all",
           naming = "IUPAC",
           glycan_linkage = "none",
           modification_limits = "none"
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
           possible_modifications <-  c('none', 'all', 'carboxylicacid',
                                        'sialicacid', 'phosphate',
                                        'deoxy', 'nacetyl', 'omethyl',
                                        'anhydrobridge', 'oacetyl', 'unsaturated',
                                        'alditol', 'amino','dehydrated','sulfate', 'aminopentyllinker')
           if (!all(object@modifications %in% possible_modifications))
             msg <- c(msg, paste0("valid options for 'modifications' are: ",
                                  paste0("'", possible_modifications, "'",
                                         collapse = ", "), "."))
           if(object@modification_limits != "none")
             if (!all(names(object@modification_limits) %in% possible_modifications))
               msg <- c(msg, paste0("valid options for 'modifications' are: ",
                                    paste0("'", possible_modifications, "'",
                                           collapse = ", "), ".",
                                    "modification limits must be a list with",
                                    "names as modifications"))
           if(object@modification_limits != "none")
             if (!all(names(object@modification_limits) %in% object@modifications))
               msg <- c(msg, paste0('names in modification limits do not match modifications!'))
           if (length(object@nmod_max) != 1 | any(object@nmod_max <= 0) |
               any(object@nmod_max > 3))
             msg <- c(msg, paste0("'nmod_max' has to be numeric",
                                  " of length 1 between 1 and 3."))
           if (length(object@double_sulfate) != 1 |
               class(object@double_sulfate) != "logical")
             msg <- c(msg, paste0("'double_sulfate' has to be a logical of",
                                  " length 1."))
           possible_labels <- c("none", "procainamide", "proca","procA", "ProA",
                                "2-ap","2-AP","pa","PA","2-aminopyridine",
                                "2-aa", "2-AA","aba", "ABA","2-aminobenzoic acid",
                                "2-ab", "2-AB", "ab", "AB", "2-aminobenzamide",
                                "pmp", "PMP", "1-phenyl-3-methyl-5-pyrazolone")
           if (length(object@label) != 1 | !any(object@label %in% possible_labels))
             msg <- c(msg, paste0("'label' has to be a character",
                                  " of length 1",
                                  " containing ",
                                  possible_labels))
           if (length(object@format) != 1 | !any(object@format %in% c("wide", "long")))
             msg <- c(msg, paste0("'format' has to be a character",
                                  " of length 1. the only allowed options",
                                  " are 'wide' and 'long'"))
           if (!object@ion_type %in% c("ESI", "MALDI") | length(object@ion_type) != 1)
             msg <- c(msg, paste0("'ion_type' has to be a character",
                                  " of length 1. the options are 'ESI' or 'MALDI'"))

           possible_adducts <- c("all", "H", "Na", "NH4", "K", "Cl", "CHOO", "nH")
           if (!all(object@adducts %in% possible_adducts))
             msg <- c(msg, paste0("valid options for 'adducts' are: ",
                                  paste0("'", possible_adducts, "'",
                                         collapse = ", "), "."))
           possible_namings <- c("IUPAC", "GlycoCT", "Oxford")
           if (!all(object@naming %in% possible_namings))
             msg <- c(msg, paste0("valid options for 'naming' are: ",
                                  paste0("'", possible_namings, "'",
                                         collapse = ", "), "."))
           if (!object@glycan_linkage %in% c("nglycan", "oglycan", "none"))
             msg <- c(msg, paste0("valid options for 'glycan_linkage' are: ",
                                  "none, nglycan or oglycan"))
           if (length(msg) >= 1)
             print(msg)
           else
             TRUE
         })
