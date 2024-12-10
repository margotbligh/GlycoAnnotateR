#' Annotate m/z values of fragements in MS2 spectra from glycans
#'
#' @include setClass.R
#' @include glycoPredict.R
#' @include glycoAnnotate.R
#'
#' @description \code{glycoMS2Annotate()} annotates fragments in MS2 spectra
#' from precursors annotated as glycans by \link[GlycoAnnotateR]{glycoAnnotate}.
#'
#' @slot precursorAnnotations Vector of precursor annotations including IUPAC name,
#' ion and DP, separated by a colon (:) unless otherwise specified. Assumes
#' multiple annotations (if multiple annotations assigned to one precursor) are
#' separated by a comma!
#' @slot precursorAnnotations_sep Non-default separator of IUPAC name, ion and DP
#' in `precursorAnnotations`
#' @slot ms2spectra MS2 spectra from fragmentation of annotated precursors. MUST be
#' a dataframe with an 'mz' column for fragment ions and a 'precursorAnnotations'
#' column with precursor annotations (matching those in `precursorAnnotations`)
#' @slot ion_type Ionisation type. Currently accepted ESI and MALDI. Impacts ions.
#' @slot error Numeric value - error used to create window for matching. mz values
#' will be matched against theoretical mzs +/- error.
#' @slot error_units Units for error - can be 'ppm' or 'Da'
#' @slot nmod_max  Maximum number of modifications per monomer, calculated by the number of modifications over the number of monomers (default 1).
#' Does not take into account unsaturated, alditol or dehydrated.
#' @slot double_sulfate Logical. Can monomers be double-sulfated. If \code{TRUE}, nmod_max needs to have a value of at least 2.
#' @slot label Are sugars labelled by reductive amination? Current supported labels are: "none", "procainamide","2-aminobenzoic acid",
#' "2-aminobenzamide", "1-phenyl-3-methyl-5-pyrazolone".
#' @slot dehydrations Logical. If TRUE 'dehydrated' will be included in modifications
#' to look for water losses.
#'
#' @export
#'
#' @seealso \link[GlycoAnnotateR]{glycoPredictParam}
#' @seealso \link[GlycoAnnotateR]{glycoPredict}
#' @seealso \link[GlycoAnnotateR]{glycoAnnotate}
#' @seealso \link[GlycoAnnotateR]{glycoMS2Extract}



glycoMS2Annotate <- function(precursorAnnotations,
                             precursorAnnotations_sep = NULL,
                             ms2spectra,
                             error = 3,
                             error_units = 'ppm',
                             ion_type = 'ESI',
                             polarity = 'neg',
                             nmod_max = 1,
                             double_sulfate = FALSE,
                             label = 'none',
                             dehydrations = FALSE){
  #check input validity
  if(!is.vector(precursorAnnotations)){
    stop('precursorAnnotations is not a vector')
  }
  if(!all(grepl(':', precursorAnnotations)) &
     is.null(precursorAnnotations_sep)){
    stop('No : separators in precursorAnnotations and no alternative provided')
  }
  if(!is.null(precursorAnnotations_sep)){
    if(!all(grepl(precursorAnnotations_sep, precursorAnnotations))){
      stop('No precursorAnnotations_sep in precursorAnnotations')
    }
  }
  if(class(ms2spectra) != 'data.frame'){
    stop('ms2spectra is not a data.frame')
  }
  if(error_units != 'ppm' & error_units != 'Da'){
    stop('Error units are not ppm or Da')
  }
  if(!is.numeric(error)){
    stop('Error is not numeric')
  }
  if(!all(polarity %in% c('pos', 'neg'))){
    stop('Polarity is not pos and/or neg')
  }
  if(ion_type != 'ESI' & ion_type != 'MALDI'){
    stop('Ion_type is not ESI or MALDI')
  }
  if(!'mz' %in% names(ms2spectra)){
    stop('no mz column in ms2spectra')
  }
  if(!'precursorAnnotations' %in% names(ms2spectra)){
    stop('no precursorAnnotations column in ms2spectra')
  }
  if(!any(precursorAnnotations %in% ms2spectra$precursorAnnotations)){
    stop('none of the precursorAnnotations are in the precursorAnnotations',
         'column of the ms2spectra dataframe')
  }

  #use default sep if none provided
  if(is.null(precursorAnnotations_sep)){ precursorAnnotations_sep = ':'}

  #make an empty list
  ms2_annot_list <- list()

  #fill in list with dfs for annotated spectra, with one entry per annotation
  for(i in 1:length(precursorAnnotations)){
    #get annotation
    annot = precursorAnnotations[i]

    #get dp from annotation
    if (grepl(',', annot)){
      #if multiple annotations, get maximum dp
      annots <- stringr::str_split_1(annot, ',')
      dps <- stringr::str_split_i(annots,
                                  precursorAnnotations_sep,
                                  3) %>% as.numeric()
      max_dp = max(dps)
    } else {
      #otherwise just get the dp
      max_dp = stringr::str_split_i(annot,
                                    precursorAnnotations_sep,
                                    3) %>% as.numeric()
    }

    #get modifications from annotations
    #split by separator or space
    p = paste0(' |', precursorAnnotations_sep)
    annot_split <- stringr::str_split_1(annot, pattern = p)
    #remove elements with the ion
    annot_split <- annot_split[!grepl('\\[', annot_split)] %>%
      #remove numbers
      gsub('\\d', '', .) %>%
      #make unique
      unique()
    #remove hexose (always included)
    #only keep elements with word characters
    annot_split <- annot_split[annot_split != 'Hex' &
                                 grepl('\\w',annot_split)]
    #if no modifications, set to 'none'
    if(length(annot_split) == 0){
      annot_split <- 'none'
    }

    #make sure everything is lower case
    annot_split_lc <- stringr::str_to_lower(annot_split) %>%
      sub('-', '', .)

    #get pentose option
    if(any(annot_split_lc == 'pen')){
      pent_option = T
      annot_split_lc <- annot_split_lc[!annot_split_lc == 'pen']
    } else {
      pent_option = F
    }
    
    #change deoxyhex to deoxy
    if(any(annot_split_lc == 'deoxyhex')){
      annot_split_lc[annot_split_lc == 'deoxyhex'] <- 'deoxy'
    }
    
    if(dehydrations == TRUE){
      annot_split_lc <- c(annot_split_lc, 'dehydrated')
    }

    #create parameter object
    gpp <- GlycoAnnotateR::glycoPredictParam(dp = c(1, max_dp),
                                             pent_option = pent_option,
                                             modifications = annot_split_lc,
                                             polarity = polarity,
                                             scan_range = c(0, 10000), #set the scan range wide
                                             ion_type = ion_type,
                                             nmod_max = nmod_max,
                                             double_sulfate = double_sulfate,
                                             label = label)

    #make df for annotation
    df <- ms2spectra %>%
      #filter for precursor annotation
      dplyr::filter(precursorAnnotations == !!annot)

    #annotate fragments
    df_annot <- GlycoAnnotateR::glycoAnnotate(df,
                                              param = gpp,
                                              collapse = T,
                                              error = error,
                                              error_units = error_units)

    #add to list
    ms2_annot_list[[i]] <- df_annot

  }

  #bind dfs in list into one df
  ms2spectra_annotated <- ms2_annot_list %>%
    dplyr::bind_rows()

  #return output
  return(ms2spectra_annotated)


}



