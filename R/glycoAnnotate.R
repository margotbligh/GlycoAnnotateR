#' Annotate m/z values with glycan compositions based on theoretical values
#'
#' @include setClass.R
#' @include glycoPredict.R
#'
#' @description \code{glycoAnnotate()} annotates peaks or features in MS data,
#' using EITHER a pre-generated table (`pred_table`) by 
#' \link[GlycoAnnotateR]{glycoPredict} OR
#' by generating a new table by suppling a `glycoPredictParam`object
#' to `param`.
#'
#' @export
#'
#' @slot data Dataframe containing data to be annotated. For example,
#' feature dataframe from XCMS pre-processing (LC-MS or direct inject)
#' or features from Cardinal (MALDI).
#' @slot mz_column Name of column containing m/z values.
#' @slot mzmin_column OPTIONAL: Name of column containing minimum m/z data values.
#' If supplied, will do overlap of two ranges matching. Generally only if mzmin and mzmax
#' values generated during peak picking. If not provided, mz value will be annotated
#' if within range of theoretical mz +/- error (value within range).
#' @slot mzmax_column OPTIONAL: Name of column containing maximum m/z data values.
#' If supplied, will do overlap of two ranges matching. Generally only if mzmin and mzmax
#' values generated during peak picking.If not provided, mz value will be annotated
#' if within range of theoretical mz +/- error (value within range).
#' @slot pred_table Table generated previously by \link[GlycoAnnotateR]{glycoPredict}.
#' MUST BE LONG FORMAT - select \code{format='long'} when running prediction.
#' Must provide value for this OR `param`.
#' @slot param \link[GlycoAnnotateR]{glycoPredictParam} object for generation of table
#' of theoretical mz values for annotation. Must provide value for this OR `pred_table`.
#' @slot collapse Logical. If \code{TRUE}, annotations will be 'collapsed' so that multiple
#' annotations for one mz will be in the same row, comma separated (nrow of output is in
#' this case equal to nrow of input data). If \code{FALSE} (default), it is possible
#' that rows in the input dataframe are repeated with different annotations. The
#' information on annotations is more detailed in this case. Collapsing can also be done
#' afterwards on the output using \link[GlycoAnnotateR]{glycoAnnotationsCollapse}.
#' @slot collapse_columns Columns to be pasted together before collapsing.
#' Only needed if \code{collapse=TRUE} and non-default columns wanted - default is
#' molecule name and ion. If prediction table provided to \code{pred_table} instead of
#' \code{param}, column names are required.
#' @slot error Numeric value - error used to create window for matching. mz values
#' will be matched against theoretical mzs +/- error.
#' @slot error_units Units for error - can be 'ppm' or 'Da'
#'
#' @examples
#'
#' #with prediction parameters
#' gpp <- glycoPredictParam(dp = c(1, 8), modifications = "deoxy", polarity = "pos", naming = "IUPAC")
#' annotated_data <- glycoAnnotate(data = data, param = gpp, error = 1.5, units = 'ppm', collapse = T)
#'
#' #with prediction table
#' gpp <- glycoPredictParam(dp = c(1, 8), modifications = "deoxy", polarity = "pos",  naming = "IUPAC")
#' pred_table <- glycoPredict(param = gpp)
#' annotated_data <- glycoAnnotate(data = data, pred_table = pred_table, error = 1.5, units = 'ppm', collapse = T, collapse_columns = c("IUPAC name", "ion"))
#'
#' @seealso \link[GlycoAnnotateR]{glycoPredictParam}
#' @seealso \link[GlycoAnnotateR]{glycoPredict}

glycoAnnotate <- function(data,
                          mz_column = 'mz',
                          mzmin_column = NULL,
                          mzmax_column = NULL,
                          pred_table = NULL,
                          param = NULL,
                          collapse = F,
                          collapse_columns = NULL,
                          error = 3,
                          error_units = 'ppm'){
  #run checks on input parameters
  if (!is.data.frame(data)){
    stop("data is not a dataframe!")
  }
  if (!mz_column %in% names(data)){
    stop("mz_column is not a column name in data.",
         " double check that mz_column is correct!")
  }
  if (!is.null(mzmin_column)){
    if(!mzmin_column %in% names(data)){
      stop("mzmin_column is not a column name in data.",
           " double check that mzmin_column is correct!")
    }
  }
  if (!is.null(mzmax_column)){
    if(!mzmax_column %in% names(data)){
      stop("mzmax_column is not a column name in data",
           ". double check that mzmax_column is correct!")
    }
  }
  if (!is.null(pred_table) & !is.null(param)){
    stop("pred_table and param supplied.",
         " please provide ONE type of input for annotation")
  }
  if (is.null(pred_table) & is.null(param)){
    stop("no glycoPredictParam supplied to 'param' AND no prediction table",
         " supplied to 'pred_table'. please provide one type of input for annotation")
  }
  if (!is.null(pred_table)){
    if (!is.data.frame(pred_table)){
      stop("pred_table is not a dataframe!")
    }
    if (!"mz" %in% names(pred_table)){
      stop("pred_table does not contain column named 'mz'!")
    }
  }
  if (!is.null(param) & class(param) != "glycoPredictParam"){
    stop("param is not an object of class glycoPredictParam",
         "please use 'glycoPredictParam' to construct an object.",
         "for documentation run '?glycoPredictParam'")
  }
  if (!is.numeric(error)){
    stop("error is not numeric!")
  }
  if (!error_units %in% c('ppm', "Da")){
    stop("error_units must be 'ppm' or 'Da'!")
  }
  if (length(error_units) != 1){
    stop("error_units must be 'ppm' OR 'Da'!")
  }
  if (!is.logical(collapse)){
    stop("collapse is not a logical!")
  }
  if (!is.null(collapse_columns) & !is.null(pred_table)){
    if(!all(collapse_columns %in% names(pred_table))){
      stop("collapse_columns are not column names in pred_table!")
    }
  }
  if (!is.null(collapse_columns) & is.null(pred_table)){
    message("warning: collapse_columns provided but no pred_table...",
            "these must correspond to columns in the table newly generated",
            "by glycoPredict!")
  }
  if (is.null(collapse_columns) & !is.null(pred_table) & isTRUE(collapse)){
    stop("pred_table provided and collapse is TRUE but no collapse_columns",
          " provided. please indicate which columns should be pasted together",
          " before collapsing (and collapsed) - for example the annotation name",
          " and ion column names")
  }

  if(!is.null(collapse_columns) & isFALSE(collapse)){
    message('collapse_columns provided but collapse is FALSE, no collapse',
            'will be performed')
  }

  if(!is.null(param)){
    if(param@format != "long"){
      message('change "format" to long in param!')}
  }

  #run glycoPredict
  if (!is.null(param)){
    message("Starting glycoPredict to generate possible annotations")
    pred_table <- GlycoAnnotateR::glycoPredict(param = param)

    if(isTRUE(collapse)){
      if(!is.null(collapse_columns)){
        if(!all(collapse_columns %in% names(pred_table))){
          stop("collapse_columns are not columns in the generated prediction table.",
               "either remove collapse_columns or ensure they match columns!")}

      }
    }
  }

  #generate mzmin and mzmax columns in pred_table
  if(error_units == 'ppm'){
    ppm_to_mz = function(mz, noise){
      ppm = mz / 1000000 * noise
      return(ppm)
    }
    pred_table <- pred_table %>%
      dplyr::mutate(mzmin_match = mz - ppm_to_mz(mz, error),
                    mzmax_match = mz + ppm_to_mz(mz, error))
  }
  if(error_units == 'Da'){
    pred_table <- pred_table %>%
      dplyr::mutate(mzmin_match = mz - error,
                    mzmax_match = mz + error)
  }

  #rename mz column in pred_table
  names(pred_table)[names(pred_table) == 'mz'] <- 'mz_pred'

  #run annotation
  message("Starting annotation with predictions against data")
  if(!is.null(mzmin_column) & !is.null(mzmax_column)){
    data$mzmin_match <- data[, mzmin_column]
    data$mzmax_match <- data[, mzmax_column]

    data.table::setDT(data)
    data.table::setDT(pred_table)
    data.table::setkey(pred_table, mzmin_match, mzmax_match)

    data_annot <- data.table::foverlaps(data, pred_table)
  }
  if(is.null(mzmin_column) & is.null(mzmax_column)){
    data <- data %>%
      dplyr::mutate(mzmin_match = get(mz_column),
                    mzmax_match = get(mz_column))

    data.table::setDT(data)
    data.table::setDT(pred_table)
    data.table::setkey(pred_table, mzmin_match, mzmax_match)

    data_annot <- data.table::foverlaps(data, pred_table)
  }

  #calculate mass error
  if(mz_column == 'mz'){
    data_annot <- data_annot %>%
      dplyr::mutate(mass_error_ppm = abs(mz - mz_pred)/mz_pred*1e6)
  }
  if(mz_column != 'mz'){
    data_annot <- data_annot %>%
      dplyr::mutate(mass_error_ppm = abs(get(mz_column) - mz_pred)/mz_pred*1e6)
  }

  #collapse annotations
  data.table::setDF(data_annot)
  if(isTRUE(collapse) & nrow(data_annot) > nrow(data)){
    message("Collapsing annotations")

    #add annotation column that is pasted together for collapsing
    if (is.null(collapse_columns)){
      if(length(param@naming) == 1){
        collapse_columns = c(paste(param@naming, "name"), "ion")
      }
      if(length(param@naming) > 1){
        collapse_columns = c(paste(param@naming[1], "name"), "ion")
      }
    }

    data_annot <- data_annot %>%
      dplyr::mutate(annotations = paste0(apply(data_annot[collapse_columns], 1,
                                               paste, collapse=':')))
    group_column_names <- c(setdiff(names(data_annot), names(pred_table)), 
                            'mass_error_ppm')
    group_column_names <- group_column_names[group_column_names != "annotations"]
    data_annot <- data_annot %>%
      dplyr::group_by(across(all_of(group_column_names))) %>%
      dplyr::summarise(annotations = toString(annotations)) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(across(all_of(c(group_column_names, "annotations")))) %>%
      dplyr::mutate(annotations = sub('NA:NA', NA, annotations))

  }

  #format final df
  data_annot <- data_annot %>%
    dplyr::select(!any_of(c('mzmin_match', 'mzmax_match',
                            'i.mzmin_match', 'i.mzmax_match')))


  return(data_annot)
}

#' Collapse annotated m/z values to have one row per peak/feature
#'
#' @include setClass.R
#' @include glycoPredict.R
#'
#' @description \code{glycoAnnotationsCollapse()} collapses the output of
#' \link[GlycoAnnotateR]{glycoAnnotate} in the case of multiple annotations
#' per peak or feature so that there is one row per peak/feature with
#' multiple annotations comma-separated.
#'
#' @export
#'
#' @slot annotated_data Dataframe annotated by \link[GlycoAnnotateR]{glycoAnnotate}
#' that has NOT been collapsed and has multiple annotations per peak/feature.
#' @slot collapse_columns Names of columns to be pasted together before collapsing.
#' Suggested is molecule name and ion.
#' @slot noncollapse_columns Names of columns that uniquely identify peaks and
#' that should be retained after collapsing - these are generally the column
#' names of your input dataframe before annotation.
#'
#' @examples
#' #annotate dataframe
#' gpp <- glycoPredictParam(dp = c(1, 8), modifications = "deoxy", polarity = "pos", naming = "IUPAC")
#' annotated_data <- glycoAnnotate(data = data, param = gpp, error = 1.5, units = 'ppm', collapse = F)
#'
#' #collapse multiple annotations
#' annotated_data_collapsed <- glycoAnnotationsCollapse(annotated_data = annotated_data, collapse_columns = c('IUPAC name', 'ion'), noncollapse_columns = c('mz', 'rt', 'sampleA', 'sampleB'))
#'
#' @seealso \link[GlycoAnnotateR]{glycoPredictParam}
#' @seealso \link[GlycoAnnotateR]{glycoPredict}
#' @seealso \link[GlycoAnnotateR]{glycoAnnotate}


glycoAnnotationsCollapse <- function(annotated_data,
                                collapse_columns,
                                noncollapse_columns){
  #run checks on input parameters
  if (!is.data.frame(annotated_data)){
    stop("annotated_data is not a dataframe!")
  }
  if(!all(collapse_columns %in% names(annotated_data))){
      stop("collapse_columns are not column names in annotated_data!")
  }
  nrow_distinct = dplyr::distinct(annotated_data,
                                  dplyr::across(dplyr::all_of(noncollapse_columns))) %>%
    nrow()
  nrow = nrow(annotated_data)
  if(nrow_distinct == nrow){
    message("all rows in annotated_data are distinct... no collapsing necessary :)")
  }

  #collapse annotations
  data.table::setDF(annotated_data)
  message("Collapsing annotations")

  annotated_data_collapsed <- annotated_data %>%
    dplyr::mutate(annotations = paste0(apply(annotated_data[collapse_columns], 1,
                                             paste, collapse=':'))) %>%
    dplyr::group_by(across(all_of(noncollapse_columns))) %>%
    dplyr::summarise(annotations = toString(annotations)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(across(all_of(c(noncollapse_columns, "annotations"))))

  return(annotated_data_collapsed)

}





