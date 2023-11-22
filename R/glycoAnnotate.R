#' Annotate m/z values with glycan compositions based on theoretical values
#'
#' @include setClass.R
#' @include glycoPredict.R
#'
#' @description \code{glycoAnnotate()} annotates peaks or features in MS data,
#' using either a pre-generated table by \link[GlycoAnnotateR]{glycoPredict} or
#' by generating a new table. 
#'
#' @export
#' 
#' @slot data Dataframe containing data to be annotated.
#' @slot mz_column Name of column containing m/z values.
#' @slot mzmin_column
#' @slot mzmax_column
#' @slot pred_table description
#' @slot param description
#' @slot collapse description
#' 
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
#' \code{glycoPredict()} is used to predict masses and mass to charge ratios of all theoretically 
#' possible glycans within a set of constraining parameters (defined in the 
#' \code{glycoPredictParam} object). This package was written 
#' for annotation of mass spec data (especially LC-MS) but if used for 
#' other purposes either ionisation mode and very wide scan ranges can be given. 
#' The function works by sourcing a python file and then using the function 
#' encoded in the python script.
#' 
#' @seealso 
#' glycoAnnotateR::glycoPredictParam()
#' 

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
  if (!is.null(collapse_columns) & !is.null(pred_table) & 
      !collapse_columns %in% names(pred_table)){
    stop("collapse_columns are not column names in pred_table!")
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
  
  #run glycoPredict
  if (!is.null(param)){
    message("Starting glycoPredict to generate possible annotations")
    pred_table <- GlycoAnnotateR::glycoPredict(param = param)
    
    if(!collapse_columns %in% names(pred_table))
      stop("collapse_columns are not columns in the generated prediction table.",
           "either remove collapse_columns or ensure they match columns!")
  }
  
  #generate mzmin and mzmax columns in pred_table
  if(error_units == 'ppm'){
    ppm_to_mz = function(mz, noise){
      ppm = mz / 1000000 * noise
      return(ppm)
    }
    pred_table <- pred_table %>% 
      dplyr::mutate(mzmin = mz - ppm_to_mz(mz, error),
                    mzmax = mz + ppm_to_mz(mz, error))
  }
  if(error_units == 'Da'){
    pred_table <- pred_table %>% 
      dplyr::mutate(mzmin = mz - error,
                    mzmax = mz + error)
  }
  
  #run annotation
  message("Starting annotation with predictions against data")
  if(!is.null(mzmin_column) & !is.null(mzmax_column)){
    if (mzmin_column != "mzmin"){
      names(data)[names(data) == mzmin_column] <- "mzmin"
    }
    if (mzmax_column != "mzmax"){
      names(data)[names(data) == mzmax_column] <- "mzmax"
    }
    data.table::setDT(data)
    data.table::setDT(pred_table)
    data.table::setkey(pred_table, mzmin, mzmax)
    
    data_annot <- data.table::foverlaps(data, pred_table)
  }
  if(is.null(mzmin_column) & is.null(mzmax_column)){
    data <- data %>% 
      dplyr::mutate(mzmin = get(mz_column),
                    mzmax = get(mz_column))

    data.table::setDT(data)
    data.table::setDT(pred_table)
    data.table::setkey(pred_table, mzmin, mzmax)
    
    data_annot <- data.table::foverlaps(data, pred_table)
  }
  
  #collapse annotations
  setDF(data_annot)
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
    group_column_names <- setdiff(names(data_annot), names(pred_table))
    group_column_names <- group_column_names[group_column_names != "annotations"]
    data_annot <- data_annot %>%
      dplyr::group_by(across(all_of(group_column_names))) %>%
      dplyr::summarise(annotations = toString(annotations)) %>%
      ungroup() %>%
      distinct(across(all_of(c(group_column_names, "annotations"))))

  }
  
  #format final df
  if(isFALSE(collapse)){
    data_annot <- data_annot %>% 
      select(!c('mzmin', 'mzmax'))
  }
  if('mz' %in% names(pred_table) & 'mz' %in% names(data_annot)){
    data_annot <-  data_annot %>% 
      dplyr::rename(mz_pred = mz)
    if(isFALSE(collapse)){
      data_annot <-  data_annot %>% 
        dplyr::rename(mz = `i.mz`)
    }
  }
  if(!is.null(mzmin_column) & "i.mzmin" %in% names(data_annot)){
      names(data_annot)[names(data_annot) == "i.mzmin"] <-  mzmin_column
  }
  if(!is.null(mzmax_column) & "i.mzmax" %in% names(data_annot)){
    names(data_annot)[names(data_annot) == "i.mzmax"] <-  mzmin_column
  }
  
}







