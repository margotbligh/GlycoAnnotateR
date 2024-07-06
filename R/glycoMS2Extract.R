#' MS/MS Spectra Extraction from annotated features
#'
#' @description This function extracts MS/MS spectra associated with
#' features annotated by \link[GlycoAnnotateR]{glycoAnnotate} function.

#' @param data_ms2 [MSnbase::MSnExp()], [MSnbase::OnDiskMSnExp()] or [xcms::XCMSnExp()]
#'  object with MS/MS spectra.
#' @param data_features [xcms::XCMSnExp()] with MS1 features or peaks defined by XCMS processing
#' that were annotated by \link[GlycoAnnotateR]{glycoAnnotate}.
#' @param annotations Output of \link[GlycoAnnotateR]{glycoAnnotate}. It needs to be a
#'  `data.frame` with numerical columns named "mz" and "rt" (mz column contains *m/z*
#'  of the features, not of the annotations!).
#' @param processing_level Is MS1 data processed to the level of peaks ('peaks')
#'  or have peaks been grouped inro features ('features', default)
#'
#' @return It returns a MSpectra object with all msLevel=2 spectra whose
#' precursors are the features annotated by \link[GlycoAnnotateR]{glycoAnnotate} function.
#'
#' @export
#'
#' @seealso [GlycoAnnotateR::glycoPredict()]
#' @seealso [GlycoAnnotateR::glycoPredictParam()]
#'

glycoMS2Extract = function(data_ms2, data_features, annotations, processing_level = 'features' ){

  # Check if "data_ms2" has msLevel = 2 spectra.
  if(all(table(data_ms2@featureData@data$msLevel, data_ms2@featureData@data$fileIdx)[2,]==0)){
    stop("Error: All of the files do not have MS level 2 data.")
  }
  
  #Check that processing level is provided
  if(!processing_level %in% c('peaks', 'features')){
    stop('Error: processing_level option not valid. Must be either peaks or features')
  }

  # Change "data_ms2" into "XCMSnExp"
  data_ms2 = as(data_ms2, "XCMSnExp")
  
  #For features
  if(processing_level == 'features'){
    #Check data type is correct
    if(!hasFeatures(data_features)){
      stop('Error: processing_level chosen is features but data_features has no features')
    }
    # Overwrite "data_ms2" peaks and features for the ones in "data_features"
    chromPeaks(data_ms2) = xcms::chromPeaks(data_features)
    featureDefinition = xcms::featureDefinitions(data_features)
    featureDefinition_filtered = featureDefinition[paste0(featureDefinition$mzmed,"_", featureDefinition$rtmed) %in%
                                                     paste0(annotations$mz,"_",annotations$rt),]
    featureDefinitions(data_ms2) = featureDefinition_filtered
    # Get msLevel=2 spectra that is related to defined features.
    MS2Spectra = xcms::featureSpectra(data_ms2, msLevel = 2, expandMz = 0.005)
  }
  
  #For peaks
  if(processing_level == 'peaks'){
    #Check data type is correct
    if(!hasChromPeaks(data_features)){
      stop('Error: processing_level chosen is peaks but data_features has no peaks')
    }
    #Filter chromPeaks to only those annotated
    chromPeaks_df <- chromPeaks(data_features) %>% as.data.frame()
    chromPeaks_df_filtered <- chromPeaks_df[paste0(chromPeaks_df$mz,"_", chromPeaks_df$rt) %in%
                                              paste0(annotations$mz,"_",annotations$rt),]
    # Assigned filtered chromPeaks to "data_ms2"  as peaks
    chromPeaks(data_ms2) = chromPeaks_df_filtered
    
    # Get msLevel=2 spectra that is related to peaks.
    MS2Spectra = xcms::chromPeakSpectra(data_ms2, msLevel = 2, expandMz = 0.005)
  }


  return(MS2Spectra)
}
