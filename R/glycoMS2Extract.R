#' MS/MS Spectra Extraction from annotated features
#'
#' @description This function extracts MS/MS spectra associated with
#' features annotated by \link[GlycoAnnotateR]{glycoAnnotate} function.

#' @param data_ms2 [MSnbase::MSnExp()], [MSnbase::OnDiskMSnExp()] or [xcms::XCMSnExp()]
#'  object with MS/MS spectra.
#' @param data_features [xcms::XCMSnExp()] with MS1 features defined by XCMS processing
#' that were annotated by \link[GlycoAnnotateR]{glycoAnnotate}.
#' @param annotations Output of \link[GlycoAnnotateR]{glycoAnnotate}. It needs to be a
#'  `data.frame` with numerical columns named "mz" and "rt" (mz columns contains *m/z*
#'  of the features, not of the annotations!).
#'
#' @return It returns a MSpectra object with all msLevel=2 spectra whose
#' precursors are the features annotated by \link[GlycoAnnotateR]{glycoAnnotate} function.
#'
#' @export
#'
#' @seealso [GlycoAnnotateR::glycoPredict()]
#' @seealso [GlycoAnnotateR::glycoPredictParam()]
#'

glycoMS2Extract = function(data_ms2, data_features, annotations){

  # Check if "data_ms2" has msLevel = 2 spectra.
  if(all(table(data_ms2@featureData@data$msLevel, data_ms2@featureData@data$fileIdx)[2,]==0)){
    stop("Error: All of the files do not have MS level 2 data.")
  }

  # Change "data_ms2" into "XCMSnExp"
  data_ms2 = as(data_ms2, "XCMSnExp")

  # Overwrite "data_ms2" peaks and features for the ones in "data_features"
  chromPeaks(data_ms2) = xcms::chromPeaks(data_features)
  featureDefinition = xcms::featureDefinitions(data_features)
  featureDefinition_filtered = featureDefinition[paste0(featureDefinition$mzmed,"_", featureDefinition$rtmed) %in%
                                                   paste0(annotations$mz,"_",annotations$rt),]

  featureDefinitions(data_ms2) = featureDefinition_filtered

  # Get msLevel=2 spectra that is related to defined features.
  MS2Spectra = xcms::featureSpectra(data_ms2, msLevel = 2, expandMz = 0.005)

  return(MS2Spectra)
}
