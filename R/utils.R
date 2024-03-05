#' Convert from ppm to mz
#'
#' @description Simple function to calculate from ppm to mz, for a given mz and 
#' amount of ppm.
#'
#' @slot mz *m/z* value to calculate for
#' @slot ppm the ppm amount to calculate
#' 
#' @export
#'
#' @examples
#' 
#' #calculate 2 ppm of m/z value of 259.0129
#' 
#' ppm_to_mz(259.0129, 2) 
#'

ppm_to_mz = function(mz, ppm){
  ppm_of_mz = mz / 1000000 * ppm
  return(ppm_of_mz)
}