#' Frankfurt airport precipitation data
#' 
#' @description 
#' Accumulated 06-30 hour precipitation observations and operational ECMWF
#' ensemble forecasts for Frankfurt airport, Germany. The observations are
#' airport station observations (WMO station index 10637), the forecasts are
#' gridded forecasts for the 0.25 degrees latitude/longitude box containing the
#' station. Dates range from 2007-01-01 to 2017-01-01, days with missing values
#' have been removed.
#' 
#' @usage 
#' data("rain")
#' 
#' @format
#' A data frame with 3617 rows. The first column gives the dates, the second 
#' column are the observations. The remaining columns are the ensemble forecasts
#' (high resolution HRES, 50 perturbed forecasts P1 to P50 and the control 
#' forecast CTR for the perturbed forecasts). The units of the forecasts and
#' observations are mm/m^2.
#' 
#' @source 
#' Observations: \url{http://www.ogimet.com/synops.phtml.en}
#' 
#' Forecasts: available on TIGGE
#' \url{https://confluence.ecmwf.int/display/TIGGE/TIGGE+archive}
#' 
#' @references
#' Bougeault et al. (2010) The THORPEX Interactive Grand Global Ensemble. Bull.
#' Amer. Meteor. Soc., 91, 1059-1072.
#'
#' Swinbank et al. (2016) The TIGGE project and its achievements. Bull. Amer.
#' Meteor. Soc., 97, 49-67.
"rain"