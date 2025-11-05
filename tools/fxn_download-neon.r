#' @title Download NEON Data
#' 
#' @description Download all NEON data for user-defined DPI code. Uses `neonUtilities::loadByProduct` so you must have already done the necessary authentication. Note that some arguments of `loadByProduct` are hard-coded within this function and not editable.
#' 
#' @param dpi (character) 13-digit NEON data product ID
#' @param start_yr (numeric) 4-digit year at start of desired time range
#' @param end_yr (numeric) 4-digit year at end of desired time range. Must be larger than 'start_yr'
#' @param prefix (character) text to add to start of exported files. Year added to end of file automatically
#' @param dest (character) local file path to which to export data
#' @param redownload (logical) whether to redownload data for a given year if they already are found in the path specified by 'dest'
#' @param wanted_files (character) name of list elements (returned by `neonUtilities::loadByProduct`) to export
#' @param quiet (logical) whether to print a message for each year as it is downloaded or skipped because it already exists locally
#' 
download_neon <- function(dpi = NULL, start_yr = NULL, end_yr = NULL,
                          prefix = "neon_", dest = getwd(), redownload = FALSE,
                          wanted_files = NULL, quiet = FALSE){

  # Error checks for DPI
  if(is.null(dpi) || is.character(dpi) != TRUE || length(dpi) != 1 || nchar(dpi) != 13)
    stop("'dpi' must be specified as a single 13-character vector")

  # Error checks for year arguments
  if(is.null(start_yr) || is.numeric(start_yr) != TRUE || length(start_yr) != 1 || nchar(start_yr) != 4)
    stop("'start_yr' but be a single, 4-digit number")
  if(is.null(end_yr) || is.numeric(end_yr) != TRUE || length(end_yr) != 1 || nchar(end_yr) != 4)
    stop("'end_yr' but be a single, 4-digit number")
  if(start_yr > end_yr)
    stop("'start_yr' must be before 'end_yr'")

  # Error checks for wanted files
  if(is.null(wanted_files) || is.character(wanted_files) != TRUE)
    stop("'wanted_files' must be specified as a character vector")
  
  # Warning for malformed 'redownload'
  if(is.null(redownload) || is.logical(redownload) != TRUE){
    warning("'redownload' must be specified as a logical. Coercing to TRUE")
    redownload <- TRUE }

  # Warning for malformed 'quiet'
  if(is.null(quiet) || is.logical(quiet) != TRUE){
    warning("'quiet' must be specified as a logical. Coercing to FALSE")
    quiet <- FALSE }
  
  # Loop across years in range
  for(focal_yr in start_yr:end_yr){
    # focal_yr <- start_yr
  
    # Identify any data for the focal year that we already have
    local_data <- dir(path = dest, pattern = paste0("_", focal_yr, ".csv"))
  
    # Skip this year if we already have the expected number of files
    ## Note that the expected number of local files is hard-coded
    if(all(paste0(prefix, wanted_files, "_", focal_yr, ".csv") %in% local_data) & redownload == FALSE){ 
      
      # Completion message
      message("Data for ", focal_yr, " already downloaded")
  
    # Otherwise, download the data!
    } else {
      
      # Progress message
      message("Acquiring data for ", focal_yr)
  
      # Grab start/end months for the focal year
      focal_start <- paste0(focal_yr, "-01")
      ## End year we need to be savvy about if the final year is this year (and thus incomplete)
      if(focal_yr != stringr::str_sub(Sys.Date(), start = 1, end = 4)){
        focal_end <- paste0(focal_yr, "-12")
      } else { focal_end <- stringr::str_sub(Sys.Date(), start = 1, end = 7) }
  
      # Identify relevant plant data
      focal_list <- neonUtilities::loadByProduct(dpID = dpi, 
        startdate = focal_start, 
        enddate = focal_end, 
        check.size = F)  
  
      # Export wanted files
      for(focal_wanted in wanted_files){
        # focal_wanted <- "div_1m2Data"
  
        # Actually export files
        write.csv(x = focal_list[[focal_wanted]], na = '', row.names = F,
                  file = file.path(dest, paste0(prefix, focal_wanted, "_", focal_yr, ".csv")))
  
      } # Close export loop
    } # Close download conditional
  } # Close year loop
} # Close function

# End ----
