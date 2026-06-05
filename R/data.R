#' HBomics dummy data set
#'
#' example clonal SNVs data from COLO320-DM/HSR treated with STP,
#' used to demonstrate the HBomics pipeline.
#'
#' @format A metasheet data frame with 1412 rows and 13 variables from COLO320-DM:
#' @source For tutorial purposes.
"meta_d"

#' @format A metasheet data frame with 1412 rows and 13 variables from COLO320-HSR:
#' @source For tutorial purposes.
"meta_h"

#' Stan Data List of COLO320-DM
#'
#' @format A list with elements required by the Stan model (T, C, S, A, N, etc.)
"data_list_d"

#' Stan Data List of COLO320-HSR
#'
#' @format A list with elements required by the Stan model (T, C, S, A, N, etc.)
"data_list_h"

#' Phased SNV Data from eagle2
#'
#' @format A data frame with 1412 rows and 2 variables:
"phased"