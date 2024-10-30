#' Chetty et al data on intergenerational mobility for US Census Zones.
#'
#' Variables from Table VI, and data exclude Alaska and Hawaii.
#'
#'
#' @format ## `opportunity`
#' A data frame with 693 rows and 10 variables.
#' \describe{
#'   \item{state_id}{State abbreviation}
#'   \item{mobility}{Absolute upward mobility: e_rank_b.}
#'   \item{short_commute}{Fraction short commute: frac_traveltime_lt15.}
#'   \item{single_mom}{Fraction single mothers: cs_fam_wkidsinglemom.}
#'   \item{gini}{Gini bottom 99: gini99.}
#'   \item{dropout_rate}{High school dropout rate: dropout_r.}
#'   \item{social_cap}{Social capital index: scap_ski90pcm.}
#'   \item{dropout_na}{Dummy for missing dropout rate.}
#'   \item{X}{Longitude.}
#'   \item{Y}{Latitude}

#'   ...
#' }
#' @source Chetty et al (2014) "Where is the Land of Opportunity? The Geography of Intergenerational Mobility in the United States."
#' Quarterly Journal of Economics, 129, 1553–1623.
"opportunity"
