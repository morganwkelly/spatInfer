utils::globalVariables(c("X", "Y",
                         "Clusters", "Likelihood", "R2", "est_p", "estimate", 
                         "index", "name", "sim_05", "sim_p",
                         "term", "width_ci", "wts", "SE", "pseudo_se")
                       )
# scpcR::scpc() needs geodist for great-circle distances between longitude and latitude
# coordinates, but only suggests it; importing it here makes sure it is installed.
#' @importFrom geodist geodist
NULL
