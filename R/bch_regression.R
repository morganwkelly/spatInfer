#' Carry out spatial basis regressions with BCH standard errors.
#'
#' @param fm  The equation formula. The variable of interest is the first one on the right hand side.
#' @param df  The name of the dataset use. Longitude and latitude must be named X and Y and have no missing values.
#' @param splines The dimension of the linear tensor used, from optimal_basis function.
#' @param pc_num The number of principal component to include, again from optimal_basis.
#' @param clusters Number of k-medoids clusters to use based on the placebo test.
#' @param weights The name of the weighting variable in the dataset, as a
#'   string, for example `weights="pop"`. Leave it out for an unweighted
#'   regression.
#' @param cov Defaults to BCH. It gives heteroskedasticity consistent standard errors otherwise.
#' @param k_medoids Use k-medoids clustering (PAM), the default. Set to `FALSE`
#'   to use Clara, as for large datasets: use the same setting as in the placebo
#'   test so that the clusters are the same.
#'
#' @return feols object that can be printed and exported using the modelsummary package.
#' Given that t-statistics from large cluster regressions are not immediately interpretable due to their low degrees
#' of freedom the regression table gives a confidence interval and p-value, although this is trivial to change in
#' the modelsummary output.
#' @export
#'
#' @examples
#' library(spatInfer)
#' data(opportunity)
#' # Choose a sample of 250 observations and run only 100 placebo simulations
#' # to speed things up.
#' set.seed(123)
#' opportunity=opportunity |> dplyr::slice_sample(n=250)
#' # Use the number of splines and PCs indicated by optimal_basis and
#' # the number of clusters from the placebo test.
#' ck=bch_regression(mobility~single_mothers+short_commute+gini+dropout_rate+social_cap,  opportunity,
#' splines=4,pc_num=3,
#' clusters=5,cov="BCH")
#'
#' #Now generate a regression output table.
#' modelsummary::modelsummary(list(CK=ck),
#' statistic = c("conf.int","p = {p.value}"),
#' coef_omit = c("Intercept|PC*"), #omit basis and intercept
#' gof_map = c("nobs", "r.squared"),fmt=2)

bch_regression=function(fm,df,splines,pc_num,clusters,weights=NULL,cov="BCH",k_medoids=TRUE){

  if(is.null(df$X)|is.null(df$Y))
    stop("You must have longitude and latitude variables named X and Y")
  if(sum(is.na(df$X))>0|sum(is.na(df$Y))>0)
    stop("You cannot have missing values in longitude and latitude.")

#regression weights, taken from the variable named by weights
wts=NULL
if(!is.null(weights)){
  df$wts=regression_weights(df,weights)
  wts=~wts
}

new_names=set_names(fm)
orig_name= new_names$orig_name  #stringr::str_split_1(new_names$gg[2],"\\+")[1]
df1=df |> dplyr::rename(dep_var=new_names$gg[1]) |>
  dplyr::mutate(explan_var=orig_name)
rhs=new_names$rhs

#get the principal components that minimise BIC and add them to the dataset.
pc=prin_comp(df1,splines,pc_num)
df1=cbind.data.frame(df1,pc)

# df1=df1 |>
#   rename(orig_name=explan_var)

#add pcs to regression equation
eq_est= paste0("dep_var~",orig_name,rhs)
eq_est=as.formula(paste(eq_est,paste(names(pc),collapse="+"),sep="+"))



if(cov=="BCH"){
  Coords=as.matrix(df |> dplyr::select(X,Y))
  clust_bch=make_clusters(Coords,clusters,k_medoids)    #BCH clusters, as in placebo()
  CK=fixest::feols(eq_est,
                   data=df1,
                   weights=wts,
           cluster=clust_bch,
           data.save=TRUE
  )}else{
  CK=fixest::feols(eq_est,
                   data=df1,
                   weights=wts,
           vcov="hetero",
           data.save=TRUE
  )}
return(CK)

}
