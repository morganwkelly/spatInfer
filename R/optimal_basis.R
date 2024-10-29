#' Plot BIC and R2 for different spatial bases
#'
#' `optimal_basis()` gives a plot of BIC and adjusted R-squared against number of
#' principal components for tensor bases starting at 3x3. It reports the combination
#' that minimizes BIC.
#'
#' @param fm  The regression formula.
#' @param df  The dataset used.
#' @param max_splines The maximum dimension of tensor to examine
#' @param Description Optional description of the regression used to label the output diagram.
#'
#' @return ggplot object: graph of BIC and R2 of regression of dependent variable on bases.
#' @export
#'
#' @examples
#' library(spatInfer)
#' data(opportunity)
#' #Choose a sample of 100 to speed things up.
#' set.seed(123)
#' opportunity=opportunity |> dplyr::slice_sample(n=100)
#' optimal_basis(mobility~racial_seg+single_mom,opportunity,
#'     max_splines=6, Description="Intergenerational Mobility.")



optimal_basis=function(fm,df,max_splines,Description=""){

  #########this picks the optimal number of spatial basis principal components to use by minimizing BIC

  if(max_splines>12)
    stop("The maximum number of splines you can use is 12.")

  fm=deparse1(fm)
  fm=stringr::str_replace_all(fm," ","")
  gg=stringr::str_split_1(fm,"~")
  rhs=stringr::str_split(gg[2],"\\+",n=2)
  rhs=unlist(rhs)[2]
  if(is.na(rhs)){
    rhs=""}else{
      rhs=paste0("+",rhs)
    }
  fm=paste0("dep_var~explan_var",rhs)
  fm=as.formula(fm)

  df=df |> dplyr::rename(dep_var=gg[1],
                   explan_var=stringr::str_split_1(gg[2],"\\+")[1])

  ####get principal components of tensors spline and add to dataset
  bic_results=list()
  mx=max_splines-2

############3x3 linear spline
  gm_2=mgcv::bam(dep_var~
                   te(X,Y,bs=c("bs"),
                      k=3,
                      m=1),
                 data=df,
                 discrete=T)
  pc=prcomp(model.matrix(gm_2))
  pc=cbind.data.frame(df$dep_var,pc$x)
  names(pc)[1]="dep_var"
  pc_fit=list()
  for(j in 1:(ncol(pc)-1)){
    ll=lm(dep_var~.,pc[,1:(j+1)])
    pc_fit[[j]]=data.frame(BIC=BIC(ll),R2=summary(ll)$adj.r.squared)
  }
  pc_fit=purrr::list_rbind(pc_fit) |> as.data.frame()
  pc_fit$index=1:nrow(pc_fit)
  bic_results[[1]]=  pc_fit   #data.frame(pc_fit$BIC)
  names(bic_results[[1]])=c(paste0("BIC_",3),paste0("R2_",3),"index")


bas=bic_results[[1]]


###########splines from 4x4 up: linear
  for (i in 2:mx){
    spl=i+2
  gm_2=mgcv::bam(dep_var~
             te(X,Y,bs=c("bs"),
                k=spl,
                m=1),
           data=df,
           discrete=T)
  pc=prcomp(model.matrix(gm_2))
  pc=cbind.data.frame(df$dep_var,pc$x)
  names(pc)[1]="dep_var"
  pc_fit=list()
  for(j in 1:(ncol(pc)-1)){
    ll=lm(dep_var~.,pc[,1:(j+1)])
   pc_fit[[j]]=data.frame(BIC=BIC(ll),R2=summary(ll)$adj.r.squared)
  }
  pc_fit=purrr::list_rbind(pc_fit) |> as.data.frame()
  pc_fit$index=1:nrow(pc_fit)
bic_results[[i]]=  pc_fit   #data.frame(pc_fit$BIC)
names(bic_results[[i]])=c(paste0("BIC_",i+2),paste0("R2_",i+2),"index")
}

  #bas=bic_results[[1]]
    for (i in 2:mx){
  bas=dplyr::full_join(bas,bic_results[[i]],by="index")
    }
bas=bas |> dplyr::relocate(index)

  best=bas |> dplyr::select(index,dplyr::starts_with("BIC")) |>
    tidyr::pivot_longer(-index,values_to = "BIC")|>
    dplyr::arrange(BIC) |>
    dplyr::slice(1) #|> dplyr::select(index) |> as.numeric()

  best_spline=stringr::str_split_1(best$name,"_")[2]
  best=best$index

  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#CC79A7", "#009E73", "#F0E442", "#0072B2", "#D55E00") #colorblind palette
  bic=bas |> dplyr::select(index,starts_with("BIC")) |>
    tidyr::pivot_longer(-index,values_to = "BIC") |>
    na.omit() |>
    ggplot2::ggplot(ggplot2::aes(index,BIC,color=name))+
    ggplot2::geom_line()+
    ggplot2::geom_vline(xintercept=best,linewidth=0.25,colour="red")+
    ggplot2::theme_bw()+
    ggplot2::scale_color_manual(values=cbPalette[1:mx])+
    ggplot2::theme(legend.title = ggplot2::element_blank())+
    ggplot2::labs(x="")
 r2= bas |> dplyr::select(index,dplyr::starts_with("R2")) |>
   tidyr::pivot_longer(-index,values_to = "R2") |>
   na.omit() |>
   ggplot2::ggplot(ggplot2::aes(index,R2,color=name))+
   ggplot2::geom_line()+
   ggplot2::geom_vline(xintercept=best,linewidth=0.25,colour="red")+
   ggplot2::labs(x="Number of Principal Components",y="Adj R2")+
   ggplot2::theme_bw()+
   ggplot2::scale_color_manual(values=cbPalette[1:mx])+
   ggplot2::theme(legend.title = ggplot2::element_blank())
  #spl_type=" Linear"

  return(
    patchwork::wrap_plots(bic,r2,nrow=2)+patchwork::plot_annotation(
    title=paste0(Description," BIC and R2 of 3x3 to ",max_splines,"x",max_splines," tensors"),
    subtitle=paste0("BIC is minimized by a ",best_spline ,"x",best_spline, " Linear Tensor with ",  best, " PCs."))
    )
}

