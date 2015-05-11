propInRope  <- function(stanfit,ropemin=-.05,ropemax=.05,digits=2,general=F,pars=NULL){
  
if(general){
  fit.df = data.frame(extract(stanfit)[stanfit@model_pars[stanfit@model_pars %in% pars]])
}else{
  fit.df = data.frame(extract(stanfit)[stanfit@model_pars[grep('beta_',x = stanfit@model_pars)]])
}
  
 result <-apply(fit.df,2,FUN = function(x){
    round(mean(ifelse(x > ropemin & x<ropemax,1,0)),digits)
    
  })
 
 return(result)
}