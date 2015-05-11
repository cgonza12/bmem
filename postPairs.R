
postPairs<- function(stanfit,save=F,width = 600, height=500,
                     name=stanfit,corrOnly = F,digits=2,suppress=.3){
source('~/Desktop/DissertationFinal/panelfuns.R')
fit.df = data.frame(extract(stanfit)[stanfit@model_pars[grep('beta_',x = stanfit@model_pars)]])
names(fit.df) = gsub('beta_','',names(fit.df))
if(corrOnly){
  
  r <- as.matrix(round(abs(cor(fit.df)),digits))
  r = as.data.frame(r)
  r[r<suppress] <- '--'
  
  print(r)
  sink('tempcorr.txt')
  print((r))
  sink()
}else{

if(save){
  if(is.character(name)){
    model_name = name
  }else{
    
    model_name = deparse(substitute(name))
  }
  png(paste0(model_name,'_postPairs.png'),width=width,height=height,units="px")
  pairs(na.omit(fit.df),
        lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
  dev.off() 
}else{
  options(device = "X11")
  dev.new()
pairs(na.omit(fit.df),
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
  options(device = "RStudioGD")

}
}
}