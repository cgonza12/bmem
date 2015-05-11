
mypostplot <- function(stanfit,ropemin=-.01,ropemax=.01,digits=2, ROPEprop=.05,mytitle=NULL){
  source('~/Desktop/DissertationFinal/stanmer2.R')
  postsum <- stanmer2(stanfit,digits=digits)
 
  postsum$parameter= as.factor(row.names(postsum))
  
  # row.names(lmtestout) = NULL
  names(postsum) <-c('Expectation','SD' ,'LL','UL','parameter') 
  # postsum$parameter = factor(postsum$parameter,levels(postsum$parameter)[c(1,5,6,7,2,4,3,8,11,10,9)])
  postsum = postsum[postsum$parameter!='(Intercept)',]
  
  ropedf = as.data.frame(propInRope(stanfit,general = F,ropemin = ropemin,ropemax = ropemax,digits = digits))
  names(ropedf) = 'propInrope'
  ropedf$parameter = gsub('beta_','',row.names(ropedf))
  
  postsum <- merge(postsum,ropedf)
  postsum$inrope = ifelse(postsum$propInrope<ROPEprop,1,0)
    
  title = ifelse(mytitle == NULL,deparse(substitute(stanfit),mytitle))
  
  library(ggplot2)
  p <- ggplot(postsum,aes(x=parameter,y=(Expectation),ymin=(LL),ymax=(UL)))+
    geom_hline(yintercept=0,linetype='dashed',color='red')+
    geom_hline(yintercept=ropemax,linetype='dashed',color='grey')+
    geom_hline(yintercept=ropemin,linetype='dashed',color='grey')+
    geom_point(aes(color=inrope))+
    geom_errorbar(aes(color=inrope))+
    guides(color=F)+
    ylab('Estimate')+
    xlab('')+
    ggtitle(title)+
    coord_flip()+
    theme_bw()+
    theme(text=element_text(size=20))
  
  print(p)
  print(postsum)
}


