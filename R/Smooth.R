

#naive smooth

#to-do list: (1) now res = RES only
#(2) use my WLR code; not via fda function; and add the added-stype basisfunction and use the same style as the smooth_SOF


Smooth<-function( RES,#RES$DRres #the result summary information table RES<-getSummaryInf(rdat) #need bx bxse by byse mean est se
                  StraMet='DR', #stratification method; default is 'DR'
                  Rall=NA, #A vector: Exposure range for smoothing and visualization
                  baseline=NA, #the baseline exposure value for orginal causal effect shape
                  splinestyle='Bspline' , #spline type
                  Norder=1,
                  XYmodel='0', #tell me what the true X-Y effect shape is?
                  Knots=NA, #internal knots;  c(1,2,3) or '3' (character)
                  Lambda=0,  #tuning parameter for roughness; lambda=0 means no roughness penalty
                  random_effect=TRUE,
                  getHeterQ=TRUE, #if to show the Q and trend test p-value for the differentiate function plot?
                  Plot=FALSE,  #logic; Plot=TRUE will automatically print the two ggplots
                  ylim_used=NA  #the user-defined ylim range for h'(x) plot
){
  if(!StraMet%in%c('DR','R') ){stop('please use the correct stratification method name')}
  if(StraMet=='DR'){res<-RES$DRres;weightfun<-RES$DRweightfun}
  if(StraMet=='R'){res<-RES$Rres;weightfun<-RES$Rweightfun}

  if(getHeterQ){
    Sr<-mr_ivw(mr_input(res$bx, res$bxse ,res$by, res$byse))#naive IVW under null (no effect heterogeneity)
    #mr_ivw是免疫NA的，即，即使把NA移除掉，naive IVW estimate依旧是一样的
    delta<-Sr@Estimate #一维的
    #Cochran's Q statistic
    ssigma_square <- res$byse^2 + delta^2*res$bxse^2
    QQ<-sum( ( res$by - delta*res$bx     )^2 /ssigma_square ,na.rm=TRUE      )
    #critical value:
    Qdf<-sum(!is.na( ( res$by - delta*res$bx     )^2 /ssigma_square  )) -1 #normally Ns-1 without NA stratum
    pvalue<-1-pchisq(QQ, Qdf ) #不会考虑NA的stratum
    HeterQ<-c(QQ, Qdf, round(pvalue,3)  )#; names(HeterQ)<-c('Q statistic' ,'df' ,'p-value'  )
    #RES$RHeterQ<-HeterQ

    ###meta regression trend test (using x and x^2 up to second order)
    metadat<-data.frame(  est=res$est , se=res$se, mean=res$mean     )
    if(nrow(res)==2){#注意 rma还需要估计原文中tau^2的大小，所以存在一个额外的parameter；并不是最简单的WLR,所�?Ns=2时只用Q就好
      metares<-list(pval=rep(NA,3))
      TrendValues<-rep(NA,3)
    }else{
      #metares<-rma(yi=est, sei=se, mods =cbind(mean , mean^2 ), data = metadat)
      metares<-rma(yi=est, sei=se, mods =cbind(mean ), data = metadat)
      TrendValues<-c( metares$QM , metares$QMdf[1] , metares$QMp    )
    }
  }

  SmoothRes<-list()
  #get the minimal useful information table
  mRES<-as.data.frame(cbind( res$mean,res$est , res$se  ))
  names(mRES)<-c('mean', 'est', 'se'  )
  #romove the NA values
  mRES<-na.omit(mRES)
  Ns<-nrow(mRES) #the number of stratum after removing NA values

  if(is.na(Rall)){
    #Rall<-c(min( mRES$mean)-1, max(mRES$mean)+1    )
    Rall<-c(  min(res$min  ), max(res$max  )  )
    cat('The exposure range used: ', Rall,'\n')}

  if(is.vector(Knots)){
    internal_knots<-sort(Knots)  #sort(NA) -> logical(0); 会�?�breaks这个vector自动忽略
    if( sum( !c(Rall[1]<internal_knots,internal_knots<Rall[2]) )>0 ){
      stop('make sure that the internal knots are inside the Rall')
    }
  }

  if(is.character(Knots)){
    if(as.numeric(Knots)-round(as.numeric(Knots))==0  ){
      internal_knots<-seq(Rall[1],Rall[2],len=as.numeric(Knots)+2)[-c(1,as.numeric(Knots)+2)]
      }else{
      stop('PLease use the correctly integer number of Knots')
    }
  }
  cat('The internal knots are: ', internal_knots,' (Empty means no internal knots used)','\n')




  time_obs<-seq(Rall[1],Rall[2],len=1000) #visualization points

  if(splinestyle=='Bspline'){
    Bs<-create.bspline.basis(  c(Rall[1],Rall[2])  ,
                               breaks=c(  Rall[1],  internal_knots ,  Rall[2] ) ,
                               norder=Norder )
    # Par <- fdPar( fdobj=Bs, Lfdobj=NULL, lambda=Lambda )#都默认不用penalty;这里应该没问???
    # smoothed_LACE <- smooth.basis(  mRES$mean  , mRES$est , Par ,wtvec=1/mRES$se^2 )$fd
    # dsdf<-smooth.basis(  mRES$mean  , mRES$est , Par ,wtvec=1/mRES$se^2 )$df  # = # of interior knots + order
    # thetahat<-smoothed_LACE$coefs
    # SmoothRes$thetahat<-as.vector(thetahat)

    MMM<- eval.basis(   mRES$mean , Bs      )  #function values at the  RES$mean position (mean_exposure observed position)
    #dim(MMM) #Ns [basisfunction number]


    MM<- eval.basis(  time_obs , Bs      )  #matrix of basis function values in all positions
  }else{
    #use the added-style spline ( only can be used for Noder=1 case; i.e. constant basis function)
    if(Norder!=1){stop('other spline type only currently valid for Norder=1')}

    #MMM: matrix of basis function values in the stratum exposure mean values
    #MM: matrix of basis function values in all positions
    #dim(MM) #1000 * [basisfunction number]

    MMM<-rep( 1,length(mRES$mean) )
    for(i in 1:length(internal_knots)){
      MMM<-cbind(MMM , as.numeric(mRES$mean>=internal_knots[i]) )
    }


    MM<-rep( 1,length(time_obs) )
    for(i in 1:length(internal_knots)){
      MM<-cbind(MM , as.numeric(time_obs>=internal_knots[i]) )
    }


  }



  ###use my code with SLR, no black box
  SSigma<-diag(  1/mRES$se )#standard (non-squared) sigma; i.e. Sigma^{-1/2}
  yy<-SSigma%*%mRES$est

  XX<-SSigma%*%MMM
  SLRthetahat<-solve(  t( XX)%*%XX )%*% t( XX)%*% yy  #Simple (unweighted) linear regression

  thetahat<-as.vector(SLRthetahat)
  SmoothRes$thetahat<-thetahat

  #clear! 关键在于wtvec需要使用平方量  #checked! 没毛�?
  ###接着估计tau^2???
  #直接用原始的weighted regression - LRT:  当然用SLR也行�? yy<-SSigma%*%RES$est # XX<-SSigma%*%MMM

  temp<-MMM%*%thetahat# fitted values
  dsdf<-Norder+length(internal_knots  ) #degree-of-freedom
  if(random_effect){
    if(Ns==dsdf){
      tau2<-1   #这里主�?�是考虑到完全full model的情???(即，fractional polynomial method)，此时df分母=0 容易报错
    }else{
      tau2<-max(1,sum(  (mRES$est-temp)^2/mRES$se^2     )/(  Ns- dsdf))#没有roughness penalty，那df自然是整???
      }
    cat('Random-effect is considered for smoothing and the actual random effect squared value is: ', tau2,'\n')
  }else{
    tau2<-1
  }



  #AIC
  ###STEP5: also return the AIC (for model selection purpose)
  Sigma<-diag(    mRES$se^2  )
  SSSigma<-Sigma*tau2  #times tau for random-effect setting

  varthetahat<-solve( t(MMM) %*% solve(Sigma) %*% MMM )  * tau2  #经典WLR的系数估计量的方�?

  SmoothRes$var.matrix<-varthetahat



  AIC<- 2*log(  det(  SSSigma )      ) + t(mRES$est-MMM%*%thetahat)%*%solve(SSSigma)%*%(mRES$est-MMM%*%thetahat) + 2*dsdf
  SmoothRes$AIC<-AIC
  BIC<- 2*log(  det(  SSSigma )      ) + t(mRES$est-MMM%*%thetahat)%*%solve(SSSigma)%*%(mRES$est-MMM%*%thetahat) + log(Ns)*dsdf
  SmoothRes$BIC<-BIC
  #or: use Simple Linear Regression
  #tau2<-max(1,sum( (yy-XX%*%SLRthetahat )^2 )/(Ns -dsdf))#checked 确实是一样的


  ###Following inference


  ###WLR with fda code
  # y2cMap<- smooth.basis(  mRES$mean  , mRES$est , Par ,wtvec=1/mRES$se^2 )$y2cMap
  # varthetahat<-y2cMap %*% diag(  mRES$se^2 )  %*% t(y2cMap  )* tau2  #经典WLR的系数估计量的方�?
  #
  # SmoothRes$var.matrix<-varthetahat
  summary_table<-cbind( as.vector(thetahat),
                        sqrt(as.vector(diag(varthetahat))))  #marginal s.e.s
  #依据tau2的值来判断使用Z-tes还是t-testt；毕竟tau2==1时，相当于normal variance 不�?�估计的�????
  if(tau2==1){
    summary_table<-cbind(summary_table,2*(1-pnorm(  abs(as.vector(thetahat)/sqrt(as.vector(diag(varthetahat))))  )   ) )
  }else{
    summary_table<-cbind(summary_table,2*(1-pt(  abs(as.vector(thetahat)/sqrt(as.vector(diag(varthetahat))) ), Ns-length(thetahat)   )   )      )
    }
  colnames(summary_table)<-c('est','s.e.','p-value')
  rownames(summary_table)<-paste0('Par',1:length(thetahat))
  SmoothRes$summary<-round(summary_table,3)
  #y2cMap表达式上就是solve(  t( XX)%*%XX )%*% t( XX)%*%SSigma；含义就是WLR的hat matrix
  #MM<- eval.basis(  time_obs , Bs      )  #matrix of basis function values in all positions
  #dim(MM) #1000 dsdf #df=degree-of-freedom= # of theta parameters


  ###final functional pointwise estimator and s.e.
  varM<- MM %*% varthetahat %*%  t(  MM)  #最后一向是random-effect

  pw_std_error<-sqrt( diag(varM) )  #由于是Pointwise, 只取对�?�线元素�????
  #smoothed_LACE_<-eval.fd(    time_obs ,smoothed_LACE )

  smoothed_LACE_<-MM%*%thetahat

  # plot(time_obs,smoothed_LACE_,type='l')
  # lines(time_obs , smoothed_LACE_+1.96*pw_std_error, lty=2, lwd=1        )
  # lines( time_obs, smoothed_LACE_-1.96*pw_std_error, lty=2, lwd=1          )

  ggdata<-data.frame( tps=mRES$mean, est = mRES$est,
                      estlow=mRES$est -1.96*mRES$se,
                      estup=mRES$est +1.96*mRES$se
  )
  ggsmooth<-data.frame( time_obs=time_obs,
                        smoothed_LACE_=smoothed_LACE_,
                        smooth_low=smoothed_LACE_-1.96*pw_std_error,
                        smooth_up=smoothed_LACE_+1.96*pw_std_error)
  ggreal<-data.frame( tps=time_obs,
                      est=NaN)
  ##now assign the nonlinear_dif
  if( XYmodel=='1' ){    nonlinear_dif<-function(t){ sapply(   t , function(tt){ 0}) }     }
  if( XYmodel=='2' ){    nonlinear_dif<-function(t){0.2*t}   }
  if( XYmodel=='3' ){    nonlinear_dif<-function(t){-0.2*t*(t>0)}    }
  if( XYmodel=='4' ){    nonlinear_dif<-function(t){5.0*(t>0.6) }  }  #5.0*(t>0.6)  #45*t-30.0*t^2
  if( XYmodel%in%c('1','2','3','4') ){ ggreal$est<-nonlinear_dif(time_obs)    }

  if(is.na(ylim_used)[1]){
    ylim_used <- c( min(res$est-1.96*res$se)-max(res$se),max(res$est+1.96*res$se)+max(res$se))
  }else{
    ylim_used <- ylim_used
  }

  if(getHeterQ){
    p<-ggplot(ggdata, aes(tps, est))+
      geom_point(data=ggdata,mapping=aes(x=tps,y=est),color='black',alpha=1)+
      geom_errorbar(data=ggdata,mapping=aes(x=tps,ymin=estlow,ymax=estup),color='black',width = 0.1)+
      geom_line(data=ggsmooth,mapping=aes(x=time_obs,y=smoothed_LACE_),color='black',alpha=1)+
      geom_line(data=ggsmooth,mapping=aes(x=time_obs,y=smooth_low),color='black',alpha=1,linetype='dashed')+
      geom_line(data=ggsmooth,mapping=aes(x=time_obs,y=smooth_up),color='black',alpha=1,linetype='dashed')+
      geom_line(data=ggreal,mapping=aes(x=tps,y=est),color='blue',alpha=1)+
      geom_hline( yintercept=0, color='grey',linetype='dashed')+
      labs(x='Exposure level',y='LACE estimates',
           subtitle =paste0("Q value/df/p-value: ",round(HeterQ[1],3),'/',HeterQ[2],'/',HeterQ[3],
                            '\n','Trend test value/df/p-value: ', round(TrendValues[1],3), '/',TrendValues[2],'/',round(TrendValues[3],6),
                            '\n','specific order term pvalues: ', round(metares$pval[2],3), ',', round(metares$pval[3],3)  ) )+
      theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      coord_cartesian(ylim =ylim_used ) #0.5这个距�?�有时候会很大
    #updates: now added the Q test and Trend test results
    SmoothRes$Q_test<-c(  'Q value'=round(HeterQ[1],3) , 'df'=HeterQ[2] , 'p-value'=HeterQ[3]  )
    SmoothRes$Trend_test<-c(  'Trend test value'=round(TrendValues[1],3) , 'df'=TrendValues[2] , 'p-value'=round(TrendValues[3],6) )
  }else{
    p<-ggplot(ggdata, aes(tps, est))+
      geom_point(data=ggdata,mapping=aes(x=tps,y=est),color='black',alpha=1)+
      geom_errorbar(data=ggdata,mapping=aes(x=tps,ymin=estlow,ymax=estup),color='black',width = 0.1)+
      geom_line(data=ggsmooth,mapping=aes(x=time_obs,y=smoothed_LACE_),color='black',alpha=1)+
      geom_line(data=ggsmooth,mapping=aes(x=time_obs,y=smooth_low),color='black',alpha=1,linetype='dashed')+
      geom_line(data=ggsmooth,mapping=aes(x=time_obs,y=smooth_up),color='black',alpha=1,linetype='dashed')+
      geom_line(data=ggreal,mapping=aes(x=tps,y=est),color='blue',alpha=1)+
      geom_hline( yintercept=0, color='grey',linetype='dashed')+
      labs(x='Exposure level',y='LACE estimates' )+
      theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      coord_cartesian(ylim = ylim_used )
  }

  if(Plot){print(p)}
  SmoothRes$p<-p








  ###SLR with my code
  # varSLRthetahat<-solve(  t( XX)%*%XX )*tau2
  # varM<- MM %*% varSLRthetahat %*%  t(  MM)
  # pw_std_error<-sqrt( diag(varM) )  #由于是Pointwise, 只取对�?�线元素�????
  # smoothed_LACE_<-eval.fd(    time_obs ,smoothed_LACE )
  # plot(time_obs,smoothed_LACE_,type='l')
  # lines(time_obs , smoothed_LACE_+1.96*pw_std_error, lty=2, lwd=1        )
  # lines( time_obs, smoothed_LACE_-1.96*pw_std_error, lty=2, lwd=1          )
  # #已经check；确实是一模一???

  ###For original effect shape
  if(is.na(baseline)){baseline<-mRES$mean[1]}
  baseline_used<-time_obs[sum(time_obs<=baseline)]
  cat(time_obs[1],time_obs[1000],baseline,'\n')
  cat('The actual basline value used:',baseline_used ,'\n')

  if(    (time_obs[1]-baseline)*(time_obs[1000]-baseline)>0     ){
    stop('please use the proper baseline value so it is inside the exposure range')
    }

  #DD: the integration matrix for the basisfunctions with given baseline exposrue level
  DD<-apply( MM, 2, cumsum )*(Rall[2]-Rall[1] )/(1000-1)  #\int_{Rall[1]}^x theta(s)ds where x=time_obs
  BaseInt<-matrix(rep( DD[sum(time_obs<=baseline),],1000  ), byrow=TRUE,nrow=1000 ) #int_{Rall[1]}^{baseline_used~+} theta(s)ds   #保证baseline_used那一行DD???0
  DD<-DD-BaseInt

  hest<-DD%*%thetahat  #h(x) pointwise est

  varD<- DD %*% varthetahat %*%  t(  DD)  #最后一项是random-effect
  pw_std_error<-sqrt( diag(varD) )  #由于是Pointwise, 只取对�?�线元素�????

  # plot(time_obs,hest,type='l', xlab='exposure level'  , ylab='h(x)')
  # lines(time_obs , hest+1.96*pw_std_error, lty=2, lwd=1        )
  # lines( time_obs, hest-1.96*pw_std_error, lty=2, lwd=1          )
  # abline(  h=0  ,lty=2, col ='blue' )
  # abline(  v=baseline_used  ,lty=2, col ='blue' )
  #checked, 两�?�DD的构造方法的可�?�化结果几乎一模一样！

  ggsmooth<-data.frame( time_obs=time_obs,
                        hest=hest,
                        smooth_low=hest-1.96*pw_std_error,
                        smooth_up=hest+1.96*pw_std_error)
  hp<-ggplot(ggdata, aes(tps, est))+
    geom_line(data=ggsmooth,mapping=aes(x=time_obs,y=hest),color='black',alpha=1)+
    geom_line(data=ggsmooth,mapping=aes(x=time_obs,y=smooth_low),color='black',alpha=1,linetype='dashed')+
    geom_line(data=ggsmooth,mapping=aes(x=time_obs,y=smooth_up),color='black',alpha=1,linetype='dashed')+
    geom_hline( yintercept=0, color='grey',linetype='dashed')+
    geom_vline( xintercept=baseline_used, color='grey',linetype='dashed')+
    labs(x='Exposure level',y='Effect')+
    theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    coord_cartesian(ylim = c( min(hest)-max(pw_std_error),max(hest)+max(pw_std_error)) )  #1 或�?0.5这些距�?�都太大�?
  if(Plot){print(hp)}
  SmoothRes$hp<-hp

  return(SmoothRes)

}


##example


# dat<-getDat( IVtype='cont', ZXmodel='A',XYmodel='2' ) #get a toy data
# rdat<-Stratify(dat)  #Do stratification on the data
# RES<-getSummaryInf( rdat,XYmodel='2')
#
# #fractional polynomial method with df=2
# smooth_res<-Smooth(RES$DRres,Norder=3,baseline=0)  #smoothing
# smooth_res$summary #B-spline basisfunctions parameters and s.e. p-values
# smooth_res$p #h'(x) function
# smooth_res$hp  #h(x) function
#
#
#
# #piecewise linear method
# cutting_values<-(RES$DRres$mean[-1] + head( RES$DRres$mean,-1) )/2
#
#
# smooth_res<-Smooth(RES$DRres,Norder=1,baseline=0,Knots=cutting_values)  #smoothing
# smooth_res$summary #B-spline basisfunctions parameters and s.e. p-values
# smooth_res$p #h'(x) function
# smooth_res$hp










