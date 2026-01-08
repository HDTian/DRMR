

### Rubin rule resampling 




DR_RubinRule<-function(data,  # the complete one-sample individual-level data
                       times=100,
                       Ns=10, # DR stratification: the number of final strata
                       SoP=NA, # DR stratification: the size of pre-stratum
                       onExposure= TRUE ,  # DR stratification: stratify on X or M?
                       family_used='gaussian' , # stratum-specific MR estimation: response model
                       covariate= FALSE  # stratum-specific MR estimation: adjust covariate?
                       ){
  
  RES<-list()
  
  EST_ZX<-c();EST_ZY<-c() #dim(EST) 10*mtimes
  SE_ZX<-c();SE_ZY<-c() #dim(SE) 10*mtimes
  MEANX<-c()  #dim(MEANX) 10*mtimes
  mtimes<-times # the rubin rule times
  for(m in 1:mtimes){
    dat_used<-dat[-sample(1:nrow(dat), 10)  ,  ] # random remove 10 individuals
    
    # stratification via DRMR and get the stratum-specific (est and s.e.) and Xmean
    rdat<-Stratify(dat_used,onExposure = onExposure, Ns=Ns, SoP=SoP )#onExposure = FALSE: stratify on M (the G*E-adjusted exposure)
    
    # get the stratum-specific statistics (only cat the black message at m==1)
    if (m == 1) {
      SRES<-getSummaryInf(rdat, family_used=family_used,covariate=covariate,onlyDR = TRUE) #Stratification RESults
    } else { 
      tmp <- NULL
      invisible(capture.output(
        tmp <- getSummaryInf( rdat, family_used = family_used,covariate   = covariate,onlyDR      = TRUE),
        file = NULL
      ))
      SRES <- tmp
     }  
    
    cat(  paste0(m,'-') )
    
    EST_ZX<-cbind( EST_ZX, SRES$DRres$bx );EST_ZY<-cbind( EST_ZY, SRES$DRres$by )
    SE_ZX<-cbind(SE_ZX, SRES$DRres$bxse  );SE_ZY<-cbind(SE_ZY, SRES$DRres$byse  )
    MEANX<-cbind( MEANX , SRES$DRres$mean)  # for trend test (if applicable)
  }
  ###Rubin's Rule
  #for ZX association ---------
  Qbar_ZX<-apply( EST_ZX , 1, mean )
  TT1_ZX<-apply( SE_ZX^2 , 1, mean  ) #i.e. Ubar  #a vector
  TT2_ZX<-(1+ 1/mtimes)*(  apply(  ( EST_ZX  -   Qbar_ZX  )^2 ,1 , sum    )/(mtimes-1)  ) #a vector 
  TT_ZX<-TT1_ZX+TT2_ZX#a vector
  Fdf_ZX<-(mtimes -1 )*(1+ TT1_ZX/TT2_ZX)^2   #the second df for F statistic (越大，使得sqrt(F_1,df))越接近Z statistic
  #qf(0.95, 1,Fdf_ZX )#still a vector
  
  #for ZY association -------------
  Qbar_ZY<-apply( EST_ZY , 1, mean )
  TT1_ZY<-apply( SE_ZY^2 , 1, mean  ) #i.e. Ubar  #a vector
  TT2_ZY<-(1+ 1/mtimes)*(  apply(  ( EST_ZY  -   Qbar_ZY  )^2 ,1 , sum    )/(mtimes-1)  ) #a vector 
  TT_ZY<-TT1_ZY+TT2_ZY#a vector
  Fdf_ZY<-(mtimes -1 )*(1+ TT1_ZY/TT2_ZY)^2   #the second df for F statistic (越大，使得sqrt(F_1,df))越接近Z statistic
  #qf(0.95, 1,Fdf_ZY )#still a vector
  
  #stratum specific MR est: simply ignore the uncertainty of the estimated ZX association
  MRest<-Qbar_ZY/Qbar_ZX
  ##MR estimand CI_low and CI_up based on F test statistic with ignoreing the estimated ZX association
  ggdata<-data.frame( X = apply( MEANX ,1 ,mean)  , 
                      Est = MRest    , 
                      CI_low = (  Qbar_ZY- sqrt( qf(0.95, 1,Fdf_ZY ) )*sqrt( TT_ZY )   )/Qbar_ZX, 
                      CI_up =  (   Qbar_ZY+ sqrt( qf(0.95, 1,Fdf_ZY ) )*sqrt( TT_ZY )  )/Qbar_ZX 
  ) # sqrt( qf(0.95, 1,Fdf_ZY ) )  ~= 1.96 if the Fdf_ZY is large
  
  ##also store some other metrics (possible for some trend test)
  ggdata$ZXEst<-Qbar_ZX;ggdata$ZYEst<-Qbar_ZY
  ggdata$TT_ZX<-TT_ZX;ggdata$TT_ZY<-TT_ZY
  ggdata$Fdf_ZX<-Fdf_ZX;ggdata$Fdf_ZY<-Fdf_ZY
  
  RES$ggdata<-ggdata  # the stratum specific data: mean exposure, LACE, 95% CI, stratum-specific ZX ZY association
  
  ##基于Multiple imputation后的stratum specific结果算Q test result ------------------------------------
  #here each stratum-specific ZX ZY associations are assumed to be Gaussian distributed
  #that is, Qbar_ZX ~ N(  Q_ZX , sqrt(TT_ZX) ^2 ) where Q_ZX represents the LACE estimand
  Bx<-Qbar_ZX;Bxse<-sqrt(TT_ZX); By<-Qbar_ZY;Byse<-sqrt(TT_ZY)
  #先算effect via naive IVW under null (no effect difference) - i.e. ignore the ZX estimate uncertainty
  Sr<-mr_ivw(mr_input(Bx, (1:length(Bx)) , By, Byse))  # bxse is not important
  delta<-Sr@Estimate #一维的
  #Cochran's Q statistic -------------------
  ssigma_square <- Byse^2 + delta^2*Bxse^2
  QQ<-sum( ( By - delta*Bx     )^2 /ssigma_square      )
  pvalue<-1-pchisq(QQ, length(  By) - 1)
  RES$Qres<- c(QQ, length(  By) - 1 ,pvalue)  # Q statistic value; Q test pvalue
  names(  RES$Qres ) <- c('Qstatistic', 'df','pvalue')
  
  
  ##基于Multiple imputation后的stratum specific结果算trend test result ------------------------------------
  Est<-Qbar_ZY/Qbar_ZX  # stratum-specific MR estimate
  se_Est <- sqrt( TT_ZY/Qbar_ZX^2 + (Qbar_ZY^2)*TT_ZX/Qbar_ZX^4 )  # Delta method with 2rd order 
  X = apply( MEANX ,1 ,mean) 
  metares <- rma(yi = Est, sei = se_Est, mods = X)  # default is method = "REML" - ie random effect model
  RES$Trendres<-c( metares$QM , metares$QMdf[1] , metares$QMp    ) # coeffcient test value; df , and pvalue
  names(  RES$Trendres ) <- c('Trend_statistic', 'df','pvalue')
  
  # LACE plot (use ggdata) ------------------
  p<-ggplot(ggdata, aes(x = X, y = Est)) + geom_point(size = 2) +
    geom_errorbar( aes(ymin = CI_low , ymax = CI_up), width = 0 ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs( x = "Mean exposure within stratum", y = "Local causal effect estimate",) + theme_classic()
  RES$p <- p
  
  
  return(RES)
  
}


### output:
# RES$ggdata : stratum-specific statistics (after Rubin Rule stabilization)
# RES$Qres : Cochran's Q heterogeneity test
# RES$Trendres : meta regression trend test test wrt the exposure mean (ie test whether the exposure mean has no association with the stratum-specific MR estimates)
# RES$p: LACE plot





