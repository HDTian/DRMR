




getSummaryInf<-function(rdat, #rdat #rdat<-Stratify(dat)   i.e. stratified individual dataset
                        covariate=FALSE, #whether or not adjust covariates?
                        target=FALSE, #if TRUE, will calculate the target causal effect for each stratum
                        XYmodel='1', #tell me what the true X-Y effect shape is? only applicable when target=TRUE
                        bxthre=1e-5, #bx threshold, under this value the bx is regarded as 0 or NA (so no MR ests)
                        getHeterQ=TRUE,  #get the heterogeneity Q statistics?
                        onlyDR=FALSE #only show the DR stratification results?
  ){

  #dat check
  if( is.null(rdat$Y) ){stop('No Y detected!: The column of the outcome should be named as Y') }

  #covariates position index
  cov_pos<-colnames(rdat)%in%paste0( 'C', 1:ncol(rdat) )
  if(  covariate&(sum(cov_pos)==0)  ){
    stop('if you wish to adjust the covariates, name the covariates as C1, C2, C3...')
  }

  if(covariate){
    cat(  'the covariates to adjuste are:', colnames(rdat)[cov_pos], '\n' )
  }


  Ns<-length(  levels(factor(rdat$DRstratum))  )

  stratify_on_X_or_not<-sum(rdat$M==rdat$X)==nrow(rdat)#judge whether the variable being stratified is the exposure X; TRUE menas yes the stratification is on the exposure


  RES<-list()


  if(target&( !stratify_on_X_or_not ) ){
    cat('The target effect is still provided but note that the stratified variable is not the exposure but your variable M', '\n',
        'You may wish to try target=FALSE or stratify on the exposure instead','\n')
  }


  if( target&stratify_on_X_or_not&(!XYmodel%in%c('1','2','3'))    ){
    if( exists('dif')){
      print('The (derivative) causal function is')
      print( dif )
      }else{ stop('Pleas define the effect function (dif) or use correct XYmodel index')   }
  }

  ##now assign the nonlinear_dif
  if( !XYmodel%in%c('1','2','3') ){   nonlinear_dif<-dif      }
  if( XYmodel=='1' ){    nonlinear_dif<-function(t){ sapply(   t , function(tt){ 0}) }     }
  if( XYmodel=='2' ){    nonlinear_dif<-function(t){0.2*t}   }
  if( XYmodel=='3' ){    nonlinear_dif<-function(t){-0.2*t*(t>0)}    }


  if( is.null(rdat$X_) ){ rdat$X_<-rdat$X  } #Target effect is the expected stratum-specific effect without any error (either confounding or coarsened error)





  if(!onlyDR){
    ##Residual stratification----------------------------------------------------------
    ##---------------------------------------------------------------------------------
    Size<-c()
    Minf<-c()  #M (either exposure or covariate) information
    MRfitting<-c()
    Target<-c()
    for(i in 1:Ns){
      selected_dat<-rdat[ rdat$Rstratum==i,  ]
      ##Stratum size
      Size<-c( Size , nrow(selected_dat))

      ##M information (min, mean, max)
      Minf<-rbind( Minf, as.numeric(summary(selected_dat$M   ))[c(1,2,4,5,6)]   )


      if( length(table(selected_dat$Z))==1  ){#如果Z都是单值的，谈何IV analysis?
        MRfitting<-rbind(MRfitting, c(NA,NA,NA,NA,NA,NA,NA)); Target<-c(Target ,   NA)
      }else{
        ##naive IVW MR est
        if(!covariate){ #no covariate adjustment
          fitGX<-lm(    selected_dat$X~  selected_dat$Z  );bGX<-summary(fitGX)$coef[-1,1] ;seGX<- summary(fitGX)$coef[-1,2] ; bx<-as.numeric( bGX  ); bxse<-as.numeric(  seGX)
          fitGY<-lm(    selected_dat$Y~  selected_dat$Z  );bGY<-summary(fitGY)$coef[-1,1] ;seGY<-summary(fitGY)$coef[-1,2] ; by<-as.numeric( bGY  ); byse<-as.numeric(  seGY)
          Fvalue<- summary(fitGX)$coef[2,3]^2    #IV strength - F value
        }else{ #with covariate adjustment
          fitGX<-lm(    selected_dat$X~  selected_dat$Z + as.matrix( selected_dat[,cov_pos] )  );bGX<-summary(fitGX)$coef[2,1] ;seGX<- summary(fitGX)$coef[2,2] ; bx<-as.numeric( bGX  ); bxse<-as.numeric(  seGX)
          fitGY<-lm(    selected_dat$Y~  selected_dat$Z + as.matrix( selected_dat[,cov_pos] ) );bGY<-summary(fitGY)$coef[2,1] ;seGY<-summary(fitGY)$coef[2,2] ; by<-as.numeric( bGY  ); byse<-as.numeric(  seGY)
          Fvalue<- summary(fitGX)$coef[2,3]^2    #IV strength - F value
        }

        if( length(bx)==0  ){#即，bx是numeric(0)；注意这并没有说明存在null G-X的情况，而是因为Z是单一值的
          MRfitting<-rbind(MRfitting, c(NA,NA,NA,by,byse,NA,NA))
        }else{
          if(abs(bx)<bxthre ){  #即，bx过小，可能产生inf项目，所以此时直接用NA代替
            MRfitting<-rbind(MRfitting, c(Fvalue,NA,NA,by,byse,NA,NA))
          }else{
            MRres<-mr_ivw(mr_input(bx, bxse, by, byse)) #checked! same as the 2SLS;    易理解，IVW唯一的变数在于random-effect;而single IV的情况random-effect肯定超过1；因为此时SEE=0呀
            MRfitting<-rbind(MRfitting, c(Fvalue,bx,bxse,by,byse,MRres@Estimate,MRres@StdError))
          }
        }



        ##Target causal effect (只有M==X时，算这个才有意义)
        #注意target是按照X求还是X_求？还是按照X_求吧
        if(target){
          #my definition
          At<-function(t){
            cov(     selected_dat$Z    , as.numeric(  selected_dat$X_ >= t  ) )/var( selected_dat$Z )
          }
          RRR<-c(min(selected_dat$X_  ),   max( selected_dat$X_    ) )  #RRR range for FDA
          time_obs<- seq( RRR[1],RRR[2], len=1000 )#length越大越精准
          At_<-function(t){sapply(t, function(t){ At(t)} ) } #At函数可积化
          Aobs<- At_(   time_obs )
          #B-spline smoothing
          Bs<-create.bspline.basis(  c(0,RRR[2]-RRR[1])  , nbasis=54 )  #50个interior knots (m默认均匀分布); order=4 by default
          Par <- fdPar( fdobj=Bs, Lfdobj=NULL, lambda=0 )
          smoothed_fd <- smooth.basis(  time_obs-RRR[1]  , Aobs  , Par  )$fd
          smoothed_At<-eval.fd(smoothed_fd, seq(0,RRR[2]-RRR[1] , len=1000 )   )
          #最终就是获得At这条曲线在time_obs点上的值，方便积分而已
          #其实使用FDA没有太大必要，

          hdx<-nonlinear_dif( time_obs      ) #h'(x)在time_obs上的值
          times<-smoothed_At*hdx  #At*hdx 在time_obs上的值
          intres<-int.simpson2( time_obs ,times, equi = TRUE, method = NULL)  #integration result
          #intres<-int.simpson2( seq(0,RRR[2]-RRR[1] , len=1000 ) ,times, equi = TRUE, method = NULL)  #integration result
          #intres<-int.simpson2( time_obs ,Aobs*hdx , equi = TRUE, method = NULL)

          fitGX_<-lm(    selected_dat$X_~  selected_dat$Z  );bx_<-as.numeric(summary(fitGX_)$coef[-1,1])
          targetvalue<-intres/bx_  #target effect  (除以bx是因为之前的At分母并不是cov(Z,X)) #bx_是依据X_算的；即，可能coarsened
          Target<-c(Target ,   targetvalue)
        }else{Target<-c(Target ,   NA)}

      }




    }

    res<-cbind( Size,Minf, MRfitting , Target)
    res<-as.data.frame(res)
    colnames(res)<-c('size', 'min','1q','mean', '3q', 'max','Fvalue','bx','bxse','by','byse','est','se','target')
    rownames(res)<-paste0('Stratum ', 1:Ns  )

    RES$Rres<-res

    if(getHeterQ){
      Sr<-mr_ivw(mr_input(res$bx, res$bxse ,res$by, res$byse))#naive IVW under null (no effect heterogeneity)
      #mr_ivw是免疫NA的，即，即使把NA移除掉，naive IVW estimate依旧是一样的
      delta<-Sr@Estimate #一维的
      #Cochran's Q statistic
      ssigma_square <- res$byse^2 + delta^2*res$bxse^2
      QQ<-sum( ( res$by - delta*res$bx     )^2 /ssigma_square ,na.rm=TRUE      )
      #critical value:
      Qdf<-sum(!is.na( ( res$by - delta*res$bx     )^2 /ssigma_square  ))    -1#normally Ns-1 without NA stratum
      pvalue<-1-pchisq(QQ, Qdf ) #不会考虑NA的stratum
      HeterQ<-c(QQ, Qdf, pvalue  ); names(HeterQ)<-c('Q statistic' ,'df' ,'p-value'  )
      RES$RHeterQ<-HeterQ
    }
  }


  ##Doubly-ranked stratification----------------------------------------------------------
  ##---------------------------------------------------------------------------------
  Size<-c()
  Minf<-c()  #M (either exposure or covariate) information
  MRfitting<-c()
  Target<-c()
  for(i in 1:Ns){
    selected_dat<-rdat[ rdat$DRstratum==i,  ]
    ##Stratum size
    Size<-c( Size , nrow(selected_dat))

    ##M information (min, mean, max)
    Minf<-rbind( Minf, as.numeric(summary(selected_dat$M   ))[c(1,2,4,5,6)]   )
    if( length(table(selected_dat$Z))==1  ){#如果Z都是单值的，谈何IV analysis?
      MRfitting<-rbind(MRfitting, c(NA,NA,NA,NA,NA,NA,NA)); Target<-c(Target ,   NA)
    }else{
      ##naive IVW MR est
      if(!covariate){ #no covariate adjustment
        fitGX<-lm(    selected_dat$X~  selected_dat$Z  );bGX<-summary(fitGX)$coef[-1,1] ;seGX<- summary(fitGX)$coef[-1,2] ; bx<-as.numeric( bGX  ); bxse<-as.numeric(  seGX)
        fitGY<-lm(    selected_dat$Y~  selected_dat$Z  );bGY<-summary(fitGY)$coef[-1,1] ;seGY<-summary(fitGY)$coef[-1,2] ; by<-as.numeric( bGY  ); byse<-as.numeric(  seGY)
        Fvalue<- summary(fitGX)$coef[2,3]^2    #IV strength - F value
      }else{ #with covariate adjustment
        fitGX<-lm(    selected_dat$X~  selected_dat$Z + as.matrix( selected_dat[,cov_pos] )  );bGX<-summary(fitGX)$coef[2,1] ;seGX<- summary(fitGX)$coef[2,2] ; bx<-as.numeric( bGX  ); bxse<-as.numeric(  seGX)
        fitGY<-lm(    selected_dat$Y~  selected_dat$Z + as.matrix( selected_dat[,cov_pos] ) );bGY<-summary(fitGY)$coef[2,1] ;seGY<-summary(fitGY)$coef[2,2] ; by<-as.numeric( bGY  ); byse<-as.numeric(  seGY)
        Fvalue<- summary(fitGX)$coef[2,3]^2    #IV strength - F value
      }
      if( length(bx)==0  ){#即，bx是numeric(0)；说明存在null G-X的情况
        MRfitting<-rbind(MRfitting, c(NA,NA,NA,by,byse,NA,NA))
      }else{
        if(abs(bx)<bxthre ){  #即，bx过小，可能产生inf项目，所以此时直接用NA代替
          MRfitting<-rbind(MRfitting, c(Fvalue,NA,NA,by,byse,NA,NA))
        }else{
          MRres<-mr_ivw(mr_input(bx, bxse, by, byse))
          MRfitting<-rbind(MRfitting, c(Fvalue,bx,bxse,by,byse,MRres@Estimate,MRres@StdError))
        }
      }

      ##Target causal effect (只有M==X时，算这个才有意义)
      #注意target是按照X求还是X_求？还是按照X_求吧
      if(target){
        #my definition
        At<-function(t){
          cov(     selected_dat$Z    , as.numeric(  selected_dat$X_ >= t  ) )/var( selected_dat$Z )
        }
        RRR<-c(min(selected_dat$X_  ),   max( selected_dat$X_    ) )  #RRR range for FDA
        time_obs<- seq( RRR[1],RRR[2], len=1000 )#length越大越精准
        At_<-function(t){sapply(t, function(t){ At(t)} ) } #At函数可积化
        Aobs<- At_(   time_obs )
        #B-spline smoothing
        Bs<-create.bspline.basis(  c(0,RRR[2]-RRR[1])  , nbasis=54 )  #50个interior knots (m默认均匀分布); order=4 by default
        Par <- fdPar( fdobj=Bs, Lfdobj=NULL, lambda=0 )
        smoothed_fd <- smooth.basis(  time_obs-RRR[1]  , Aobs  , Par  )$fd
        smoothed_At<-eval.fd(smoothed_fd, seq(0,RRR[2]-RRR[1] , len=1000 )   )
        #最终就是获得At这条曲线在time_obs点上的值，方便积分而已
        #其实使用FDA没有太大必要，

        hdx<-nonlinear_dif( time_obs      ) #h'(x)在time_obs上的值
        times<-smoothed_At*hdx  #At*hdx 在time_obs上的值
        intres<-int.simpson2( time_obs ,times, equi = TRUE, method = NULL)  #integration result
        #intres<-int.simpson2( seq(0,RRR[2]-RRR[1] , len=1000 ) ,times, equi = TRUE, method = NULL)  #integration result
        #intres<-int.simpson2( time_obs ,Aobs*hdx , equi = TRUE, method = NULL)

        fitGX_<-lm(    selected_dat$X_~  selected_dat$Z  );bx_<-as.numeric(summary(fitGX_)$coef[-1,1])
        targetvalue<-intres/bx_  #target effect  (除以bx是因为之前的At分母并不是cov(Z,X)) #bx_是依据X_算的；即，可能coarsened
        Target<-c(Target ,   targetvalue)
      }else{Target<-c(Target ,   NA)}
    }


  }

  res<-cbind( Size,Minf, MRfitting , Target)
  res<-as.data.frame(res)
  colnames(res)<-c( 'size', 'min', '1q','mean','3q','max','Fvalue','bx','bxse','by','byse','est','se','target')
  rownames(res)<-paste0('Stratum ', 1:Ns  )

  RES$DRres<-res

  if(getHeterQ){
    Sr<-mr_ivw(mr_input(res$bx, res$bxse ,res$by, res$byse))#naive IVW under null (no effect heterogeneity)
    #mr_ivw是免疫NA的，即，即使把NA移除掉，naive IVW estimate依旧是一样的
    delta<-Sr@Estimate #一维的
    #Cochran's Q statistic
    ssigma_square <- res$byse^2 + delta^2*res$bxse^2
    QQ<-sum( ( res$by - delta*res$bx     )^2 /ssigma_square ,na.rm=TRUE      )
    #critical value:
    Qdf<-sum(!is.na( ( res$by - delta*res$bx     )^2 /ssigma_square  ))    -1#normally Ns-1 without NA stratum
    pvalue<-1-pchisq(QQ, Qdf ) #不会考虑NA的stratum
    HeterQ<-c(QQ, Qdf, pvalue  ); names(HeterQ)<-c('Q statistic' ,'df' ,'p-value'  )
    RES$DRHeterQ<-HeterQ
  }


  return(RES)
}

#Warning message:
#  In summary.lm(fitGX) : essentially perfect fit: summary may be unreliable
#Z肯定是多值的(否则不会跑到这一条代码)，这其实代表着X是一个值(往往是coarsened导致的)；此时bx为0，
#代码已经设计了此时会把bx=NA by不用管 MRest=NA MRse=NA; 不用管target，因为target是针对X_，而X_是连续的


###examples
# getSummaryInf(rdat, target=FALSE, bxthre=1e-5)
#
#
#
# dat<-getDat( IVtype='cont', ZXmodel='C',XYmodel='1' ) #get a toy data
# rdat<-Stratify(dat)  #Do stratification on the data
# RES<-getSummaryInf( rdat, target=TRUE, bxthre=1e-5, XYmodel='1',getHeterQ=TRUE)

