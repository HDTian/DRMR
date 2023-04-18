




getSummaryInf<-function(rdat, #rdat #rdat<-Stratify(dat)   i.e. stratified individual dataset
                        target=TRUE, #if TRUE, will calculate the target causal effect for each stratum
                        XYmodel='1', #tell me what the true X-Y effect shape is? only applicable when target=TRUE #XYmodel='0' : use self-defined effect shape, in global environment must define:  dif()
                        bxthre=1e-5, #bx threshold, under this value the bx is regarded as 0 or NA (so no MR ests)
                        getHeterQ=TRUE,  #get the heterogeneity Q statistics?
                        onlyDR=FALSE #only show the DR stratification results?
  ){

  #dat check
  if( is.null(rdat$Y) ){stop('No Y detected!: The column of the outcome should be named as Y') }


  Ns<-length(  levels(factor(rdat$DRstratum))  )
  stratify_on_X_or_not<-sum(rdat$M==rdat$X)==nrow(rdat)


  RES<-list()


  if(target&( !stratify_on_X_or_not ) ){
    stop('The target is not applicable as the stratified variable is not the exposure itself;
        Try target=FALSE or stratify on the exposure instead')
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


      if( length(table(selected_dat$Z))==1  ){#���Z���ǵ�ֵ�ģ�̸��IV analysis?
        MRfitting<-rbind(MRfitting, c(NA,NA,NA,NA,NA,NA)); Target<-c(Target ,   NA)
      }else{
        ##naive IVW MR est
        fitGX<-lm(    selected_dat$X~  selected_dat$Z  );bGX<-summary(fitGX)$coef[-1,1] ;seGX<- summary(fitGX)$coef[-1,2] ; bx<-as.numeric( bGX  ); bxse<-as.numeric(  seGX)
        fitGY<-lm(    selected_dat$Y~  selected_dat$Z  );bGY<-summary(fitGY)$coef[-1,1] ;seGY<-summary(fitGY)$coef[-1,2] ; by<-as.numeric( bGY  ); byse<-as.numeric(  seGY)
        if( length(bx)==0  ){#����bx��numeric(0)���Ⲣû��˵������null G-X�������������ΪZ�ǵ�һֵ��
          MRfitting<-rbind(MRfitting, c(NA,NA,by,byse,NA,NA))
        }else{
          if(abs(bx)<bxthre ){  #����bx��С�����ܲ���inf��Ŀ�����Դ�ʱֱ����NA����
            MRfitting<-rbind(MRfitting, c(NA,NA,by,byse,NA,NA))
          }else{
            MRres<-mr_ivw(mr_input(bx, bxse, by, byse))
            MRfitting<-rbind(MRfitting, c(bx,bxse,by,byse,MRres@Estimate,MRres@StdError))
          }
        }



        ##Target causal effect (ֻ��M==Xʱ���������������)
        #ע��target�ǰ���X����X_�󣿻��ǰ���X_���
        if(target){
          #my definition
          At<-function(t){
            cov(     selected_dat$Z    , as.numeric(  selected_dat$X_ >= t  ) )/var( selected_dat$Z )
          }
          RRR<-c(min(selected_dat$X_  ),   max( selected_dat$X_    ) )  #RRR range for FDA
          time_obs<- seq( RRR[1],RRR[2], len=1000 )#lengthԽ��Խ��׼
          At_<-function(t){sapply(t, function(t){ At(t)} ) } #At�����ɻ���
          Aobs<- At_(   time_obs )
          #B-spline smoothing
          Bs<-create.bspline.basis(  c(0,RRR[2]-RRR[1])  , nbasis=54 )  #50��interior knots (mĬ�Ͼ��ȷֲ�); order=4 by default
          Par <- fdPar( fdobj=Bs, Lfdobj=NULL, lambda=0 )
          smoothed_fd <- smooth.basis(  time_obs-RRR[1]  , Aobs  , Par  )$fd
          smoothed_At<-eval.fd(smoothed_fd, seq(0,RRR[2]-RRR[1] , len=1000 )   )
          #���վ��ǻ��At����������time_obs���ϵ�ֵ��������ֶ���
          #��ʵʹ��FDAû��̫���Ҫ��

          hdx<-nonlinear_dif( time_obs      ) #h'(x)��time_obs�ϵ�ֵ
          times<-smoothed_At*hdx  #At*hdx ��time_obs�ϵ�ֵ
          intres<-int.simpson2( time_obs ,times, equi = TRUE, method = NULL)  #integration result
          #intres<-int.simpson2( seq(0,RRR[2]-RRR[1] , len=1000 ) ,times, equi = TRUE, method = NULL)  #integration result
          #intres<-int.simpson2( time_obs ,Aobs*hdx , equi = TRUE, method = NULL)

          fitGX_<-lm(    selected_dat$X_~  selected_dat$Z  );bx_<-as.numeric(summary(fitGX_)$coef[-1,1])
          targetvalue<-intres/bx_  #target effect  (����bx����Ϊ֮ǰ��At��ĸ������cov(Z,X)) #bx_������X_��ģ���������coarsened
          Target<-c(Target ,   targetvalue)
        }else{Target<-c(Target ,   NA)}

      }




    }

    res<-cbind( Size,Minf, MRfitting , Target)
    res<-as.data.frame(res)
    colnames(res)<-c('size', 'min','1q','mean', '3q', 'max','bx','bxse','by','byse','est','se','target')
    rownames(res)<-paste0('Stratum ', 1:Ns  )

    RES$Rres<-res

    if(getHeterQ){
      Sr<-mr_ivw(mr_input(res$bx, res$bxse ,res$by, res$byse))#naive IVW under null (no effect heterogeneity)
      #mr_ivw������NA�ģ�������ʹ��NA�Ƴ�����naive IVW estimate������һ����
      delta<-Sr@Estimate #һά��
      #Cochran's Q statistic
      ssigma_square <- res$byse^2 + delta^2*res$bxse^2
      QQ<-sum( ( res$by - delta*res$bx     )^2 /ssigma_square ,na.rm=TRUE      )
      #critical value:
      Qdf<-sum(!is.na( ( res$by - delta*res$bx     )^2 /ssigma_square  ))    -1#normally Ns-1 without NA stratum
      pvalue<-1-pchisq(QQ, Qdf ) #���ῼ��NA��stratum
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
    if( length(table(selected_dat$Z))==1  ){#���Z���ǵ�ֵ�ģ�̸��IV analysis?
      MRfitting<-rbind(MRfitting, c(NA,NA,NA,NA,NA,NA)); Target<-c(Target ,   NA)
    }else{
      ##naive IVW MR est
      fitGX<-lm(    selected_dat$X~  selected_dat$Z  );bGX<-summary(fitGX)$coef[-1,1] ;seGX<- summary(fitGX)$coef[-1,2] ; bx<-as.numeric( bGX  ); bxse<-as.numeric(  seGX)
      fitGY<-lm(    selected_dat$Y~  selected_dat$Z  );bGY<-summary(fitGY)$coef[-1,1] ;seGY<-summary(fitGY)$coef[-1,2] ; by<-as.numeric( bGY  ); byse<-as.numeric(  seGY)
      if( length(bx)==0  ){#����bx��numeric(0)��˵������null G-X�����
        MRfitting<-rbind(MRfitting, c(NA,NA,by,byse,NA,NA))
      }else{
        if(abs(bx)<bxthre ){  #����bx��С�����ܲ���inf��Ŀ�����Դ�ʱֱ����NA����
          MRfitting<-rbind(MRfitting, c(NA,NA,by,byse,NA,NA))
        }else{
          MRres<-mr_ivw(mr_input(bx, bxse, by, byse))
          MRfitting<-rbind(MRfitting, c(bx,bxse,by,byse,MRres@Estimate,MRres@StdError))
        }
      }

      ##Target causal effect (ֻ��M==Xʱ���������������)
      #ע��target�ǰ���X����X_�󣿻��ǰ���X_���
      if(target){
        #my definition
        At<-function(t){
          cov(     selected_dat$Z    , as.numeric(  selected_dat$X_ >= t  ) )/var( selected_dat$Z )
        }
        RRR<-c(min(selected_dat$X_  ),   max( selected_dat$X_    ) )  #RRR range for FDA
        time_obs<- seq( RRR[1],RRR[2], len=1000 )#lengthԽ��Խ��׼
        At_<-function(t){sapply(t, function(t){ At(t)} ) } #At�����ɻ���
        Aobs<- At_(   time_obs )
        #B-spline smoothing
        Bs<-create.bspline.basis(  c(0,RRR[2]-RRR[1])  , nbasis=54 )  #50��interior knots (mĬ�Ͼ��ȷֲ�); order=4 by default
        Par <- fdPar( fdobj=Bs, Lfdobj=NULL, lambda=0 )
        smoothed_fd <- smooth.basis(  time_obs-RRR[1]  , Aobs  , Par  )$fd
        smoothed_At<-eval.fd(smoothed_fd, seq(0,RRR[2]-RRR[1] , len=1000 )   )
        #���վ��ǻ��At����������time_obs���ϵ�ֵ��������ֶ���
        #��ʵʹ��FDAû��̫���Ҫ��

        hdx<-nonlinear_dif( time_obs      ) #h'(x)��time_obs�ϵ�ֵ
        times<-smoothed_At*hdx  #At*hdx ��time_obs�ϵ�ֵ
        intres<-int.simpson2( time_obs ,times, equi = TRUE, method = NULL)  #integration result
        #intres<-int.simpson2( seq(0,RRR[2]-RRR[1] , len=1000 ) ,times, equi = TRUE, method = NULL)  #integration result
        #intres<-int.simpson2( time_obs ,Aobs*hdx , equi = TRUE, method = NULL)

        fitGX_<-lm(    selected_dat$X_~  selected_dat$Z  );bx_<-as.numeric(summary(fitGX_)$coef[-1,1])
        targetvalue<-intres/bx_  #target effect  (����bx����Ϊ֮ǰ��At��ĸ������cov(Z,X)) #bx_������X_��ģ���������coarsened
        Target<-c(Target ,   targetvalue)
      }else{Target<-c(Target ,   NA)}
    }


  }

  res<-cbind( Size,Minf, MRfitting , Target)
  res<-as.data.frame(res)
  colnames(res)<-c( 'size', 'min', '1q','mean','3q','max','bx','bxse','by','byse','est','se','target')
  rownames(res)<-paste0('Stratum ', 1:Ns  )

  RES$DRres<-res

  if(getHeterQ){
    Sr<-mr_ivw(mr_input(res$bx, res$bxse ,res$by, res$byse))#naive IVW under null (no effect heterogeneity)
    #mr_ivw������NA�ģ�������ʹ��NA�Ƴ�����naive IVW estimate������һ����
    delta<-Sr@Estimate #һά��
    #Cochran's Q statistic
    ssigma_square <- res$byse^2 + delta^2*res$bxse^2
    QQ<-sum( ( res$by - delta*res$bx     )^2 /ssigma_square ,na.rm=TRUE      )
    #critical value:
    Qdf<-sum(!is.na( ( res$by - delta*res$bx     )^2 /ssigma_square  ))    -1#normally Ns-1 without NA stratum
    pvalue<-1-pchisq(QQ, Qdf ) #���ῼ��NA��stratum
    HeterQ<-c(QQ, Qdf, pvalue  ); names(HeterQ)<-c('Q statistic' ,'df' ,'p-value'  )
    RES$DRHeterQ<-HeterQ
  }


  return(RES)
}

#Warning message:
#  In summary.lm(fitGX) : essentially perfect fit: summary may be unreliable
#Z�϶��Ƕ�ֵ��(���򲻻��ܵ���һ������)������ʵ������X��һ��ֵ(������coarsened���µ�)����ʱbxΪ0��
#�����Ѿ�����˴�ʱ���bx=NA by���ù� MRest=NA MRse=NA; ���ù�target����Ϊtarget�����X_����X_��������


###examples
# getSummaryInf(rdat, target=FALSE, bxthre=1e-5)
#
#
#
# dat<-getDat( IVtype='cont', ZXmodel='C',XYmodel='1' ) #get a toy data
# rdat<-Stratify(dat)  #Do stratification on the data
# RES<-getSummaryInf( rdat, target=TRUE, bxthre=1e-5, XYmodel='1',getHeterQ=TRUE)
