



#######################################################################################
#This is a guide script to reproduce the simulation results of the doubly-ranked paper#
#######################################################################################



library(tidyverse)
library(MendelianRandomization)
library(facetscales) #remotes::install_github("zeehio/facetscales@archived")
library(fda)
library(fda.usc)
library(xtable)



#devtools::install_github("HDTian/DRMR")
#library(DRMR)

#path<-'...' #define your desired path

##continuous or binary IV---------------------------------------------------------------
##--------------------------------------------------------------------------------------
theIVtype<-'cont'
#theIVtype<-'bi'
Nsim<-1000
for(theZXmodel in c('A','B','C','D')){
  print(theZXmodel)
  REst<-c() ; DREst<-c()
  RTarget<-c() ;  DRTarget<-c()
  RMSE<-c() ; DRMSE<-c()
  Rcov<-c()  ; DRcov<-c()
  for(j in c('1','2','3')){
    print(j)
    for(i in 1:Nsim){
      dat<-getDat( IVtype=theIVtype, ZXmodel=theZXmodel,XYmodel=j ) ;rdat<-Stratify(dat)
      RES<-getSummaryInf( rdat,XYmodel=j)
      Rres<-RES$Rres ; DRres<-RES$DRres
      REst<-cbind(REst, Rres$est); DREst<-cbind(DREst, DRres$est)
      RTarget<-cbind( RTarget, Rres$target ) ;DRTarget<-cbind( DRTarget, DRres$target )
      RMSE<-cbind(RMSE, (Rres$est-Rres$target)^2)  ; DRMSE<-cbind(DRMSE, (DRres$est-DRres$target)^2)
      Rcov<-cbind(Rcov,  abs(  Rres$est-Rres$target )<1.96*Rres$se);DRcov<-cbind(DRcov,  abs(  DRres$est-DRres$target )<1.96*DRres$se)
    }
  }

  ##table
  MCtable<-cbind( c( mean(as.vector(RMSE[,1:Nsim]), na.rm=TRUE  ), mean(as.vector(RMSE[,Nsim+1:Nsim]) , na.rm=TRUE ) ,mean(as.vector(RMSE[,2*Nsim+1:Nsim]) , na.rm=TRUE )   ),
                  c( mean(as.vector(Rcov[,1:Nsim]) , na.rm=TRUE ), mean(as.vector(Rcov[,Nsim+1:Nsim]) , na.rm=TRUE ) ,mean(as.vector(Rcov[,2*Nsim+1:Nsim]) , na.rm=TRUE )   ),
                  c( mean(as.vector(DRMSE[,1:Nsim]) , na.rm=TRUE ), mean(as.vector(DRMSE[,Nsim+1:Nsim]) , na.rm=TRUE ) ,mean(as.vector(DRMSE[,2*Nsim+1:Nsim]) , na.rm=TRUE )   ),
                  c( mean(as.vector(DRcov[,1:Nsim]) , na.rm=TRUE ), mean(as.vector(DRcov[,Nsim+1:Nsim]) , na.rm=TRUE ) ,mean(as.vector(DRcov[,2*Nsim+1:Nsim]), na.rm=TRUE  )   )
  )  #MC: MSE and Coverage

  ##visualization
  ggdata<-data.frame(Stratum=1:10 ,
                     ests=c(as.vector(REst),as.vector(DREst)) ,
                     Stratification =c(  rep('Residual method', 3*10*Nsim  ),rep('Doubly-ranked method', 3*10*Nsim )),
                     XYmodel= rep(  paste0(theZXmodel,1:3) , each=10*Nsim   )
  )
  Mean<-function(x){mean(x,na.rm=TRUE)  }
  RTargetmean<-c( apply(RTarget[,1:Nsim] , 1, Mean), apply(RTarget[,Nsim+1:Nsim] , 1, Mean) ,apply(RTarget[,2*Nsim+1:Nsim] , 1, Mean)   )
  DRTargetmean<-c( apply(DRTarget[,1:Nsim] , 1, Mean), apply(DRTarget[,Nsim+1:Nsim] , 1, Mean) ,apply(DRTarget[,2*Nsim+1:Nsim] , 1, Mean)   )
  ggtrue<- data.frame(   Stratum=1:10  ,
                         tragets=c(  RTargetmean,DRTargetmean ) ,
                         Stratification =c(  rep('Residual method', 3*10 ),rep('Doubly-ranked method', 3*10 )),
                         XYmodel= rep( paste0(theZXmodel,1:3) , each=10  )
  )

  ###storage
  write.csv(MCtable  , paste0(path,theIVtype,'\\',"MCtable_",theZXmodel  , '.csv'), row.names=F)
  write.csv(ggdata  , paste0(path,theIVtype,'\\',"ggdata_",theZXmodel  , '.csv'), row.names=F)
  write.csv(ggtrue  , paste0(path,theIVtype,'\\',"ggtrue_",theZXmodel  , '.csv'), row.names=F)
}

##multiple instruments---------------------------------------------------------
#------------------------------------------------------------------------------
theIVtype<-'high-dim'
Nsim<-1000
for(theZXmodel in c('E','F','G','H')){
  print(theZXmodel)
  REst<-c() ; DREst<-c()
  RTarget<-c() ;  DRTarget<-c()
  RMSE<-c() ; DRMSE<-c()
  Rcov<-c()  ; DRcov<-c()
  for(j in c('1','2','3')){
    print(j)
    for(i in 1:Nsim){
      dat<-getDat( IVtype=theIVtype, ZXmodel=theZXmodel,XYmodel=j ) ;rdat<-Stratify(dat)
      RES<-getSummaryInf( rdat,XYmodel=j)
      Rres<-RES$Rres ; DRres<-RES$DRres
      REst<-cbind(REst, Rres$est); DREst<-cbind(DREst, DRres$est)
      RTarget<-cbind( RTarget, Rres$target ) ;DRTarget<-cbind( DRTarget, DRres$target )
      RMSE<-cbind(RMSE, (Rres$est-Rres$target)^2)  ; DRMSE<-cbind(DRMSE, (DRres$est-DRres$target)^2)
      Rcov<-cbind(Rcov,  abs(  Rres$est-Rres$target )<1.96*Rres$se);DRcov<-cbind(DRcov,  abs(  DRres$est-DRres$target )<1.96*DRres$se)
    }
  }

  ##table
  MCtable<-cbind( c( mean(as.vector(RMSE[,1:Nsim]), na.rm=TRUE  ), mean(as.vector(RMSE[,Nsim+1:Nsim]) , na.rm=TRUE ) ,mean(as.vector(RMSE[,2*Nsim+1:Nsim]) , na.rm=TRUE )   ),
                  c( mean(as.vector(Rcov[,1:Nsim]) , na.rm=TRUE ), mean(as.vector(Rcov[,Nsim+1:Nsim]) , na.rm=TRUE ) ,mean(as.vector(Rcov[,2*Nsim+1:Nsim]) , na.rm=TRUE )   ),
                  c( mean(as.vector(DRMSE[,1:Nsim]) , na.rm=TRUE ), mean(as.vector(DRMSE[,Nsim+1:Nsim]) , na.rm=TRUE ) ,mean(as.vector(DRMSE[,2*Nsim+1:Nsim]) , na.rm=TRUE )   ),
                  c( mean(as.vector(DRcov[,1:Nsim]) , na.rm=TRUE ), mean(as.vector(DRcov[,Nsim+1:Nsim]) , na.rm=TRUE ) ,mean(as.vector(DRcov[,2*Nsim+1:Nsim]), na.rm=TRUE  )   )
  )  #MC: MSE and Coverage

  ##visualization
  ggdata<-data.frame(Stratum=1:10 ,
                     ests=c(as.vector(REst),as.vector(DREst)) ,
                     Stratification =c(  rep('Residual method', 3*10*Nsim  ),rep('Doubly-ranked method', 3*10*Nsim )),
                     XYmodel= rep(  paste0(theZXmodel,1:3) , each=10*Nsim   )
  )
  Mean<-function(x){mean(x,na.rm=TRUE)  }
  RTargetmean<-c( apply(RTarget[,1:Nsim] , 1, Mean), apply(RTarget[,Nsim+1:Nsim] , 1, Mean) ,apply(RTarget[,2*Nsim+1:Nsim] , 1, Mean)   )
  DRTargetmean<-c( apply(DRTarget[,1:Nsim] , 1, Mean), apply(DRTarget[,Nsim+1:Nsim] , 1, Mean) ,apply(DRTarget[,2*Nsim+1:Nsim] , 1, Mean)   )
  ggtrue<- data.frame(   Stratum=1:10  ,
                         tragets=c(  RTargetmean,DRTargetmean ) ,
                         Stratification =c(  rep('Residual method', 3*10 ),rep('Doubly-ranked method', 3*10 )),
                         XYmodel= rep( paste0(theZXmodel,1:3) , each=10  )
  )

  ###storage
  write.csv(MCtable  , paste0(path,theIVtype,'\\',"MCtable_",theZXmodel  , '.csv'), row.names=F)
  write.csv(ggdata  , paste0(path,theIVtype,'\\',"ggdata_",theZXmodel  , '.csv'), row.names=F)
  write.csv(ggtrue  , paste0(path,theIVtype,'\\',"ggtrue_",theZXmodel  , '.csv'), row.names=F)
}




##Table and visualization--------------------------------------------------
##-------------------------------------------------------------------------

for( theIVtype in c( 'cont', 'bi' , 'high-dim')){
  MCTABLE<-c()
  if(theIVtype!='high-dim'){
    theZXmodelcandidates<-c('A','B','C','D') }else{
      theZXmodelcandidates<-c('E','F','G','H')
    }

  for( theZXmodel in theZXmodelcandidates  ){
    MCtable<-read.csv(  paste0(path, theIVtype,"\\MCtable_",theZXmodel  , '.csv'),header=T, na.strings ="?")
    MCTABLE<-rbind( MCTABLE , MCtable )
    ggdata<-read.csv(  paste0(path, theIVtype,"\\ggdata_",theZXmodel  , '.csv'),header=T, na.strings ="?")
    ggtrue<-read.csv(  paste0(path, theIVtype,"\\ggtrue_",theZXmodel  , '.csv'),header=T, na.strings ="?")
    ggdata[,2]<-as.numeric(ggdata[,2])
    ggtrue[,2]<-as.numeric(ggtrue[,2])

    if(theZXmodel=='A'){
      scales_y <- list(
        `A1` = scale_y_continuous(limits = c(-1, 1)),
        `A2` = scale_y_continuous(limits = c(-1, 1)),
        `A3` = scale_y_continuous(limits = c(-1, 0.5))
      )
    }
    if(theZXmodel=='B'){
      scales_y <- list(
        `B1` = scale_y_continuous(limits = c(-0.5, 0.5)),
        `B2` = scale_y_continuous(limits = c(-1, 0.6)),
        `B3` = scale_y_continuous(limits = c(-0.6, 0.5))
      )
    }
    ####similar for 'C', 'D', 'E', ...

    p <- ggplot(ggdata, aes(Stratum, ests))+
      geom_boxplot(aes(group = cut_width(Stratum, 0.25)))+
      labs(x='Strata Rank',y='LACE estimates')+
      geom_point( data=ggtrue, mapping=aes( x=Stratum, y=tragets    ),col='red'  )+
      geom_hline(yintercept = 0,linetype='dashed')+
      facet_grid_sc(rows=vars(XYmodel),cols=vars(Stratification) , scales = list(y = scales_y))+
      scale_x_continuous(breaks = c(1, 2, 3,4,5 ,6, 7,8,9,10      ) )+
      theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    p
    ggsave(paste0(theZXmodel,'.eps' ),
           plot = p ,
           path=paste0(path, theIVtype,'\\'),
           height = 6, width = 8, units = "in",limitsize=TRUE)
  }
  print(theIVtype)
  print(  xtable(MCTABLE,digits=3) )
}


##Different U-Y relationship-----------------------------------------------------------------
##-------------------------------------------------------------------------------------------
Nsim<-1000
DREst<-c()
DRTarget<-c()
DRMSE<-c()
DRcov<-c()  #coverage rate
for(theUYmodel in c('U1','U2','U3')){
  print(theUYmodel)
  for(j in c('1','2','3')){
    print(j)
    for(i in 1:Nsim){  #Nsim  simulation times
      dat<-getDatU( IVtype='cont',XYmodel=j,UYmodel=theUYmodel ) ;rdat<-Stratify(dat)
      RES<-getSummaryInf( rdat,XYmodel=j,onlyDR=TRUE)
      DRres<-RES$DRres
      DREst<-cbind(DREst, DRres$est)
      DRTarget<-cbind( DRTarget, DRres$target )
      DRMSE<-cbind(DRMSE, (DRres$est-DRres$target)^2)#实际上是SE: squared errors
      DRcov<-cbind(DRcov,  abs(  DRres$est-DRres$target )<1.96*DRres$se)
    }
  }
}
dim(DREst)

Mean<-function(x){mean(x,na.rm=TRUE)  }
MSEvalues<-c();Cvalues<-c()
for(i in 1:9){
  MSEvalues<-c(MSEvalues,   Mean( as.vector(DRMSE[,(i-1)*Nsim+1:Nsim])  )  )
  Cvalues<-c( Cvalues , Mean( as.vector(DRcov[,(i-1)*Nsim+1:Nsim])  )  )
}


##table
MCtable<-cbind(matrix( MSEvalues,3,3  ) , matrix(Cvalues,3,3) )   #MC: MSE and Coverage
MCtable<-MCtable[,c(1,4,2,5,3,6)]
xtable(MCtable,digits=3)

##visualization
ggdata<-data.frame(Stratum=1:10 ,
                   ests=as.vector(DREst) ,
                   UYmodel =rep(  c('U1','U2','U3') , each=3*10*Nsim   ),
                   XYmodel= rep(  c('1','2','3') , each=10*Nsim   )
)
DRTargetmean<-c()
for(i in 1:9){
  DRTargetmean<-c(DRTargetmean, apply(DRTarget[,(i-1)*Nsim+1:Nsim ] , 1, Mean) )
}

ggtrue<- data.frame(   Stratum=1:10  ,
                       tragets=DRTargetmean ,
                       UYmodel =rep(  c('U1','U2','U3') , each=3*10   ),
                       XYmodel= rep(  c('1','2','3') , each=10   )
)
scales_y <- list(
  `1` = scale_y_continuous(limits = c(-1, 1)),
  `2` = scale_y_continuous(limits = c(-1, 1)),
  `3` = scale_y_continuous(limits = c(-1, 1)))

p <- ggplot(ggdata, aes(Stratum, ests))+
  geom_boxplot(aes(group = cut_width(Stratum, 0.25)))+
  labs(x='Stratum',y='LACE estimates')+
  geom_point( data=ggtrue, mapping=aes( x=Stratum, y=tragets    ),col='red'  )+
  geom_hline(yintercept = 0,linetype='dashed')+
  facet_grid_sc(rows=vars(XYmodel),cols=vars(UYmodel) , scales = list(y = scales_y))+
  scale_x_continuous(breaks = c(1, 2, 3,4,5 ,6, 7,8,9,10      ) )+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p





##linearity testing--------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------
devtools::install_github("JonSulc/PolyMR")
library(PolyMR)

Nsim<-1000

#one specific example only - other examples can be considered as well
IVtype_used<-'bi'
ZXmodel_used<-'A'
XYmodel_used<-'2'

set.seed(1123)
pvalues<-c()
for(i in 1:Nsim){
    dat<-getDat( IVtype=IVtype_used, ZXmodel=ZXmodel_used,XYmodel=XYmodel_used )
    rdat<-Stratify(dat); RES<-getSummaryInf( rdat,XYmodel='1')
    pvalues<-c(pvalues, c(RES$DRHeterQ[3],RES$RHeterQ[3] ) )
    polymr_results <- polymr( exposure=dat$X , outcome=dat$Y, genotypes=cbind(dat$Z) )
    pvalues<-c(pvalues,polymr_results$polymr$pval_linear_model)
}
Pvalues<-matrix(pvalues, 3,Nsim  )
print(   apply(Pvalues<0.05 , 1,mean)   )







