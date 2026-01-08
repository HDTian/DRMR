###

#creat toy data with ZXmodel (only applicable for IV type='cont' or 'bi')
#but with different UYmodel

getDatU<-function(N=10000,
                 IVtype='cont',   #only for 'cont' 'bi' and only ZXmodel='A'
                 UYmodel='U1',  #'U1' 'U2' 'U3'
                 XYmodel='1', #'1' '2' '3'
                 printRR=FALSE  #print the RR results (ie instrument strength) of this data?
){

  if(!IVtype%in%c( 'cont', 'bi'  ) ){stop('please use the correct IV type')}
  if(IVtype=='cont'){
    Z<-rnorm( N,0,0.5   )
  }
  if(IVtype=='bi'){
    Z<-rbinom( N,1,0.3   )
  }




  U1<-rnorm(N,0,1)  #epsilon_X
  U2<-rnorm(N,0,1)   #U epsilon 同作
  U3<-rnorm(N,0,1)   #epsilon_Y



  ##Z-X model - only ZXmodel='A' i.e. h_X<-function(Z){ 0.5*Z  }

  X<-0.5*Z+U2+U1


  #if(!exists('X_')){X_<-X} #不会使用到coarsend exposure X-Y model




  ##X-Y model


  if(!XYmodel%in%c( '1','2','3','4' ) ){stop('please use the correct XYmodel index')}
  if(XYmodel=='1'){
    nonlinear<-function(t){0*t} #constant effect; i.e. linear effect
    nonlinear_dif<-function(t){ sapply(   t , function(tt){ 0}) }
  }
  if(XYmodel=='2'){
    nonlinear<-function(t){0.1*t^2}
    nonlinear_dif<-function(t){0.2*t}
  }
  if(XYmodel=='3'){
    nonlinear<-function(t){-0.1*(t)^2*(t>0) }
    nonlinear_dif<-function(t){-0.2*t*(t>0)}
  }
  # if(XYmodel=='4'){
  #   nonlinear<-function(t){0.3*t^2}
  #   nonlinear_dif<-function(t){0.6*t}
  # }

  if(!UYmodel%in%c( 'U1','U2','U3' ) ){stop('please use the correct UYmodel index')}
  if(UYmodel=='U1'){ Y<-nonlinear(X)+U2^2+U3  }
  if(UYmodel=='U2'){ Y<-nonlinear(X)+abs(U2)+U3  }
  #if(UYmodel=='U3'){ Y<-nonlinear(X)+U2-U1+U3  }
  if(UYmodel=='U3'){ Y<-nonlinear(X)+U2+U1^2+2*U1*U2+U3  }




  fitX<-lm( X~Z   )
  RR<-as.numeric(summary(fitX)$r.squared[1])
  if(printRR){ print(  paste0('Coefficient of determination: ',round(as.numeric(RR),3)) )   }


  dat<-data.frame( Z=Z,X=X,Y=Y,U=U2 )
  return(  dat )
}

# dat<-getDatU(IVtype='cont', XYmodel='1',UYmodel='U1',printRR=TRUE )



