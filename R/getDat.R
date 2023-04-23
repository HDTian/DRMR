###aaa

###

###

getDat<-function(N=10000,
                 IVtype='cont',   #'cont' 'bi' 'high-dim'
                 ZXmodel='A',  #A B C D (for IVtype='cont' or 'bi') E F G H (for IVtype='high-dim')
                 XYmodel='1',
                 printRR=FALSE  #print the RR results (ie instrument strength) of this data?
){

  if(!IVtype%in%c( 'cont', 'bi' ,'high-dim' ) ){stop('please use the correct IV type')}
  if(IVtype=='cont'){
    Z<-rnorm( N,0,0.5   )
  }
  if(IVtype=='bi'){
    Z<-rbinom( N,1,0.3   )
  }
  if(IVtype=='high-dim'){
    G_1<-rbinom( N,1,0.3   ) ; G_2<-rbinom( N,1,0.3   ) ;G_3<-rbinom( N,1,0.3   )
  }



  U1<-rnorm(N,0,1)
  U2<-rnorm(N,0,1)   #U epsilon Í¬×÷
  U3<-rnorm(N,0,1)



  ##Z-X model

  if(!ZXmodel%in%c( 'A', 'B', 'BC', 'C','D','E','F','G','H' ) ){stop('please use the correct ZXmodel index')}

  if( (IVtype %in% c('cont', 'bi')  )&(ZXmodel%in%c('E','F','G','H')) ){
    stop('please use the appropriate combination of the IVtype and the ZXmodel')
  }
  if( (IVtype == 'high-dim'  )&(ZXmodel%in%c('A','B','C','D')) ){
    stop('please use the appropriate combination of the IVtype and the ZXmodel')
  }



  if(ZXmodel=='A'){
    h_X<-function(Z){ 0.5*Z  }
    X<-h_X(Z)+U2+U1            }
  if(ZXmodel=='B'){
    h_X<-function(Z){   0.5*Z+  (2*Z^3+2)*(Z<(-1))     }
    X<-h_X(Z)+U2+U1            }
  if(ZXmodel=='BC'){
    X<-exp(0.5*Z+U2+U1)            }      #strictly positive exposure
  if(ZXmodel=='C'){
    h_X<-function(Z,U){ -10+ (1.5+0.4*U)*(Z+5) }
    X<-h_X(Z,U2)+U2 +U1        }
  if(ZXmodel=='D'){
    h_X<-function(Z){ 0.5*Z  }
    X<-h_X(Z)+U2+U1
    X_<-X  # true continuous X
    X<-round(X)      }

  if(ZXmodel=='E'){
    X<-0.2*G_1+0.4*G_2+0.6*G_3+U2+U1            }
  if(ZXmodel=='F'){
    X<-0.1*G_1+2.0*G_2*G_3+U2+U1            }
  if(ZXmodel=='G'){
    X<--10+(5+0.2*G_1+0.4*G_2+0.6*G_3)*(1.5+0.5*U2)+U2+U1        }
  if(ZXmodel=='H'){
    X<-0.2*G_1+0.4*G_2+0.6*G_3+U2+U1
    X_<-X  # true continuous X
    X<-round(X)      }

  if(ZXmodel%in%c('E','F','G','H')){ #for high-dim IV, build a single weighted gene score
    GXfit<-lm(X~ 1+G_1+G_2+G_3  )
    Z<-cbind(G_1,G_2,G_3)  %*%as.numeric(coef(GXfit))[-1]  #Z: wighted gene score
  }

  if(!exists('X_')){X_<-X}




  ##X-Y model


  if(!XYmodel%in%c( '1','2','3' ) ){stop('please use the correct XYmodel index')}
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



  Y<-nonlinear(X)+U2+U3
  R<-c( min(X) , max(X )  )
  vals<-nonlinear_dif( seq( R[1],R[2], len=100  )  )
  RC<-c(min(vals)-0.5,max(vals)+0.5  )

  fitX<-lm( X~Z   )
  RR<-as.numeric(summary(fitX)$r.squared[1])
  if(printRR){ print(  paste0('Coefficient of determination: ',round(as.numeric(RR),3)) )   }


  dat<-data.frame( Z=Z,X=X,Y=Y,U=U2, X_=X_ )
  return(  dat )
}

# rdat<-getDat(IVtype='high-dim', ZXmodel='G',printRR=TRUE )



