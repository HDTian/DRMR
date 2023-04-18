
###get Gelman-Rubin statistics


#GR����������coarsened�Ƿ�����(��Ҫ��ָ�ȶ���Ҫ��Ȼcoarsenend trend��ǿ������violaten local exchangeability)


#�����������stratum�����ɣ���Ȼ��ÿ��single--pre-stratum--individual stratumҲmake sense


getGRstats<-function(rdat, #rdat<-Stratify(dat)
                     Nc=2, #number of chains calculating GR stats
                     roundnum=3 #GR results round precision
){
  #input��rdat #rdat<-Stratify(dat)
  Np<-length(  levels(factor(rdat$pre_stratum))  )  #number of pre-stratum

  #������һ��pre-stratum size�͵����ڶ�����һ������Ĩȥ
  if( table( rdat$pre_stratum )[Np]!=table( rdat$pre_stratum )[Np-1]  ){
    rdat<-rdat[ rdat$pre_stratum!=Np  , ]
    Np<-Np-1  #�Ƴ�end pre-stratum��� # of pre-stratum
  }


  ###Build A matrix: the ready-to-go matrix
  #ע�⣺rdat���Ѿ�����˫�ź��dat��(����)matrix�ˣ���pre-stratum˳�����ӣ�ÿ��pre-stratum�ڲ���exposureҲ��˳������
  A<-matrix( rdat$X,   nrow(rdat)/Np   , Np            )  #nrow(rdat)/Np  ����SoP:  Size of Pre-stratum


  #lowfun upfun: get the low/up index for any inputted Vector
  lowfun<-function( vec  ){
    rep(      c(1,head(cumsum(as.numeric(table(vec)))+1,-1) )     ,   as.numeric(table(vec))     ) #checked
  }
  upfun<-function( vec  ){
    rep(      cumsum(as.numeric(table(vec)))    ,   as.numeric(table(vec))     )  #checked
  }

  Lows<-apply(A, 2    ,lowfun ) #the low/up index
  Ups<-apply(A, 2    ,upfun )

  #Nc: number of chains for low/up index in GR statistics analysis
  #ncol(A)=length(Lows)=Np #Np�Ѿ���updated֮�����
  Lows_<-Lows[,sort(sample( 1:Np,floor(Np/Nc)*Nc))  ]#����Ƴ���һЩ��ʹ�ÿ�������split to Nc equally-length chains
  Ups_<-Ups[,sort(sample( 1:Np,floor(Np/Nc)*Nc))  ]

  #find the position (from 1:SoP) that we need Low or Up
  rdat_<-rdat[rdat$pre_stratum==1,]
  low_p<-c(1,   head( cumsum(as.numeric(table( rdat_$DRstratum )))+1 ,-1 ) ) #low_p where p means position
  up_p<- cumsum(as.numeric(table( rdat_$DRstratum )))

  GR_low<-c() #GR_low GR_up������== Ns;ÿ��Stratum������һ����Ӧ��low��up GRֵ (in the low_p up_p position)
  GR_up<-c()
  Ns<-length(  levels(factor(rdat$DRstratum))  )
  for(i in 1:Ns){
    ###for Lows----------------------------------

    p<-low_p[i]   #p: position (from 1:SoP)

    chains<-matrix( Lows_[p,],  ncol=Nc ) #vector��matrixʱ��������
    cm<-apply(chains, 2, mean )  #chain mean
    tm<-mean( cm )   #total mean
    B<-sum((cm-tm)^2)*(  nrow(chains)/(Nc-1) ) #Between-chain variance
    cv<-apply(chains, 2, var )  #chain variance
    W<-mean(cv)  #W: Within chain variance
    V<- (nrow(chains)-1)/nrow(chains)*W + 1/nrow(chains)*B
    if(W==0){ GR<-1 }else{GR<-sqrt(   V/W ) }

    GR_low<-c(GR_low, round(GR,roundnum)  )

    ###for Ups------------------------------------

    p<-up_p[i]   #p: position (from 1:SoP)

    chains<-matrix( Ups_[p,],  ncol=Nc ) #vector��matrixʱ��������
    cm<-apply(chains, 2, mean )  #chain mean
    tm<-mean( cm )   #total mean
    B<-sum((cm-tm)^2)*(  nrow(chains)/(Nc-1) ) #Between-chain variance
    cv<-apply(chains, 2, var )  #chain variance
    W<-mean(cv)  #W: Within chain variance
    V<- (nrow(chains)-1)/nrow(chains)*W + 1/nrow(chains)*B
    if(W==0){ GR<-1 }else{GR<-sqrt(   V/W ) }

    GR_up<-c(GR_up, round(GR,roundnum)  )
  }

  GRRES<-cbind( paste0('Stratum', 1:Ns  ), GR_low, GR_up  ,pmax(GR_low,GR_up )  )
  colnames(GRRES)<-c('Stratum', 'GR_low', 'GR_up', 'maxGR' )
  return(GRRES)
}


##example
# getGRstats(rdat,Nc=2,roundnum=3)
# getGRstats(rdat,Nc=5,roundnum=3)
# getGRstats(rdat,Nc=10,roundnum=3)
# getGRstats(rdat,Nc=100,roundnum=3)



