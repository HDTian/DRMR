###Stratification

#ע�⣺ Stratify()���������û��randomness�ģ����Զ���DR mthod (Residual method��֪�ȶ��Ժܸ�)��
#�������same-value IV (eg binary IV)����ԣ�dat��֮ǰ��Ҫ��������һ��
#�������continuous IV + same-valued exposure����������ԣ���ô�ʼ������û�õģ���Ϊһ���ᰴ��IV��˳����X
#��ʱһ�ֽ�������Ǹ���same-valued exposure����һ���ر�ϸϸϸ΢���Ŷ�ֵ
#���ֲ���Ҳ������same-valued IV���
#�������Missing data��dat֮ǰ����Ҫ����Ĩ��missing data
#����뱣֤stratum������������ȣ�dat��֮ǰ����Ҫ����total sample size

#Stratify�� Residual��DR stratification���ᱻʹ��

Stratify<-function(dat,  #must be a data frame   contains $Z $X $Y $M (һ�������ˣ��ͻ������); no missing data
                   onExposure=TRUE,  #TRUE for exposure; otherwise/FALSE on $M
                   Ns=10, #the number of (final) stratum
                   SoP=NA){
  #dat check
  if(!is.data.frame(dat)){ stop('Your inputted data is not a dataframe')  }
  if( is.null(dat$Z) ){stop('No Z detected!: The column of the instrument should be named as Z') }
  if( is.null(dat$X) ){stop('No X detected!: The column of the exposure should be named as X') }

  #other logic check
  if( (!onExposure)&(is.null(dat$M)) ){
    stop('You wish to stratify on a covariate, but not defined in your data; try to defined them as $M or stratify on exposure only'  )
  }
  if( (onExposure)&(!is.null(dat$M))   ){
    print('Because you wish to stratify on the exposure and you have defined M in your data, your M information will be replaced by the exposure')
  }
  if(onExposure ){
    dat$M<-dat$X
  }

  #missing data check
  if( sum( is.na(dat) )!= 0  ){ stop( 'Missing data exists, please clean the data first'  ) }

  N<-nrow(dat)


  #ordered by Z
  dat_order<-dat[ order(dat$Z  ),  ]   #���ﲻ����random��������ͬ��ֵorder()���ǰ����ֳ���˳���ŵ�(Ҫ�����input dat���Դ���һ��)
  dat_order$ID<-1:N #IV����õ�״̬������ر�(in case, ����������Ҫ)



  ###residual stratification--------------------------------------------------------------
  ###-------------------------------------------------------------------------------------
  dat_order$residual<-dat_order$M-lm(  dat_order$M~dat_order$Z   )$fitted

  #֮ǰ��floor()��һ�㲻�ã����ܾ��Կ���Ns�����������Ի����������˼·
  #dat_order$Rstratum<-  floor( (rank( dat_order$residual,ties.method ='random' )/((N/2)+0.000000001) ) )+1 #NΪ����Ҳ��Ҫ����

  dat_order<-dat_order[  order(dat_order$residual  )  ,  ]
  dat_order$Rstratum<- sort(    rep(1:Ns, length.out=N)  )#һ�����Կ���Ns


  ###Doubly-ranked stratification----------------------------------------------------------
  ###--------------------------------------------------------------------------------------
  dat_order<-arrange(  dat_order ,ID )#�ص�֮ǰ��IV����״̬

  #֮ǰ��No<-round(N/SoP)��һ�㲻�ã����ܾ��Կ���SoP�Ĵ�С�����Ի����������˼·
  #No<-round(N/SoP)  #SoP�����ǿ��Կ��Ƶģ� No������pre-stratum������
  if(is.na(SoP)){SoP<-Ns}
  dat_order$pre_stratum<-    rep(1:(floor(N/SoP)+1), each=SoP,length.out=N)
  #(floor(N/SoP)+1)*SoP >= N ��֤�ܳ�������

  #rank twice (ie doubly-ranked)
  temp<-arrange(  dat_order, M )  #����M������һ��  #arrange()Ӧ��Ҳû�������;����M����same-valued�������
  dat_order<-arrange(  temp ,pre_stratum ) #������֤pre_strata��˳�����У�����ÿ��pre_strata�е�Ŀ������������

  dat_order$DRstratum<-as.vector( unlist(sapply( as.numeric(table( dat_order$pre_stratum )) , function(x) sort(rep(1:Ns,length.out=x)) )   ) )

  return(dat_order)
}


#�����rdat, һ�ֻ���dat������data



###3��stratification�� R DR N

#����Ns��Ŀ��Ϊ10

#N���൱��DR�а�SoP����=total sample size; ����Stratify(dat,SoP=nrow(dat))

#R��N���ƣ�ֻ����Ŀ���������residual ��������inputted data�� dat$M<-dat$M-lm(  dat$M~dat$Z   )$fitted
#Ȼ�����M��N����(�����൱��DR�а�SoP����=total sample size) ����Stratify(dat,onExposure=FALSE,SoP=nrow(dat))

##example
# dat<-getDat( IVtype='cont', ZXmodel='C',XYmodel='2' ) #get a toy data
#
# #DR
# rdat<-Stratify(dat) ; getSummaryInf( rdat,target=FALSE)
# #Naive
# rdat<-Stratify(dat,SoP=nrow(dat)) ; getSummaryInf( rdat,target=FALSE,onlyDR=TRUE)
# #Residual
# dat$M<-dat$X-lm(  dat$X~dat$Z   )$fitted
# rdat<-Stratify(dat,SoP=nrow(dat),onExposure=FALSE) ; getSummaryInf( rdat,target=FALSE,onlyDR=TRUE)
