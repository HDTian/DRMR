###Stratification

#注意： Stratify()这个函数是没有randomness的！所以对于DR method (Residual method易知稳定性很高)：
#如果想用same-value IV (eg binary IV)随机性；dat在之前就要打乱排序一次
#如果是在continuous IV + same-valued exposure下体现随机性，那么最开始打乱是没用的，因为一定会按照IV的顺序排X
#此时一种解决方案是给给same-valued exposure增加一个特别细细细微的扰动值
#这种操作也适用于same-valued IV情况
#如果存在Missing data，dat之前就需要自主抹掉missing data
#如果想保证stratum的数量绝对相等，dat在之前就需要控制total sample size

#Stratify中 Residual和DR stratification均会被使用

Stratify<-function(dat,  #must be a data frame   contains $Z $X $Y $M (一旦设置了，就会用这个); no missing data
                   onExposure=TRUE,  #TRUE for exposure; otherwise/FALSE on $M
                   Ns=10, #the number of (final) stratum
                   SoP=NA,
                   seed=NA){

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

  #permutation or not?
  if(!is.na(seed)){
    set.seed(seed)
    per_index<-sample(  1:nrow(dat),nrow(dat)  )
    dat<-dat[per_index,]  #permutation here is to induce randomness for same-valued IV
    Mvalues<-sort(as.numeric(levels( factor( dat$M )  )) )  #possible M values
    Mdif<- min( Mvalues[-1] - head(Mvalues,-1) ) #minimal difference between non-same M values
    dat$Merror<-runif(  nrow(dat), -Mdif/2, Mdif/2 ) #the M error will not affect the rank at all
    dat$Msave<-dat$M
    dat$M<-dat$M+dat$Merror #之后会减去
  }

  N<-nrow(dat)


  #ordered by Z
  dat_order<-dat[ order(dat$Z  ),  ]   #这里不存在random，对于相同的值order()就是按数字出现顺序排的(要想随机input dat可以打乱一下)
  dat_order$ID<-1:N #IV排序好的状态：方便回标(in case, 但往往不需要)



  ###residual stratification--------------------------------------------------------------
  ###-------------------------------------------------------------------------------------
  dat_order$residual<-dat_order$M-lm(  dat_order$M~dat_order$Z   )$fitted

  #之前的floor()有一点不好：不能绝对控制Ns的数量，所以还是用排序的思路
  #dat_order$Rstratum<-  floor( (rank( dat_order$residual,ties.method ='random' )/((N/2)+0.000000001) ) )+1 #N为奇数也不要紧啊

  dat_order<-dat_order[  order(dat_order$residual  )  ,  ]
  dat_order$Rstratum<- sort(    rep(1:Ns, length.out=N)  )#一定可以控制Ns


  ###Doubly-ranked stratification----------------------------------------------------------
  ###--------------------------------------------------------------------------------------
  dat_order<-arrange(  dat_order ,ID )#回到之前的IV排序状态

  #之前的No<-round(N/SoP)有一点不好：不能绝对控制SoP的大小，所以还是用排序的思路
  #No<-round(N/SoP)  #SoP是我们可以控制的！ No代表着pre-stratum的数量
  if(is.na(SoP)){SoP<-Ns}
  dat_order$pre_stratum<-    rep(1:(floor(N/SoP)+1), each=SoP,length.out=N)
  #(floor(N/SoP)+1)*SoP >= N 保证能超过就行

  #rank twice (ie doubly-ranked)
  temp<-arrange(  dat_order, M )  #按照M升序排一下  #arrange()应该也没有随机性;尤其M是在same-valued的情况下
  dat_order<-arrange(  temp ,pre_stratum ) #即，保证pre_strata按顺序排列，并且每个pre_strata中的目标量都是升序

  dat_order$DRstratum<-as.vector( unlist(sapply( as.numeric(table( dat_order$pre_stratum )) , function(x) sort(rep(1:Ns,length.out=x)) )   ) )

  if(!is.na(seed)){
    dat_order$M<-dat_order$Msave  #back to the original M
    dat_order<-dat_order[ , !(names(dat_order) %in% c('Merror','Msave'))] #remove the Merror and Msave column
  }

  return(dat_order)
}


#结果是rdat, 一种基于dat的增广data



###3种stratification： R DR N

#假设Ns都目标为10

#N就相当于DR中把SoP设置=total sample size; 即：Stratify(dat,SoP=nrow(dat))

#R和N类似，只不过目标量变成了residual ，可以在inputted data中 dat$M<-dat$M-lm(  dat$M~dat$Z   )$fitted
#然后针对M作N即可(即，相当于DR中把SoP设置=total sample size) 即：Stratify(dat,onExposure=FALSE,SoP=nrow(dat))

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

