# DRMR
**Doubly-Ranked Stratification in Mendelian Randomization**

Doubly-ranked (DR) stratification is a nonparametric, simple yet powerful method for instrument variable (IV) analysis and Mendelian randomization (MR) studies to create sub-groups, known as "strata", in which the IV assumption is satisfied. DR stratification can be applied in a wide range of IV or MR studies, including assessing homogeneity assumption, conducting nonlinear causal studies, and examining heterogeneous effects. The DRMR package can assist in generating stratifications, providing relevant results, and supporting further analysis based on stratification methods. 

This manuscript will demonstrate the fundamental concepts of DR stratification and provide a simple, step-by-step guide for applying this method. Additionally, a series of subsequent papers will be published to provide further details and insights.[^1].
[^1]: Different papers would use different writing style for illustration and interpretation. Contact me if you are confused or interested in any aspect of the DR stratification method.


Make sure to download the DRMR package

```R
devtools::install_github("HDTian/DRMR")
library(DRMR)
```



## Stratification
When examining potential non-linear causal effects, one approach is to divide the population into multiple subgroups or strata, each with varying levels of exposure. IV analysis can then be applied to each stratum, and the stratum-specific IV estimates can reveal the underlying shape of the causal effect. There are three main stratification methods:
|             | Naive stratification |Residual stratification[^R] |Doubly-ranked stratification[^DR] |
| ----------- | ----------- | ----------- | ----------- |
| Discription | Directly stratify on the exposure level.  |First build the 'residual exposure' by regressing the exposure on the instrument and do the naive-style stratification on the residual values      |First rank the instrument values to form pre-stratum, then rank the exposure values to form stratum (see schematic diagram below)     |
| Potential Issues    | The exposure is often the common effect of the instrument and the confounders. Therefore, stratifying or conditioning on the common effect may introduce collider bias[^collider], and violate the exchangeability of the strata.      | Residual stratification requires strong parametric assumptions; for exmaple, the linearity and homogeneity assumption of the instrument-exposure model |       |

[^collider]: https://academic.oup.com/ije/article/39/2/417/680407
[^R]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4222800/
[^DR]: https://www.biorxiv.org/content/10.1101/2022.06.28.497930v1.abstract
<img alt="DR" src="https://user-images.githubusercontent.com/127906571/229943653-cefbd6ad-fcd5-45a4-95b4-6530b4ea9836.png">

Use the function `getDat()` to create a toy example, where you can experiment with different types of instruments (IVtype), instrument-exposure models (ZXmodel), and exposure-outcome models (XYmodel).
```R
dat<-getDat( IVtype='cont', ZXmodel='C',XYmodel='2' )
```

  

Draw the doubly-ranked stratification (use `?Stratify` to check the parameters for stratification)
```R
rdat<-Stratify(dat) 
RES<-getSummaryInf( rdat,target=FALSE)
head(RES$DRres)
```
```R
#          size        min        1q      mean        3q        max       bx       bxse         by       byse        est         se target
#Stratum 1 1000 -15.555421 -8.660667 -7.427431 -6.035485 -1.3493712 1.066652 0.12171101 -1.4378919 0.17301791 -1.3480425 0.16220656     NA
#Stratum 2 1000 -10.439764 -6.721907 -5.691190 -4.521686 -0.7632903 1.158202 0.09373183 -1.3370736 0.10698333 -1.1544386 0.09237015     NA
#Stratum 3 1000  -9.775878 -5.546970 -4.584672 -3.578259 -0.0726714 1.241074 0.08080925 -1.0823088 0.08425849 -0.8720746 0.06789162     NA
#Stratum 4 1000  -8.834297 -4.708556 -3.717982 -2.779919  0.3667541 1.285770 0.07975024 -0.8680834 0.07767019 -0.6751469 0.06040754     NA
#Stratum 5 1000  -7.772422 -3.887295 -2.920028 -1.966694  3.0986642 1.328140 0.07838736 -0.7696176 0.06869428 -0.5794700 0.05172215     NA
#Stratum 6 1000  -7.741110 -3.107658 -2.138776 -1.195922  3.3621614 1.441501 0.07931447 -0.6188700 0.06872492 -0.4293235 0.04767596     NA
```
The function `Stratify()` will return the individual stratification index. The function `getSummaryInf()` will provide summary information based on the stratified result.

Note that the naive and residual stratification can be regarded as special scenarios of the DR stratification. Hence, you may prefer to use only the DR method..

For naive stratification, it is equivalent to applying DR with one pre-stratum:
```R
rdat<-Stratify(dat,SoP=nrow(dat)) 
getSummaryInf( rdat,target=FALSE,onlyDR=TRUE)
```
For residual stratification, it is equivalent to applying the naive stratification wrt the residual variables
```R
dat$M<-resid(lm(  dat$X~dat$Z   )) #obtain the residual variables first
rdat<-Stratify(dat,SoP=nrow(dat),onExposure=FALSE) 
getSummaryInf( rdat,target=FALSE,onlyDR=TRUE)
```


## Stratify on other covariates
If you need to stratify on a covariate (which is common in many heterogeneous effect studies[^HTE]), first define the covariate to be stratified by as `M` in `dat`, and then simply run the following code.
```R
dat$M<-covarite_vector  #added the covariate information
rdat<-Stratify(dat,onExposure=FALSE) 
getSummaryInf( rdat,target=FALSE)
```

[^HTE]: https://link.springer.com/article/10.1007/s10654-022-00879-0

## Stratify on coarsened exposures
The coarsened exposure refers to the exposure measured with discrete values that are considered as an approximation for its latent continuous values. Such values could be rounded, binned into categories, or truncated, and all of them should satisfy the rank preserving condition[^444]. The properties of such coarsened exposure enable the use of doubly-ranked stratification, but extra assessment for the degree of coarseness should be done before stratification. One way to do this is by calculating the Gelman-Rubin uniformity values of the rank index in each stratum. You can run the following code to perform this assessment.
```R
getGRstats(rdat)

```
```R
#      Stratum     GR_low  GR_up   maxGR  
#      "Stratum1"  "1"     "1.001" "1.001"
#      "Stratum2"  "1.001" "1"     "1.001"
#      "Stratum3"  "1"     "1"     "1" 
#       ...

```
small values (<1.02) indicate a small degree of coarsenness.


[^444]: See the Supplementary Text S2 of the origina doubly-ranked stratification paper


## Different types of instrument
The doubly-ranked stratification supports a wide variaty of instrument types, including dichotomous IV, discrete IV, high-dimensional IV, continuous IV, mixture IVs. For single IV, just denote its value by `Z` in `dat`. For high-dimensional IV, integrate them as a weighted score (e.g. the weigthed gene score in MR). For exmaple
```R
dat$Z<-lm(  dat$X~ dat$G1+ dat$G2 + ...   )$fitted  #G1 G2 ... are genetic variants
```

## Control the randomness and reproducibility
One feature of the stratification in the DRMR package is that the stratification function will not introduce randomness (in the sense that the same input data will always generate consistent stratification results)[^randomness]. This is beneficial in terms of transparency and reproducibility. However, since ties are broken at random for constant instrument or exposure values, the randomness of such breaks needs to be created before running `Stratify()`. For example, you can add negligible errors to the variable column containing constant values to induce randomness.
```R
set.seed(1) #track 
dat$Z<-dat$Z+rnorm(nrow(dat),0,1e-10 )  #when the instrument is discrete-valued
dat$X<-dat$X+rnorm(nrow(dat),0,1e-10 )  #when the exposure is discrete (eg coarsened)
Stratify(dat) #then run the DR stratification
 ```
[^randomness]: I considered whether to have the function perform break permutation automatically, but I decided against it as it might result in the user losing control over the variation of results caused by the randomness introduced internally (e.g., if they tend not to set a seed). Although these variations often do not affect any conclusions, I removed the randomness from inside the function to ensure transparency and reproducibility.
       


## Variable transformation
You may consider transformed exposures or covariates[^trans]. The transformed variable could be helpful in terms of interpretation. If you use a bijective transformation, the DR stratification will not be changed. You can simply transform you variable wrt `dat` by using `dat$X<-f(dat$X)` where `f()` is your transformation function`   
[^trans]: One example is the log-transformation: https://www.medrxiv.org/content/10.1101/2022.10.26.22280570v2 
## Smoothing
There are several methods for smoothing the stratification estimation (the stratum-specific estimates are called LACE), including the fractional polynomial method and the piecewise linear method[^smoothing].

To use the fractional polynomial method (e.g. with degree 2), run
```R
smooth_res<-Smooth(RES$DRres,Norder=3,baseline=0) #RES is the result returned by getSummaryInf()
```
If you wish to visualize the results, simply run `smooth_res$p ` or `smooth_res$hp`
![Rplot](https://user-images.githubusercontent.com/127906571/232224279-908bc7f7-8e38-4871-b4f5-0f19077d569a.png)

![Rplot02](https://user-images.githubusercontent.com/127906571/232224348-0605d2a0-9ff5-4c1a-beb6-c81f8c2c262b.png)




To use the piecewise linear method, run:
```R
cutting_values<-(RES$DRres$mean[-1] + head( RES$DRres$mean,-1) )/2
smooth_res<-Smooth(RES$DRres,Norder=1,baseline=0,Knots=cutting_values)
```
[^smoothing]: https://onlinelibrary.wiley.com/doi/full/10.1002/gepi.22041

## Real application
Assume you now have the real samples `dat` ready for stratificaton; simply run the code below
```R       
rdat<-Stratify(dat)
RES<-getSummaryInf(rdat,target=FALSE)
smooth_res<-Smooth(RES$DRres,Norder=3,baseline=0)
```      
then you can use the information in `RES` and `smooth_res` to build the results you desire[^further].

[^further]: Feel free to contact me if you have any further questions or need additional information that the current function cannot provide. I would be happy to discuss your study further and provide any assistance you may need.

