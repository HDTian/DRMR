# DRMR
**Doubly-Ranked Stratification in Mendelian Randomization**

*Latest updated: Apr/04/2025*

Doubly-ranked (DR) stratification is a nonparametric, simple yet powerful method for instrument variable (IV) analysis and Mendelian randomization (MR) studies to create sub-groups, known as "strata", in which the IV assumption is satisfied. DR stratification can be applied in a wide range of IV or MR studies, including assessing homogeneity assumption, conducting nonlinear causal studies, and heterogeneous effects studies. The DRMR package can assist in generating stratifications, providing relevant results, and supporting further analysis based on stratification results (also well connected with [RFQT](https://github.com/HDTian/RFQT) and [SSS](https://github.com/HDTian/SSS) packages). 

By using the DRMR package, you will activate stratification-based analysis, allowing you:  
✅ drawing nonlinear IV/MR analysis  
✅ estimate the exposure-outcome effect curve (supported by [SSS](https://github.com/HDTian/SSS))  
✅ detect and estimate the effect change-point (supported by [SSS](https://github.com/HDTian/SSS))  
✅ estimate individual and heterogeneous treatment effect (supported by [RFQT](https://github.com/HDTian/RFQT))  
✅ find the important effect drivers (supported by [RFQT](https://github.com/HDTian/RFQT))  
✅ evaluating homogeneity conditions



You can install the DRMR package in R:
```R
devtools::install_github("HDTian/DRMR")
```
```R
library(DRMR)
```
you can also try the [SUMnlmr](https://github.com/amymariemason/SUMnlmr) package, where the DR method was embedded, but SUMnlmr currently does not support effect heterogeneity and effect change-point analysis 

Related papers:  
[Relaxing parametric assumptions for non-linear Mendelian randomization using a doubly-ranked stratification method](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010823) (DR stratification)  
[Violation of the Constant Genetic Effect Assumption Can Result in Biased Estimates for Non-Linear Mendelian Randomization](https://karger.com/hhe/article/88/1/79/862215/Violation-of-the-Constant-Genetic-Effect) (why DR stratification outperform)  
[Stratification-based Instrumental Variable Analysis Framework for Nonlinear Effect Analysis](https://arxiv.org/html/2507.07349v1) (Nonlinear MR)  
[A data-adaptive method for investigating effect heterogeneity with high-dimensional covariates in Mendelian randomization](https://link.springer.com/article/10.1186/s12874-024-02153-1)
(high-dim effect heterogeneity)

⚠️ One-sample individual-level data is required. This may not be a drawback, as DRMR provides much more informative results that are not attainable through any summary-data method  
⚠️ All effect heterogeneity analyses can be robust to unreliable DR stratification. Even if DR stratification is unreliable, your effect heterogeneity results (those from RFQT) are typically not affected. The usage of DRMR and RFQT is not determined by whether or not DR stratification is reliable


The remaining part will demonstrate the fundamental concepts of DR stratification and provide a simple, step-by-step guide for applying this method for nonlinear Mendelian randomization. You may want to skip, and directly see  
👉 [Real application](#real-application) for the guidance with your real data  
👉 [FAQs](#faqs) for the frequently asked questions (remember you can always contact me for specific suggestions)

## Stratification
When examining potential non-linear causal effects, one approach is to divide the population into multiple subgroups or strata, each with varying levels of exposure. IV analysis can then be applied to each stratum, and the stratum-specific IV estimates can reveal the underlying shape of the causal effect. There are three main stratification methods:
|             | Naive stratification |Residual stratification |Doubly-ranked stratification |
| ----------- | ----------- | ----------- | ----------- |
| Discription | Directly stratify on the exposure level.  |First build the 'residual exposure' by regressing the exposure on the instrument and do the naive-style stratification on the residual values      |First rank the instrument values to form pre-stratum, then rank the exposure values to form stratum (see schematic diagram below)     |
| Potential Issues    | The exposure is often the common effect of the instrument and the confounders. Therefore, stratifying or conditioning on the common effect may introduce collider bias[^collider], and violate the exchangeability of the strata.      | Residual stratification requires strong parametric assumptions; for exmaple, the linearity and homogeneity assumption of the instrument-exposure model |       |

[^collider]: https://academic.oup.com/ije/article/39/2/417/680407

<img alt="DR" src="https://user-images.githubusercontent.com/127906571/229943653-cefbd6ad-fcd5-45a4-95b4-6530b4ea9836.png">

Use the function `getDat()` to create a toy example, where you can experiment with different types of instruments (`IVtype`), instrument-exposure models (`ZXmodel`), and exposure-outcome models (`XYmodel`).
```R
dat<-getDat( IVtype='cont', ZXmodel='C',XYmodel='2' )
```

  

Draw the doubly-ranked stratification (use `?Stratify` to check the parameters for stratification)
```R
rdat<-Stratify(dat) 
RES<-getSummaryInf(rdat)
head(RES$DRres)
```
```R
#          size        min        1q      mean        3q        max    Fvalue        bx       bxse         by       byse        est         se target
#Stratum 1 1000 -14.635373 -8.665547 -7.374266 -6.053520 -1.1173800  69.95673 0.9898284 0.11834371 -1.4172479 0.16323707 -1.4318117 0.16491451     NA
#Stratum 2 1000 -12.093102 -6.755072 -5.666040 -4.554825 -1.0455971 159.92760 1.2086197 0.09557140 -1.3118679 0.10666678 -1.0854266 0.08825504     NA
#Stratum 3 1000  -9.489875 -5.611448 -4.575575 -3.524639 -0.3745066 241.20073 1.3153286 0.08469250 -1.2018580 0.08321997 -0.9137321 0.06326934     NA
#Stratum 4 1000  -8.000145 -4.732373 -3.720508 -2.668669 -0.2082876 306.00529 1.4035171 0.08023304 -1.0607026 0.07773271 -0.7557461 0.05538423     NA
#Stratum 5 1000  -7.260873 -3.926916 -2.943308 -1.926862  1.3822956 407.00318 1.5382303 0.07624695 -0.8423861 0.07007546 -0.5476333 0.04555590     NA
#Stratum 6 1000  -6.791493 -3.238730 -2.146673 -1.129673  2.9942412 453.88022 1.6360871 0.07679550 -0.6558742 0.06643903 -0.4008797 0.04060849     NA
```
The function `Stratify()` will return the individual stratification index. The function `getSummaryInf()` will provide summary information based on the stratified result.

Note that the naive and residual stratification can be regarded as special scenarios of the DR stratification. Hence, you may prefer to use only the DR method.

For naive stratification, it is equivalent to applying DR with one pre-stratum:
```R
rdat<-Stratify(dat,SoP=nrow(dat)) 
getSummaryInf( rdat,onlyDR=TRUE)
```
For residual stratification, it is equivalent to applying the naive stratification wrt the residual variables
```R
dat$M<-resid(lm(  dat$X~dat$Z   )) #obtain the residual variables first
rdat<-Stratify(dat,SoP=nrow(dat),onExposure=FALSE) 
getSummaryInf( rdat,onlyDR=TRUE)
```


## Stratify on other covariates
If you need to stratify on a covariate (which is common in many heterogeneous effect studies[^HTE]), first define the covariate to be stratified by as `M` in `dat`, and then simply run the following code.
```R
dat$M<-covarite_vector  #added the covariate information
rdat<-Stratify(dat,onExposure=FALSE) 
getSummaryInf( rdat)
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


## Different types of instruments
The doubly-ranked stratification supports a wide variety of instrument types, including dichotomous IV, discrete IV, high-dimensional IV, continuous IV, mixture IVs. For single IV, just denote its value by `Z` in `dat`. For high-dimensional IV, integrate them as a weighted score (e.g. the weighted gene score in MR). For exmaple
```R
dat$Z<-lm(  dat$X~ dat$G1+ dat$G2 + ...   )$fitted  #G1 G2 ... are genetic variants
```

## Other outcome types
If your outcome is not or cannot be transformed to a continuous outcome with a simple linear model, you may need to use other models[^outcometype]. For example, if your outcome is a binary, you can adjust the fitting model using the argument `family_used`
```R
getSummaryInf( rdat,family_used='binomial')
```
[^outcometype]: Some outcome types are tricky; contact me to discuss the most appropriate model

## Covariate adjustment
If you need to adjust for covariates, add them to your data `dat` or `rdat` with the name `C1`, `C2`, etc. Ensure that the data type is correct as the fitting will follow the rules of the `lm` function[^lm]. To enable covariate adjustment, set the argument `covariate=TRUE`.
```R
getSummaryInf( rdat,covariate=TRUE)
```

[^lm]: All character variables will be treated as factor variables, and corresponding dummy variables will be created to fit in `lm`. If a continuous variable is in character form, it should be transformed into a numerical format before being used in the analysis.

## Control the randomness and reproducibility
One feature of the stratification in the DRMR package is that the stratification function will not introduce randomness (in the sense that the same input data will always generate the same stratification results)[^randomness]. This is beneficial in terms of transparency and reproducibility. However, since ties are broken at random for constant instrument or exposure values, the randomness of such breaks may be considered in some analyses. You can use the argument `seed` to create and track the randomness. For example:  
```R
Stratify(dat)
Stratify(dat,seed=1)
Stratify(dat,seed=2)
 ```
[^randomness]: I considered whether to have the function perform break permutation automatically, but I decided against it as it might result in the user losing control over the variation of results caused by the randomness introduced internally (e.g., if they tend not to set a seed). Although these variations often do not affect any conclusions, I removed the randomness from inside the function to ensure transparency and reproducibility.
       


## Variable transformation
You may consider transformed exposures or covariates[^trans]. The transformed variable could be helpful in terms of interpretation. If you use a bijective transformation, the DR stratification will not be changed. You can simply transform you variable wrt `dat` by using `dat$X<-f(dat$X)` where `f()` is your transformation function`   
[^trans]: One example is the log-transformation: https://www.medrxiv.org/content/10.1101/2022.10.26.22280570v2 
## Smoothing[^smooth]
There are several methods for smoothing the stratification estimation (the stratum-specific estimates are called LACE), including the fractional polynomial method and the piecewise linear method[^smoothing].
[^smooth]: The necessity of smoothing depends on the purpose of your study; if you are doing counterfactual prediction, smoothing is well worth trying; if you are doing inference, smoothing is usually not necessary.

To use the fractional polynomial method (e.g. of degree 2 with first and second order) with the doubly-ranked stratification, run
```R
smooth_res<-Smooth(RES,Norder=3) #RES is the result returned by getSummaryInf()
```
If you wish to visualize the results, simply run `smooth_res$p ` or `smooth_res$hp`
![Rplot](https://user-images.githubusercontent.com/127906571/232224279-908bc7f7-8e38-4871-b4f5-0f19077d569a.png)

![Rplot02](https://user-images.githubusercontent.com/127906571/232224348-0605d2a0-9ff5-4c1a-beb6-c81f8c2c262b.png)




To use the piecewise linear method with the doubly-ranked stratification, run:
```R
cutting_values<-(RES$DRres$mean[-1] + head( RES$DRres$mean,-1) )/2
smooth_res<-Smooth(RES,Norder=1,Knots=cutting_values)
```
[^smoothing]: https://onlinelibrary.wiley.com/doi/full/10.1002/gepi.22041


## Advanced nonlinear effect analysis
We are now able to draw complete stratification-based nonlinear effect analysis, including estimating any complicated (eg, jump) effect shape, and investigating the potential effect change-points/thresholds. Details are given in the new extended package [SSS](https://github.com/HDTian/SSS). We will make this available soon. 


## Real application
The DRMR package follows a one-argument style, which means that in principle all you need to provide is a one-sample data frame. Before using the package, ensure that the input data is in the form of a dataframe:
```R
is.data.frame(your_prepared_dat)
dat<-as.data.frame(your_prepared_dat)
```
The DRMR package does not require you to specify instruments, exposures, or outcomes as arguments. The functions will automatically draw stratification and perform related analysis based on the column names of the inputted dataframe. It is important to ensure that the column names match the required naming conventions for the function to work properly.:
- The instrument: `Z` (see the 'Different types of instrument' section if your instrument is complex)
- The exposure: `X`
- The outcome: `Y`
- (if applicable) the covariate to stratify: `M`
- (if applicable) the covariate to adjust: `C1`, `C2`, `C3`, ...

Check the names of your data frame
```R
colnames(dat)
```

Always be careful about the value type (e.g. for each variable do you need numeric values or character values?)
```R
apply(dat, 2, is.numeric)
apply(dat, 2, is.character)
```
Then clean your data (e.g. values modification, values rescale, imputation, removing missing data, etc)

Assume you now have the samples in a data frame `dat` with correctly named variables (`Z`, `X`, `Y`, etc.), you can simply run the following code:
```R       
rdat<-Stratify(dat)
RES<-getSummaryInf(rdat)
smooth_res<-Smooth(RES,Norder=3)
```      
then you can use the information in `RES` and `smooth_res` (try `?Stratify`, `?getSummaryInf` and `?Smooth` to check the details) to build the results you desire[^further].

[^further]: Feel free to contact me if you have any further questions or need additional information that the current function cannot provide. I would be happy to discuss your study further and provide any assistance you may need.

## FAQs
(Q1) How to selected the appropriate value of "Norder" in "Smooth" function? Should it based on the AIC or BIC of the model (smooth_res)?<br>
(A1) 'Norder' is a value controlling how flexible your fitted curve (ie effect shape) can be. There are many ways to determine it - even relying on your domain knowledge. AIC/BIC given here is also for 'Norder' selection, so you can choose your 'Norder' by, say, choosing the smallest AIC or BIC value.

(Q2) What is the difference between "smooth_res1$hp" and "smooth_res1$p"?<br>
(A2) Both 'smooth_res1$hp' and 'smooth_res1$p' returns you the fitted curve (ie effect shape). $hp returns you the original effect shape, while $p returns you the differentiated effect shape. So technically, $hp is the integral curve of $p. I personally prefer to use $p as the differentiated effect shape plot directly reflects the effect size at the corresponding local exposure value point

(Q3) How should I interpret the result when the exposure is a continuous variable and the outcome is a binary variable?<br>
(A3) When your exposure is continuous and your outcome is binary and your model is the logistic model. The returned LACE estimate is targeted to the log odds ratio - so, say, if you wish to interpret OR (and its 95% CI) you can simply adopt an exponential transformation to them.

(Q4) How to determine if the association between them is linear or non-linear?<br>
(A4) The left upper corner of your 'smooth_res1$p' figure provides two kinds of the tests on the effect shape. The null hypothesis is that the effect is linear, so any rejection (e.g. pvalue<0.05) should give you the evidence to conclude that the effect is nonlinear. You can choose any one of them, but typically, the trend test (at the second line) is of high power and may be recommended.
