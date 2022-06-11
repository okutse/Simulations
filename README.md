---
title: "The Role of Covariate Adjustment in Enhancing Efficiency: A Simulation Approach"
#subtitle: "Homework 1"
author: "Amos Okutse"
date: "  11 June, 2022 "
header-includes:
- \usepackage{pdflscape}
- \newcommand{\blandscape}{\begin{landscape}}
- \newcommand{\elandscape}{\end{landscape}}
- \usepackage{fvextra}
- \usepackage{float}
- \usepackage{wrapfig}
- \usepackage{amsmath}
- \usepackage{float}
- \usepackage{graphicx}
- \usepackage{setspace}
- \usepackage[font=singlespacing]{caption} #can change font here for caption here!!
- \DefineVerbatimEnvironment{Highlighting}{Verbatim}{breaklines, commandchars=\\\{\}}
- \onehalfspacing
fontsize: 10pt
output:
#  bookdown::pdf_document2:
#    latex_engine: xelatex
#    toc: true
#    number_sections: true
#    keep_md: true
#bibliography: references.bib
#csl: APA.csl
#citation_package: natbib
#link-citations: yes
#editor_options: 
#  chunk_output_type: console
  bookdown::html_document2:
    #mathjax: https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-AMS_CHTML.js
    #self_contained: no
    theme: readable
    toc: true
    toc_float: false
    toc_depth: 5
    number_sections: yes
    code_folding: hide #code is hidden by default
    keep_md: true
    highlight: pygments
bibliography: references.bib
csl: APA.csl
citation_package: natbib
link-citations: yes
---







<center>
**Reading Time: 24 minute(s)**
</center>
<br>

## Introduction

In this analysis, we examine the role of covariate adjustment in enhancing efficiency in statistical analysis in Randomized Controlled Trials (RCT) using the case of a linear regression model. In particular, we examine the role that covariate adjustment plays in enhancing efficiency when there is no correlation between the baseline covariates $X's$ and the outcome variable $Y$ and when we have moderate to high correlation between the baseline covariates adjusted for and the outcome variables and the subsequent efficiency improvements on the estimate of the treatment effect in a RCT analyzed using the multiple linear linear regression model. First, the simulation procedures begin by exploring how the estimates of the treatment effect when no baseline covariate adjustment is carried out and how these compare to the estimates of the treatment effect when we adjust for baseline covariates that have a weak to a high correlation with the outcome variable of interest. All simulation procedures in this analyses assume the null hypothesis that baseline covariate adjustment enhances efficiency and precision of the estimate of the treatment effect in a randomized controlled trial and this efficiency increases when the baseline covariates adjusted for in the model have a high association with the outcome variable of interest ($Y$). This is a finding that has been highlighted by @Steingrimsson2017 who note that "In randomized clinical trials with baseline variables that are prognostic for the primary outcome, there is potential to improve precision and reduce sample size by appropriately adjusting for these variables."

## The Data Generating Process

In these analyses, we generate a sample of $n = 500$ observations and simulate $N = 1000$ datasets consisting of four (4) variables, $Y, ~X_1, ~X_2, ~\textrm{and}, ~A$ where $Y$ is the outcome variable, $X_1~ \textrm{and}~ X_2$ are the baseline covariates, and $A$ is the treatment variable denoting the effect of the treatment being studied. We first fix $n = 500$ observations and assign $A$ such that $A = 0$ for the first $\frac{n}{2}$ observations and $A = 1$ for the last $\frac{n}{2}$ observations. We define our outcome model as $Y_i = \alpha + \beta A + \epsilon_i$ where $i$ defines the index of the observation and $\epsilon \sim N(0, ~\sigma^2_{\epsilon})$. The simulated data matrix is such that:

\def\A{
\begin{bmatrix}
    Y\\
    X_1\\
    X_2\\
\end{bmatrix}}

\def\B{
\begin{bmatrix}
\alpha + \beta A \\
0 \\
0\\
\end{bmatrix}
}

\def\C{
\begin{bmatrix}
\sigma^2_e & \tau_1 & \tau_2 \\
\tau_1 & \sigma^2_1 & \sigma^2_2 \\
\tau_2 & \tau_1 & \sigma^2_2\\
\end{bmatrix}
}

\begin{equation}
\label{eqn:datamatrix}
\A ~\sim~ N\left( \B, \C \right) 
\end{equation}

The parameters of the outcome model are pre-specified as $\alpha=0, ~\beta=1, ~\sigma^2_e, ~\sigma^2_j=1$, and $\rho_j=\left(0, ~0.3, ~0.8\right)$. 

## The Simulated Models

We simulate both unadjusted and varied adjusted models to understand the role of adjusting for baseline covariates that are of prognostic value in improving the precision associated with the estimate of the treatment effect using multiple linear regression model, while varying the correlation between these baseline covariates and the outcome variable, $Y$. In this analysis, we simulate the following scenarios:
\begin{equation}
E(Y|A)=\alpha+\beta A
(\#eq:model-one)
\end{equation}

\begin{equation}
E(Y|A, \textbf{X})=\alpha+\beta A + \theta X_1
(\#eq:model-two)
\end{equation}

<!--\begin{equation}
E(Y|A, \textbf{X})=\alpha+\beta A + \theta \frac{X_1}{X_2}
(\#eq:model-three)
\end{equation} -->

\begin{equation}
E(Y|A, \textbf{X})=\alpha+\beta A + \theta_1 X_1+ \theta_2 X_2 + \theta_3 \left(\frac{X_1}{X_2}\right)
(\#eq:model-three)
\end{equation}

\begin{equation}
E(Y|A, \textbf{X})=\alpha+\beta A + \theta_1 \left(\frac{X_1}{X_2}\right) + \theta_2 \left( X_1 - \bar{X_1} \right)^2 + \theta_3 \left(X_2-\bar{X_2}\right)^2
(\#eq:model-four)
\end{equation}

\begin{equation}
E(Y|A, \textbf{X})=\alpha+\beta A + \theta_1 X_1 + \theta_2 X_2+ \theta_3 \frac{X_1}{X_2} + \theta_4 (X_1-\bar{X_1})^2 + \theta_5 (X_2-\bar{X_2})^2
(\#eq:model-five)
\end{equation}


```r
# set the seed for reproducibility
set.seed(1)
```

## Simulations

### Simulation setting 1: a correctly specified un-adjusted linear regression model

In the first simulation setting, we use an un-adjusted linear regression model where we do not adjust for the baseline covariates. The simulated regression model is highlighted in Equation \@ref(eq:model-one). We first create a function `model.one()` which takes in `nobs` as the number of observations, `a` as the pre-specified value of the parameter, $\alpha$, `b` as the pre-specified value of the parameter, $\beta$, and `nreps` as the number of replications (denoting the number of simulated datasets). The `model.one()` function employs the `helper.fun()` function which uses the simulated vector, $\epsilon$ from the normal distribution, the vector, `trt` denoting the treatment variable, the outcome model, $Y$ and returns the estimates of the regression model based on $n = 500$ observations. This helper function is then replicated within `model.one()` to create 1000 simulations of the results spit out from the helper function which are then returned as the average across all simulated datasets using functions provided by the `broom` and `tidyverse` packages. In Table \@ref(tab:sim-one), we highlight the regression estimates based on the un-adjusted linear model and the associated standard errors. In all analyses here, we consider small standard errors to be indicative of an estimate that is very precise (hence a very efficient statistical model) relative to large standard errors.


```r
#setting one:
model.one<-function(nobs=500, a=0, b=1, nreps = 1000){
  
helper.fun=function(){
e=rnorm(nobs)
trt<-c()
for (i in 1:nobs){
  if (i<=nobs/2){
    trt[i]<-0
  }else{
    trt[i]<-1
  }
}
y=a+b*(trt)+e
mod<-lm(y~as.factor(trt))
return(mod)
}

#create 1000 replicates/ simulated datasets
reps=replicate(nreps, helper.fun(), simplify = FALSE)
intercept=reps %>% map_df(tidy) %>% dplyr::filter(term == "(Intercept)")%>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
intercept
trt2<-reps %>% map_df(tidy) %>% dplyr::filter(term == "as.factor(trt)1")%>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#prepare the table of the regression estimates to return
tbl=data.frame(rbind(intercept, trt2))
row.names(tbl)<-c("Intercept", "Treatment1")
return(tbl)
}
```


```r
#with no covariate adjustment
one<-model.one() %>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the unadjusted regression model in simulation setting one", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)
one
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:sim-one)Regression estimates and standard errors from the unadjusted regression model in simulation setting one</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> Std. Error </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Intercept </td>
   <td style="text-align:right;"> 0.0006 </td>
   <td style="text-align:right;"> 0.0632 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Treatment1 </td>
   <td style="text-align:right;"> 0.9978 </td>
   <td style="text-align:right;"> 0.0894 </td>
  </tr>
</tbody>
</table>

Based on the results from the un-adjusted model, we note that the expected value of the outcome is 0.0026 which is very close to our pre-specified value of the parameter $\alpha$. The standard error associated with this estimate is 0.0633. On the other hand, the estimate of the treatment effect is 1.0008 whereas the estimate of the uncertainty associated with this estimate is 0.0895. In the next section, we explore the effect of adjusting for baseline covariates that are (un)related to the outcome variable being investigated and compare the resulting estimate of the treatment effect to that in simulation setting one.

### Simulation setting 2: a correctly specified adjusted linear regression model

In this simulation setting, we use an adjusted model which is also correctly specified as in Equation \@ref(eq:model-two) and we assume three cases. In the first case, we assume that the baseline covariates adjusted for in the model have no correlation with the outcome whereas in the second case, we assume that the baseline covariates are fairly or highly correlated with the outcome variable of interest in our clinical trial. In other words, we assume that the correlation, `rho` is either 0, 0.3 or 0.8 to understand not only the role that covariate adjustment plays in enhancing efficiency (or precision) but also the role that adjusting for covariates that are of prognostic value plays in statistical modeling (here, prognostic baseline covariates are defined as those that have some association with the outcome variable). In the modeling, we use a similar approach as in the first simulation setting but allow for the correlation between the baseline covariates and the outcome to vary by simulating these based on the multivariate normal distribution. 

In particular, we simulate the treatment variable and then generate the variance-covariance matrix such that the diagonal elements denote the variance of the variables whereas the off-diagonal elements denote the correlation between the variables. The outcome variable $Y$ is then simulated conditional on the treatment variable $A$ and the baseline covariate $X_1$ by shifting the mean of the distribution. 1000 datasets are then simulated and a linear regression model fitted to each of these datasets. The estimates of the parameters and the standard errors are then returned by the function as the average of the estimates from the 1000 simulated datasets. 


```r
#setting two:
#set.seed(1)
model.two<-function(nobs = 500, nvars = 2, nreps = 1000, rho = 0, a = 0, b = 1, c = 1){
  helper.fun<-function(){
    
#create the treatment variable where n=250 is trt = 0 and 1 elsewhere
trt<-c()
for (i in 1:nobs){
  if (i<=nobs/2){
    trt[i]<-0
  }else{
    trt[i]<-1
  }
}

#generate the variance-covariance matrix
sigma<-matrix(rho, ncol = nvars, nrow = nvars)

#set the values of the leading diagonal to be 1 for variances
diag(sigma)<-1
mat<-mvrnorm(n = nobs, mu = c(0, 0), Sigma = sigma) # shift the mean of the outcome variable: replace y=z and then y=z+a+b*trt

#assign x1 as the 2nd column and z as col 1
x1 = mat[, 2]; z = mat[, 1]

#define new y conditional on values of trt and x
y = z + a + b*trt + c * x1
mat.new<-data.frame(y, trt, x1)

#run the lm on these data 
lm.mod<-lm(y ~ factor(trt) + x1, data = mat.new)
return(lm.mod)
  }
  
#replicate n=1000 times
sims<-replicate(nreps, helper.fun(), simplify = FALSE )

#average intercept
intercept=sims %>% map_df(tidy) %>% dplyr::filter(term == "(Intercept)") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#average treatment
trt2=sims %>% map_df(tidy) %>% dplyr::filter(term == "factor(trt)1") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#average x1
x=sims %>% map_df(tidy) %>% dplyr::filter(term == "x1") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
rslts=data.frame(rbind(intercept, trt2, x))

#clean up the names for the output
row.names(rslts)<-c("Intercept", "Treatment1", "x")
return(rslts)
}
```

#### Case 1: Adjusting for non-prognostic baseline covariates

Table \@ref(tab:sim-two) presents a summary of the regression estimates based on the adjusted model with no correlation between the adjusted baseline covariates and the outcome variable. This results highlight the implication of adjusting for baseline covariates that are not prognostic.


```r
## Case 1: adjusting for non-prognostic baseline covariates: when correlation between x and y is 0
model.two()%>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting two: adjusting for a non-prognostic baseline covariate", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:sim-two)Regression estimates and standard errors from the adjusted regression model in simulation setting two: adjusting for a non-prognostic baseline covariate</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> Std. Error </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Intercept </td>
   <td style="text-align:right;"> 0.0005 </td>
   <td style="text-align:right;"> 0.0633 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Treatment1 </td>
   <td style="text-align:right;"> 0.9963 </td>
   <td style="text-align:right;"> 0.0895 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x </td>
   <td style="text-align:right;"> 1.0019 </td>
   <td style="text-align:right;"> 0.0448 </td>
  </tr>
</tbody>
</table>

The estimate of the treatment effect based on the adjusted model when $r = 0$ is 0.9956 and the corresponding standard error is 0.0897. We notice a very substantial correspondence between the estimate of the treatment effect and the associated estimate of the standard error between the un-adjusted model and the model adjusting for baseline covariates unrelated to the outcome variable. As expected, we notice that adjusting for variables that are not prognostic has the effect of resulting in salient increases in the standard errors and hence a decrease in power. A similar observation was highlighted elsewhere [@Kahan2014]. In their paper, @Kahan2014 noted that adjusting for baseline covariates that are not prognostic in RCT might have the effect of reducing power since:

>"... each continuous or binary baseline covariate uses a ‘degree of freedom’, which effectively reduces the sample size, meaning that there is less information with which to estimate the treatment effect (in cases where the covariate actually is prognostic, the benefits of the prognostic ability outweigh any loss of information, and power will be increased despite the loss of a degree of freedom)."

#### Case 2: Adjusting for fairly prognostic baseline covariates

In Table \@ref(tab:model-twob) we highlight the results from the regression model where we adjust for a baseline covariate, $X_1$ which is fairly prognostic (that is the case where we have a fair association ($r = 0.3$) between this variable and the outcome variable. 


```r
## Case 2: adjusting for a fairly prognostic variable (r=0.3)
model.two(rho = 0.3)%>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting two: adjusting for a fairly prognostic variable", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:model-twob)Regression estimates and standard errors from the adjusted regression model in simulation setting two: adjusting for a fairly prognostic variable</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> Std. Error </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Intercept </td>
   <td style="text-align:right;"> 0.0007 </td>
   <td style="text-align:right;"> 0.0603 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Treatment1 </td>
   <td style="text-align:right;"> 0.9980 </td>
   <td style="text-align:right;"> 0.0852 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x </td>
   <td style="text-align:right;"> 1.2997 </td>
   <td style="text-align:right;"> 0.0427 </td>
  </tr>
</tbody>
</table>

We note that when we adjust for a baseline covariate that is fairly prognostic, we experience prominent gains in the power associated with the estimate of the treatment effect. In particular, we note that while the estimate of the treatment effect is 0.999 (SE = 0.0895) based on the un-adjusted model, the estimate of the treatment effect based on the model adjusting for a fairly prognostic variable is 1.0003 (SE = 0.0895). 

#### Case 3: Adjusting for (a highly) prognostic baseline covariates

Table \@ref(tab:model-twoc) highlights the regression estimates and standard errors from the adjusted regression model in simulation setting two: adjusting for a prognostic variable. We simulated the baseline covariate adjusted for, such that it is very highly associated with the baseline covariate ($r = 0.8$). 


```r
## Case 3: adjusting for a prognostic variables (r=0.8)
model.twoc = model.two(rho = 0.8)
model.twoc%>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting two: adjusting for a prognostic variable", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:model-twoc)Regression estimates and standard errors from the adjusted regression model in simulation setting two: adjusting for a prognostic variable</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> Std. Error </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Intercept </td>
   <td style="text-align:right;"> -0.0004 </td>
   <td style="text-align:right;"> 0.0379 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Treatment1 </td>
   <td style="text-align:right;"> 1.0016 </td>
   <td style="text-align:right;"> 0.0536 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x </td>
   <td style="text-align:right;"> 1.8003 </td>
   <td style="text-align:right;"> 0.0269 </td>
  </tr>
</tbody>
</table>

Table \@ref(tab:model-twoc) highlights the estimates based on the regression model adjusting for prognostic baseline covariates. In particular, we note that while the estimate of the treatment effect is 1.0011 (SE = 0.0895) based on the un-adjusted model, the estimate of the treatment effect based on the model adjusting for a fairly prognostic variable is 1.0006 (SE = 0.0894) whereas the estimate of the treatment effect based on the model adjusting for a prognostic baseline variable is 1.0016 (SE = 0.0536). Notice the reduction in the estimate of the standard error alluding to the very substantial gain in power we experience when we adjust for prognostic baseline covariates from case 2 relative to case 1 in this second simulation setting. This gains have been reported in length elsewhere [@Kahan2014;@Steingrimsson2017]. In the subsequent simulation settings, we explore the implication of adjusting for both prognostic and non-prognostic baseline covariates when we have a mis-specified statistical model. 

### Simulation setting 3: a mis-specified adjusted linear regression model

In this simulation we explore the implication of both prognostic and non-prognostic baseline covariate adjustment in the context of an incorrectly specified linear regression model. This setting is based on Equation \@ref(eq:model-three) and the data generating process described previously is used. We explore the implication of adjusting for (non)prognostic baseline covariates in a RCT by allowing the correlation between the baseline covariates and the outcome variable to vary at different levels. To ensure that our model does not suffer from multicollinearity, we set the correlation between the baseline covariates to 0 in the variance-covariance matrix simulated from the multivariate normal distribution and analyze three cases as in simulation setting 2. In the first case, we adjust for baseline covariates that are not prognostic whereas in the second and third cases we adjust for baseline covariates that are fairly and highly prognostic, respectively.  


```r
#setting 3: adjusting for x1, x2, and x1/x2
#set.seed(1)

model.three<-function(nobs = 500, nvars = 3, nreps = 1000, rho = 0, a = 0, b = 1, theta1 = 1, theta2 = 1, theta3 = 1){
  helper.fun<-function(){
#create the treatment variable where n=250 is treatment = 0 and 1 elsewhere
trt<-c()
for (i in 1:nobs){
  if (i<=nobs/2){
    trt[i]<-0
  }else{
    trt[i]<-1
  }
}

#generate the variance-covariance matrix
sigma<-matrix(rho, ncol = nvars, nrow = nvars)

#set the values of the leading diagonal to be 1 for variances
diag(sigma)<-1
sigma[2,3] = 0; sigma[3,2] = 0
mat<-mvrnorm(n = nobs, mu = c(0, 0, 0), Sigma = sigma, tol = 0.1) #noted idea;replace y=z and then y=z+a+b*trt

#assign x1 as the 2nd column and z as col 1
x1 = mat[, 2]; z = mat[, 1]; x2 = mat[, 3]; x = (x1/x2)

#define new y conditional on values of treatment and X's
y = z + a + b*trt + theta1*x1 + theta2*x2 + theta3 * x
mat.new<-data.frame(y, trt, x, x1, x2)
#run the lm on these data 
lm.mod<-lm(y ~ factor(trt) + x1 + x2 + x, data = mat.new)
#print(lm.mod)
return(lm.mod)
  }
  
#replicate n=1000 times
sims<-replicate(nreps, helper.fun(), simplify = FALSE )

#average intercept
intercept=sims %>% map_df(tidy) %>% dplyr::filter(term == "(Intercept)") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#average treatment
trt2=sims %>% map_df(tidy) %>% dplyr::filter(term == "factor(trt)1") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
#average x1
x1=sims %>% map_df(tidy) %>% dplyr::filter(term == "x1") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
#average x2
x2=sims %>% map_df(tidy) %>% dplyr::filter(term == "x2") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
#average x=x1/x2
x=sims %>% map_df(tidy) %>% dplyr::filter(term == "x") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
rslts=data.frame(rbind(intercept, trt2, x1, x2, x))
row.names(rslts)<-c("Intercept", "Treatment2", "x1", "x2", "x=(x1/x2)")
return(rslts)
}
```

#### Case 1: Adjusting for non-prognostic baseline covariates in a mis-specified model with main effects

Table \@ref(tab:sim-three) presents a summary of the regression estimates based on the adjusted model with no correlation between the adjusted baseline covariates and the outcome variable. This results highlight the implication of adjusting for baseline covariates that are not prognostic.


```r
## Case 1: adjusting for non-prognostic baseline covariates: when correlation between x and y is 0
model.three()%>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting three: adjusting for a non-prognostic baseline covariate", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE) %>% 
  add_footnote("These results are based on an incorrect model specification")
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:sim-three)Regression estimates and standard errors from the adjusted regression model in simulation setting three: adjusting for a non-prognostic baseline covariate</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> Std. Error </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Intercept </td>
   <td style="text-align:right;"> -0.0028 </td>
   <td style="text-align:right;"> 0.0634 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Treatment2 </td>
   <td style="text-align:right;"> 1.0008 </td>
   <td style="text-align:right;"> 0.0897 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x1 </td>
   <td style="text-align:right;"> 0.9987 </td>
   <td style="text-align:right;"> 0.0450 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x2 </td>
   <td style="text-align:right;"> 1.0022 </td>
   <td style="text-align:right;"> 0.0448 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x=(x1/x2) </td>
   <td style="text-align:right;"> 0.9998 </td>
   <td style="text-align:right;"> 0.0020 </td>
  </tr>
</tbody>
<tfoot>
<tr>
<td style = 'padding: 0; border:0;' colspan='100%'><sup>a</sup> These results are based on an incorrect model specification</td>
</tr>
</tfoot>
</table>

The estimate of the treatment effect based on the adjusted mis-specified model (Equation \@ref(eq:model-three)) when $r = 0$ is 0.9974 and the corresponding standard error is 0.0897 whereas the estimate of the treatment effect is 0.9979 (SE = 0.0892) based on the un-adjusted model (Equation \@ref(eq:model-one)). We, similarly, notice a very substantial correspondence between the estimate of the treatment effect and the associated estimate of the standard error between the un-adjusted model and the model adjusting for baseline covariates unrelated to the outcome variable. Similarly, we observe that adjusting for variables that are not prognostic has the effect of resulting in salient increases in the standard errors and hence a decrease in power.

#### Case 2: Adjusting for fairly prognostic baseline covariates in a mis-specified model with main effects

In Table \@ref(tab:model-threeb) we highlight the results from the regression model where we adjust for baseline covariates, which are fairly prognostic (that is the case where we have a fair association ($r = 0.3$) between these variables and the outcome variable. 


```r
## Case 2: adjusting for a fairly prognostic variable (r=0.3) in a mis-specified model
model.threeb = model.three(rho = 0.3)
model.threeb %>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting three: adjusting for a fairly prognostic variable", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)%>% 
  add_footnote("These results are based on an incorrect model specification")
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:model-threeb)Regression estimates and standard errors from the adjusted regression model in simulation setting three: adjusting for a fairly prognostic variable</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> Std. Error </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Intercept </td>
   <td style="text-align:right;"> 0.0023 </td>
   <td style="text-align:right;"> 0.0573 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Treatment2 </td>
   <td style="text-align:right;"> 0.9985 </td>
   <td style="text-align:right;"> 0.0811 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x1 </td>
   <td style="text-align:right;"> 1.3001 </td>
   <td style="text-align:right;"> 0.0407 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x2 </td>
   <td style="text-align:right;"> 1.3008 </td>
   <td style="text-align:right;"> 0.0406 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x=(x1/x2) </td>
   <td style="text-align:right;"> 1.0000 </td>
   <td style="text-align:right;"> 0.0018 </td>
  </tr>
</tbody>
<tfoot>
<tr>
<td style = 'padding: 0; border:0;' colspan='100%'><sup>a</sup> These results are based on an incorrect model specification</td>
</tr>
</tfoot>
</table>

The estimate of the treatment effect is 1.0026 (SE = 0.0893) based on the correctly specified un-adjusted model whereas the estimate of the treatment effect based on the mis-specified model adjusting for fairly prognostic variables is 0.9985 (SE = 0.0811). 

#### Case 3: Adjusting for (a highly) prognostic baseline covariates in a mis-specified model with main effects

Table \@ref(tab:model-threec) highlights the regression estimates and standard errors from the adjusted regression model in simulation setting three: adjusting for prognostic variables. We simulated the baseline covariates adjusted for, such that they are very highly associated with the baseline covariates ($r = 0.8$). 


```r
## Case 3: adjusting for a prognostic variables (r=0.8)
model.threec = model.two(rho = 0.8)
model.threec%>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting three: adjusting for a prognostic variable", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)%>% 
  add_footnote("These results are based on an incorrect model specification")
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:model-threec)Regression estimates and standard errors from the adjusted regression model in simulation setting three: adjusting for a prognostic variable</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> Std. Error </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Intercept </td>
   <td style="text-align:right;"> 0.0006 </td>
   <td style="text-align:right;"> 0.0380 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Treatment1 </td>
   <td style="text-align:right;"> 0.9988 </td>
   <td style="text-align:right;"> 0.0537 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x </td>
   <td style="text-align:right;"> 1.7995 </td>
   <td style="text-align:right;"> 0.0270 </td>
  </tr>
</tbody>
<tfoot>
<tr>
<td style = 'padding: 0; border:0;' colspan='100%'><sup>a</sup> These results are based on an incorrect model specification</td>
</tr>
</tfoot>
</table>

Table \@ref(tab:model-threec) highlights the estimates based on the regression model adjusting for prognostic baseline covariates. In particular, we note that while the estimate of the treatment effect is 1 (SE = 0.0893) based on the un-adjusted model, the estimate of the treatment effect based on the model adjusting for a fairly prognostic variable is 0.9985 (SE = 0.0811) whereas the estimate of the treatment effect based on the model adjusting for a strong prognostic baseline variable is 0.9988 (SE = 0.0537). We similarly notice the reduction in the estimate of the standard error alluding to the very substantial gain in power we experience when we adjust for prognostic baseline covariates from case 2 relative to case 1 in this third simulation setting even when the model specification is wrong. 

### Simulation setting 4: a mis-specified adjusted linear regression model

The outcome model in simulation setting five is as highlighted in Equation \@ref(eq:model-four) where we adjust for $x_1/x_2$, $(x_1-\bar{x}_1)^2$, and $(x_2-\bar{x_2})^2$ which are functions of the prognostic baseline covariates $X_1$ and $X_2$. Note that the model excludes the main effects of the prognostic baseline covariates. In this simulation setting, we under-specify the model such that the implication of covariate adjustment in such a case is explored. In particular, we expect that the results should conform to those obtained in case 1 from both simulation setting one and two where we have adjusted for baseline covariates unrelated to the outcome variable since the variables adjusted for are not directly related to the outcome variable, $Y$. The estimates of the standard error on the treatment variable from this model should highlight reductions in power and hence reductions in efficiency. Even though we vary the level of correlation between the main covariates ($X_1$ and $X_2$) we expect these correlations to have minimal effects (if any) on the covariates that we adjust for since these are extrapolations of the variables that are directly related to the outcome. In the first case, we adjust for baseline covariates that are not prognostic whereas in the second and third cases we adjust for baseline covariates that are fairly and highly prognostic, respectively. 


```r
#Scenario 4: mis-specified model adjusting for x1/x2, (x1-x1bar)^2, (x2-x2bar)^2
model.four<-function(nobs = 500, nvars = 3, nreps = 1000, rho = 0, a = 0, b = 1, theta1 = 1, theta2 = 1, theta3 = 1){
  helper.fun<-function(){
#create the treatment variable where n=250 is trt = 0 and 1 elsewhere
trt<-c()
for (i in 1:nobs){
  if (i<=nobs/2){
    trt[i]<-0
  }else{
    trt[i]<-1
  }
}

#generate the variance-covariance matrix
sigma<-matrix(rho, ncol = nvars, nrow = nvars)

#set the values of the leading diagonal to be 1 for variances
diag(sigma)<-1
sigma[2,3] = 0; sigma[3,2] = 0
mat<-mvrnorm(n = nobs, mu = c(0, 0, 0), Sigma = sigma, tol = 0.1) #noted idea;replace y=z and then y=z+a+b*trt

#assign x1 as the 2nd column and z as col 1
x1 = mat[, 2]; z = mat[, 1]; x2 = mat[, 3]; x = (x1/x2)
x1dev=(x1-mean(x1))^2; x2dev=(x2-mean(x2))^2

#define new y conditional on values of trt and x
y = z + a + b*trt + theta1*x + theta2*x1dev +theta3*x2dev
mat.new<-data.frame(y, trt, x, x1dev, x2dev)
#run the lm on these data 
lm.mod<-lm(y ~ factor(trt) + x + x1dev + x2dev, data = mat.new)
#print(lm.mod)
return(lm.mod)
  }
  
#replicate n=1000 times
sims<-replicate(nreps, helper.fun(), simplify = FALSE )

#average intercept
intercept=sims %>% map_df(tidy) %>% dplyr::filter(term == "(Intercept)") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#average treatment
trt2=sims %>% map_df(tidy) %>% dplyr::filter(term == "factor(trt)1") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#average x=x1/x2
x=sims %>% map_df(tidy) %>% dplyr::filter(term == "x") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
#average x1dev
x1dev=sims %>% map_df(tidy) %>% dplyr::filter(term == "x1dev") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
#average x2dev
x2dev=sims %>% map_df(tidy) %>% dplyr::filter(term == "x2dev") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
rslts=data.frame(rbind(intercept, trt2, x, x1dev, x2dev))
row.names(rslts)<-c("Intercept", "trt2", "x1/x2", "(x1-x1bar)^2", "(x2-x2bar)^2")
return(rslts)
}
```

#### Case 1: Adjusting for non-prognostic baseline covariates in a mis-specified model with no main effects

Table \@ref(tab:sim-three) presents a summary of the regression estimates based on the adjusted model with no correlation between the adjusted baseline covariates and the outcome variable. This results highlight the implication of adjusting for baseline covariates that are not prognostic in an under-specified model.


```r
## Case 1: adjusting for non-prognostic baseline covariates: when correlation between x and y is 0
model.four()%>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting four: adjusting for non-prognostic baseline covariates", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE) %>% 
  add_footnote("These results are based on an incorrect model specification. The model is under-specified since the main effects of the prognostic baseline variables are not included.")
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:sim-four)Regression estimates and standard errors from the adjusted regression model in simulation setting four: adjusting for non-prognostic baseline covariates</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> Std. Error </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Intercept </td>
   <td style="text-align:right;"> -0.0053 </td>
   <td style="text-align:right;"> 0.0778 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> trt2 </td>
   <td style="text-align:right;"> 1.0031 </td>
   <td style="text-align:right;"> 0.0897 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x1/x2 </td>
   <td style="text-align:right;"> 1.0000 </td>
   <td style="text-align:right;"> 0.0020 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (x1-x1bar)^2 </td>
   <td style="text-align:right;"> 1.0003 </td>
   <td style="text-align:right;"> 0.0322 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (x2-x2bar)^2 </td>
   <td style="text-align:right;"> 1.0015 </td>
   <td style="text-align:right;"> 0.0321 </td>
  </tr>
</tbody>
<tfoot>
<tr>
<td style = 'padding: 0; border:0;' colspan='100%'><sup>a</sup> These results are based on an incorrect model specification. The model is under-specified since the main effects of the prognostic baseline variables are not included.</td>
</tr>
</tfoot>
</table>

The estimate of the treatment effect based on the adjusted mis-specified model (Equation \@ref(eq:model-four)) when $r = 0$ is 1.0028 and the corresponding standard error is 0.0896 whereas the estimate of the treatment effect is 0.9977 (SE = 0.0894) based on the un-adjusted model (Equation \@ref(eq:model-one)). We, similarly, notice a very substantial correspondence between the estimate of the treatment effect and the associated estimate of the standard error between the un-adjusted model and the model adjusting for baseline covariates unrelated to the outcome variable. Similarly, adjusting for variables that are not prognostic has the effect of resulting in salient increases in the standard errors and hence a decrease in power.

#### Case 2: Adjusting for fairly prognostic baseline covariates in a mis-specified model with no main effects

In Table \@ref(tab:model-fourbb) we highlight the results from the regression model where we adjust for baseline covariates, which are fairly prognostic (that is the case where we have a fair association ($r = 0.3$) between these variables and the outcome variable. 


```r
## Case 2: adjusting for a fairly prognostic variable (r=0.3) in a mis-specified model
model.fourb = model.four(rho = 0.3)
model.fourb %>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting four: adjusting for a fairly prognostic variable", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)%>% 
  add_footnote("These results are based on an incorrect model specification. The model is under-specified since the main effects of the prognostic baseline variables are not included.")
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:model-fourbb)Regression estimates and standard errors from the adjusted regression model in simulation setting four: adjusting for a fairly prognostic variable</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> Std. Error </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Intercept </td>
   <td style="text-align:right;"> -0.0018 </td>
   <td style="text-align:right;"> 0.0778 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> trt2 </td>
   <td style="text-align:right;"> 1.0025 </td>
   <td style="text-align:right;"> 0.0896 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x1/x2 </td>
   <td style="text-align:right;"> 1.0000 </td>
   <td style="text-align:right;"> 0.0020 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (x1-x1bar)^2 </td>
   <td style="text-align:right;"> 1.0012 </td>
   <td style="text-align:right;"> 0.0321 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (x2-x2bar)^2 </td>
   <td style="text-align:right;"> 0.9996 </td>
   <td style="text-align:right;"> 0.0321 </td>
  </tr>
</tbody>
<tfoot>
<tr>
<td style = 'padding: 0; border:0;' colspan='100%'><sup>a</sup> These results are based on an incorrect model specification. The model is under-specified since the main effects of the prognostic baseline variables are not included.</td>
</tr>
</tfoot>
</table>

The estimate of the treatment effect is 1.0007 (SE = 0.0895) based on the correctly specified un-adjusted model whereas the estimate of the treatment effect based on the mis-specified model adjusting for fairly prognostic variables is 1.0025 (SE = 0.0896). 

#### Case 3: Adjusting for (a highly) prognostic baseline covariates in a mis-specified model with no main effects

Table \@ref(tab:model-fourcc) highlights the regression estimates and standard errors from the adjusted regression model in simulation setting three: adjusting for prognostic variables. Baseline covariates adjusted for are simulated such that they are very highly associated with the baseline covariates ($r = 0.8$). 


```r
## Case 3: adjusting for a prognostic variables (r=0.8)
model.fourc = model.four(rho = 0.8)
model.fourc%>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting four: adjusting for a prognostic variable", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)%>% 
  add_footnote("These results are based on an incorrect model specification. The model is under-specified since the main effects of the prognostic baseline variables are not included.")
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:model-fourcc)Regression estimates and standard errors from the adjusted regression model in simulation setting four: adjusting for a prognostic variable</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> Std. Error </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Intercept </td>
   <td style="text-align:right;"> 0.0008 </td>
   <td style="text-align:right;"> 0.0802 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> trt2 </td>
   <td style="text-align:right;"> 0.9984 </td>
   <td style="text-align:right;"> 0.0924 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x1/x2 </td>
   <td style="text-align:right;"> 1.0000 </td>
   <td style="text-align:right;"> 0.0021 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (x1-x1bar)^2 </td>
   <td style="text-align:right;"> 0.9991 </td>
   <td style="text-align:right;"> 0.0319 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (x2-x2bar)^2 </td>
   <td style="text-align:right;"> 1.0014 </td>
   <td style="text-align:right;"> 0.0320 </td>
  </tr>
</tbody>
<tfoot>
<tr>
<td style = 'padding: 0; border:0;' colspan='100%'><sup>a</sup> These results are based on an incorrect model specification. The model is under-specified since the main effects of the prognostic baseline variables are not included.</td>
</tr>
</tfoot>
</table>

Table \@ref(tab:model-fourcc) highlights the estimates based on the regression model adjusting for prognostic baseline covariates in an under-specified model (Equation \@ref(eq:model-four)). We note that while the estimate of the treatment effect is 0.9974 (SE = 0.0894) based on the un-adjusted model, the estimate of the treatment effect based on the model adjusting for a fairly prognostic variable is 1.0025 (SE = 0.0896) whereas the estimate of the treatment effect based on the model adjusting for a strong prognostic baseline variable is 0.9984 (SE = 0.0924). *We observe that adjusting for baseline covariates not directly related to the outcome of interest in an unspecified model has the effect of decreasing power.*  

### Simulation setting 5: a mis-specified adjusted linear regression model

The outcome model in simulation setting five is as highlighted in Equation \@ref(eq:model-five) where we adjust for $X_1$, $X_2$, $X_1/X_2$, $(X_1-\bar{X}_1)^2$, and $(X_2-\bar{X_2})^2$. In this simulation setting, we adjust for both baseline covariates, as well as, $X_1/X_2$, $(X_1-\bar{X}_1)^2$, and $(X_2-\bar{X_2})^2$ which are derived from $X_1$ and $X_2$. In this simulation setting we notice that even though we adjust for baseline covariates for which we can control how much association they have with the outcome variable, we also adjust for covariates derived from $X_1$ and $X_2$ and which have unknown distributions and which are not directly related to the outcome variable. Similarly and as in the other simulation settings, we control for the correlation between the main baseline covariates ($X_1$ and $X_2$) begining with case 1 where we assume a no association, the second case where we assume a fair correlation, and the third case where we assume that these two covariates are strongly prognostic. We expect that including the main effects of the baseline covariates should result in increases in power with increasing correlation.


```r
# Scenario 5: Adjust for x1/x2, (x1-x1bar)^2, (x2-x2bar)^2, x1, x2
model.five<-function(nobs = 500, nvars = 3, nreps = 1000, rho = 0, a = 0, b = 1, theta1 = 1, theta2 = 1, theta3 = 1, theta4 =1, theta5 = 1){
helper.fun<-function(){
#create the treatment variable where n=250 is trt = 0 and 1 elsewhere
trt<-c()
for (i in 1:nobs){
  if (i<=nobs/2){
    trt[i]<-0
  }else{
    trt[i]<-1
  }
}

#generate the variance-covariance matrix
sigma<-matrix(rho, ncol = nvars, nrow = nvars)

#set the values of the leading diagonal to be 1 for variances
diag(sigma)<-1
sigma[2,3] = 0; sigma[3,2] = 0
mat<-mvrnorm(n = nobs, mu = c(0, 0, 0), Sigma = sigma, tol = 0.1) 
#assign x1 as the 2nd column and z as col 1
x1 = mat[, 2]; z = mat[, 1]; x2 = mat[, 3]; x = (x1/x2)
x1dev=(x1-mean(x1))^2; x2dev=(x2-mean(x2))^2

#define new y conditional on values of trt and x
y = z + a + b*trt + theta1*x1 + theta2*x2 + theta3*x + theta4*x1dev +theta5*x2dev
mat.new<-data.frame(y, trt, x1, x2, x, x1dev, x2dev)
#run the lm on these data 
lm.mod<-lm(y ~ factor(trt) + x1 + x2 + x + x1dev + x2dev, data = mat.new)
#print(lm.mod)
return(lm.mod)
  }
  
#replicate n=1000 times
sims<-replicate(nreps, helper.fun(), simplify = FALSE )

#average intercept
intercept=sims %>% map_df(tidy) %>% dplyr::filter(term == "(Intercept)") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#average treatment
trt2=sims %>% map_df(tidy) %>% dplyr::filter(term == "factor(trt)1") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#average x1
x1=sims %>% map_df(tidy) %>% dplyr::filter(term == "x1") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#average x2
x2=sims %>% map_df(tidy) %>% dplyr::filter(term == "x2") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#average x=x1/x2
x=sims %>% map_df(tidy) %>% dplyr::filter(term == "x") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
#average x1dev
x1dev=sims %>% map_df(tidy) %>% dplyr::filter(term == "x1dev") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
#average x2dev
x2dev=sims %>% map_df(tidy) %>% dplyr::filter(term == "x2dev") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
rslts=data.frame(rbind(intercept, trt2, x1, x2, x, x1dev, x2dev))
row.names(rslts)<-c("Intercept", "trt2", "x1", "x2", "x1/x2", "(x1-x1bar)^2", "(x2-x2bar)^2")
return(rslts)
}
```

#### Case 1: Adjusting for non-prognostic baseline covariates in a mis-specified model with no main effects

Table \@ref(tab:sim-five) presents a summary of the regression estimates based on the adjusted model with no correlation between the adjusted baseline covariates and the outcome variable. This results highlight the implication of adjusting for baseline covariates that are not prognostic in a mis-specified model.


```r
## Case 1: adjusting for non-prognostic baseline covariates: when correlation between x and y is 0
model.five()%>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting five: adjusting for non-prognostic baseline covariates", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE) %>% 
  add_footnote("These results are based on an incorrect model specification. The model is specified to include the main effects of the prognostic baseline variables.")
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:sim-five)Regression estimates and standard errors from the adjusted regression model in simulation setting five: adjusting for non-prognostic baseline covariates</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> Std. Error </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Intercept </td>
   <td style="text-align:right;"> 0.0026 </td>
   <td style="text-align:right;"> 0.0781 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> trt2 </td>
   <td style="text-align:right;"> 0.9993 </td>
   <td style="text-align:right;"> 0.0899 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x1 </td>
   <td style="text-align:right;"> 0.9999 </td>
   <td style="text-align:right;"> 0.0452 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x2 </td>
   <td style="text-align:right;"> 0.9990 </td>
   <td style="text-align:right;"> 0.0450 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x1/x2 </td>
   <td style="text-align:right;"> 1.0001 </td>
   <td style="text-align:right;"> 0.0020 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (x1-x1bar)^2 </td>
   <td style="text-align:right;"> 0.9993 </td>
   <td style="text-align:right;"> 0.0322 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (x2-x2bar)^2 </td>
   <td style="text-align:right;"> 1.0009 </td>
   <td style="text-align:right;"> 0.0321 </td>
  </tr>
</tbody>
<tfoot>
<tr>
<td style = 'padding: 0; border:0;' colspan='100%'><sup>a</sup> These results are based on an incorrect model specification. The model is specified to include the main effects of the prognostic baseline variables.</td>
</tr>
</tfoot>
</table>

The estimate of the treatment effect based on the adjusted mis-specified model (Equation \@ref(eq:model-five)) when $r = 0$ is 1.0042 and the corresponding standard error is 0.0898 whereas the estimate of the treatment effect is 0.9985 (SE = 0.0894) based on the un-adjusted model (Equation \@ref(eq:model-one)). We, similarly, notice a very substantial correspondence between the estimate of the treatment effect and the associated estimate of the standard error between the un-adjusted model and the model adjusting for baseline covariates unrelated to the outcome variable. Additionally, adjusting for variables that are not prognostic has the effect of resulting in salient increases in the standard errors and hence a decrease in power.

#### Case 2: Adjusting for fairly prognostic baseline covariates in a mis-specified model with no main effects

In Table \@ref(tab:model-fivebb) we highlight the results from the regression model where we adjust for baseline covariates, which are fairly prognostic (that is the case where we have a fair association ($r = 0.3$) between these variables and the outcome variable. 


```r
## Case 2: adjusting for a fairly prognostic variable (r=0.3) in a mis-specified model
model.fiveb = model.five(rho = 0.3)
model.fiveb %>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting five: adjusting for a fairly prognostic variable", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)%>% 
  add_footnote("These results are based on an incorrect model specification. The model is specified to include the main effects of the prognostic baseline variables.")
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:model-fivebb)Regression estimates and standard errors from the adjusted regression model in simulation setting five: adjusting for a fairly prognostic variable</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> Std. Error </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Intercept </td>
   <td style="text-align:right;"> -0.0008 </td>
   <td style="text-align:right;"> 0.0706 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> trt2 </td>
   <td style="text-align:right;"> 1.0012 </td>
   <td style="text-align:right;"> 0.0813 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x1 </td>
   <td style="text-align:right;"> 1.3004 </td>
   <td style="text-align:right;"> 0.0409 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x2 </td>
   <td style="text-align:right;"> 1.2994 </td>
   <td style="text-align:right;"> 0.0408 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x1/x2 </td>
   <td style="text-align:right;"> 1.0001 </td>
   <td style="text-align:right;"> 0.0019 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (x1-x1bar)^2 </td>
   <td style="text-align:right;"> 0.9992 </td>
   <td style="text-align:right;"> 0.0293 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (x2-x2bar)^2 </td>
   <td style="text-align:right;"> 1.0009 </td>
   <td style="text-align:right;"> 0.0293 </td>
  </tr>
</tbody>
<tfoot>
<tr>
<td style = 'padding: 0; border:0;' colspan='100%'><sup>a</sup> These results are based on an incorrect model specification. The model is specified to include the main effects of the prognostic baseline variables.</td>
</tr>
</tfoot>
</table>

The estimate of the treatment effect is 1.0017 (SE = 0.0892) based on the correctly specified un-adjusted model whereas the estimate of the treatment effect based on the mis-specified model adjusting for fairly prognostic variables is 1.0012 (SE = 0.0813). 

#### Case 3: Adjusting for (a highly) prognostic baseline covariates in a mis-specified model with no main effects

Table \@ref(tab:model-fivecc) highlights the regression estimates and standard errors from the adjusted regression model in simulation setting three: adjusting for prognostic variables. Baseline covariates adjusted for are simulated such that they are very highly associated with the baseline covariates ($r = 0.8$). 


```r
## Case 3: adjusting for a prognostic variables (r=0.8)
model.fivec = model.five(rho = 0.8)
model.fivec %>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting five: adjusting for a prognostic variable", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)%>% 
  add_footnote("These results are based on an incorrect model specification. The model is specified to include the main effects of the prognostic baseline variables.")
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:model-fivecc)Regression estimates and standard errors from the adjusted regression model in simulation setting five: adjusting for a prognostic variable</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> Std. Error </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Intercept </td>
   <td style="text-align:right;"> 0.0000 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> trt2 </td>
   <td style="text-align:right;"> 1.0000 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x1 </td>
   <td style="text-align:right;"> 1.7071 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x2 </td>
   <td style="text-align:right;"> 1.7071 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x1/x2 </td>
   <td style="text-align:right;"> 1.0000 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (x1-x1bar)^2 </td>
   <td style="text-align:right;"> 1.0000 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (x2-x2bar)^2 </td>
   <td style="text-align:right;"> 1.0000 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
<tfoot>
<tr>
<td style = 'padding: 0; border:0;' colspan='100%'><sup>a</sup> These results are based on an incorrect model specification. The model is specified to include the main effects of the prognostic baseline variables.</td>
</tr>
</tfoot>
</table>

Table \@ref(tab:model-fivecc) highlights the estimates based on the regression model adjusting for prognostic baseline covariates in an under-specified model (Equation \@ref(eq:model-five)). We note that while the estimate of the treatment effect is 0.9978 (SE = 0.0895) based on the un-adjusted model, the estimate of the treatment effect based on the model adjusting for a fairly prognostic variable is 1.0012 (SE = 0.0813) whereas the estimate of the treatment effect based on the model adjusting for a strong prognostic baseline variable is 1 (SE = 0). *We observe that adjusting for baseline covariates not directly related to the outcome of interest in an unspecified model has the effect of increasing power especially when the model also includes the main effects of the baseline covariates directly related to the outcome variable $X_1$ and $X_2$.*  

## Conclusions 

In this post, we have explored the implications of adjusting for both prognostic and non-prognostic baseline covariates in RCTs using a linear regression model. We performed simulations in five different settings exploring the impact of adjusting for both prognostic and non-prognostic baseline covariates when we have correctly and incorrectly specified statistical modeling. In the first simulation setting we have a correctly specified model with no covariate adjustment whereas the second simulation setting considers a correctly specified linear regression model that also adjusts for both prognostic and non-prognostic baseline covariates. We notice that the standard error on the treatment variable is much higher in the un-adjusted model relative to the adjusted model. In particular, adjusting for non-prognostic baseline covariates has the effect of resulting in a higher SE alluding to a reduction in power. On the other hand, when the baseline covariates adjusted for are prognostic we experience salient gains in statistical power as seen through reductions in the SE and these efficiency gains increase with increasing correlation between the baseline covariates and the outcome variable.

In simulation settings three, four, and five, we explore the implication of covariate adjustment when we have mis-specified models. We observed that when we adjust for non-prognostic baseline covariates in a mis-specified model, we similarly experience reductions in power. On the other hand, even in a mis-specifed model, we experience efficiency gains when we adjust for prognostic baseline covariates and these efficiency gains increase depending on the level of prognosis of the associated baseline covariates adjusted. In these simulations, we assume that we have no missing data. In the next post, we explore the implication of covariate adjustment in the case of a continuous outcome variable and missing data. 

## References

<div id="refs"></div>

## Appendices

### Appendix 1A: Code


```r
knitr::opts_chunk$set(
	echo = TRUE,
	message = FALSE,
	warning = FALSE,
	fig.align = 'center',
	fig.pos = 'H',
	dpi = 350,
	tidy.opts = list(width.cutoff = 80, tidy = TRUE)
)
# function to install missing packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, repos='http://cran.rstudio.com/')
  sapply(pkg, require, character.only = TRUE)
}
#install.packages('package_name', dependencies=TRUE, repos='http://cran.rstudio.com/')
packages =c( "tidyverse","knitr", "kableExtra","skimr", "MatchIt", "RItools","optmatch", "ggplot2", "tufte", "tufterhandout", "plotly", "snowfall", "rstan", "gridExtra", "knitr", "gtsummary", "data.table", "GGally", "MASS", "broom")
ipak(packages)
bytes <- file.size("README.Rmd")
words <- bytes/10
minutes <- words/200
# set the seed for reproducibility
set.seed(1)
#setting one:
model.one<-function(nobs=500, a=0, b=1, nreps = 1000){
  
helper.fun=function(){
e=rnorm(nobs)
trt<-c()
for (i in 1:nobs){
  if (i<=nobs/2){
    trt[i]<-0
  }else{
    trt[i]<-1
  }
}
y=a+b*(trt)+e
mod<-lm(y~as.factor(trt))
return(mod)
}

#create 1000 replicates/ simulated datasets
reps=replicate(nreps, helper.fun(), simplify = FALSE)
intercept=reps %>% map_df(tidy) %>% dplyr::filter(term == "(Intercept)")%>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
intercept
trt2<-reps %>% map_df(tidy) %>% dplyr::filter(term == "as.factor(trt)1")%>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#prepare the table of the regression estimates to return
tbl=data.frame(rbind(intercept, trt2))
row.names(tbl)<-c("Intercept", "Treatment1")
return(tbl)
}
#with no covariate adjustment
one<-model.one() %>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the unadjusted regression model in simulation setting one", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)
one
#setting two:
#set.seed(1)
model.two<-function(nobs = 500, nvars = 2, nreps = 1000, rho = 0, a = 0, b = 1, c = 1){
  helper.fun<-function(){
    
#create the treatment variable where n=250 is trt = 0 and 1 elsewhere
trt<-c()
for (i in 1:nobs){
  if (i<=nobs/2){
    trt[i]<-0
  }else{
    trt[i]<-1
  }
}

#generate the variance-covariance matrix
sigma<-matrix(rho, ncol = nvars, nrow = nvars)

#set the values of the leading diagonal to be 1 for variances
diag(sigma)<-1
mat<-mvrnorm(n = nobs, mu = c(0, 0), Sigma = sigma) # shift the mean of the outcome variable: replace y=z and then y=z+a+b*trt

#assign x1 as the 2nd column and z as col 1
x1 = mat[, 2]; z = mat[, 1]

#define new y conditional on values of trt and x
y = z + a + b*trt + c * x1
mat.new<-data.frame(y, trt, x1)

#run the lm on these data 
lm.mod<-lm(y ~ factor(trt) + x1, data = mat.new)
return(lm.mod)
  }
  
#replicate n=1000 times
sims<-replicate(nreps, helper.fun(), simplify = FALSE )

#average intercept
intercept=sims %>% map_df(tidy) %>% dplyr::filter(term == "(Intercept)") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#average treatment
trt2=sims %>% map_df(tidy) %>% dplyr::filter(term == "factor(trt)1") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#average x1
x=sims %>% map_df(tidy) %>% dplyr::filter(term == "x1") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
rslts=data.frame(rbind(intercept, trt2, x))

#clean up the names for the output
row.names(rslts)<-c("Intercept", "Treatment1", "x")
return(rslts)
}
## Case 1: adjusting for non-prognostic baseline covariates: when correlation between x and y is 0
model.two()%>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting two: adjusting for a non-prognostic baseline covariate", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)
## Case 2: adjusting for a fairly prognostic variable (r=0.3)
model.two(rho = 0.3)%>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting two: adjusting for a fairly prognostic variable", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)
## Case 3: adjusting for a prognostic variables (r=0.8)
model.twoc = model.two(rho = 0.8)
model.twoc%>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting two: adjusting for a prognostic variable", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)

#setting 3: adjusting for x1, x2, and x1/x2
#set.seed(1)

model.three<-function(nobs = 500, nvars = 3, nreps = 1000, rho = 0, a = 0, b = 1, theta1 = 1, theta2 = 1, theta3 = 1){
  helper.fun<-function(){
#create the treatment variable where n=250 is treatment = 0 and 1 elsewhere
trt<-c()
for (i in 1:nobs){
  if (i<=nobs/2){
    trt[i]<-0
  }else{
    trt[i]<-1
  }
}

#generate the variance-covariance matrix
sigma<-matrix(rho, ncol = nvars, nrow = nvars)

#set the values of the leading diagonal to be 1 for variances
diag(sigma)<-1
sigma[2,3] = 0; sigma[3,2] = 0
mat<-mvrnorm(n = nobs, mu = c(0, 0, 0), Sigma = sigma, tol = 0.1) #noted idea;replace y=z and then y=z+a+b*trt

#assign x1 as the 2nd column and z as col 1
x1 = mat[, 2]; z = mat[, 1]; x2 = mat[, 3]; x = (x1/x2)

#define new y conditional on values of treatment and X's
y = z + a + b*trt + theta1*x1 + theta2*x2 + theta3 * x
mat.new<-data.frame(y, trt, x, x1, x2)
#run the lm on these data 
lm.mod<-lm(y ~ factor(trt) + x1 + x2 + x, data = mat.new)
#print(lm.mod)
return(lm.mod)
  }
  
#replicate n=1000 times
sims<-replicate(nreps, helper.fun(), simplify = FALSE )

#average intercept
intercept=sims %>% map_df(tidy) %>% dplyr::filter(term == "(Intercept)") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#average treatment
trt2=sims %>% map_df(tidy) %>% dplyr::filter(term == "factor(trt)1") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
#average x1
x1=sims %>% map_df(tidy) %>% dplyr::filter(term == "x1") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
#average x2
x2=sims %>% map_df(tidy) %>% dplyr::filter(term == "x2") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
#average x=x1/x2
x=sims %>% map_df(tidy) %>% dplyr::filter(term == "x") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
rslts=data.frame(rbind(intercept, trt2, x1, x2, x))
row.names(rslts)<-c("Intercept", "Treatment2", "x1", "x2", "x=(x1/x2)")
return(rslts)
}
## Case 1: adjusting for non-prognostic baseline covariates: when correlation between x and y is 0
model.three()%>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting three: adjusting for a non-prognostic baseline covariate", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE) %>% 
  add_footnote("These results are based on an incorrect model specification")
## Case 2: adjusting for a fairly prognostic variable (r=0.3) in a mis-specified model
model.threeb = model.three(rho = 0.3)
model.threeb %>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting three: adjusting for a fairly prognostic variable", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)%>% 
  add_footnote("These results are based on an incorrect model specification")
## Case 3: adjusting for a prognostic variables (r=0.8)
model.threec = model.two(rho = 0.8)
model.threec%>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting three: adjusting for a prognostic variable", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)%>% 
  add_footnote("These results are based on an incorrect model specification")
#Scenario 4: mis-specified model adjusting for x1/x2, (x1-x1bar)^2, (x2-x2bar)^2
model.four<-function(nobs = 500, nvars = 3, nreps = 1000, rho = 0, a = 0, b = 1, theta1 = 1, theta2 = 1, theta3 = 1){
  helper.fun<-function(){
#create the treatment variable where n=250 is trt = 0 and 1 elsewhere
trt<-c()
for (i in 1:nobs){
  if (i<=nobs/2){
    trt[i]<-0
  }else{
    trt[i]<-1
  }
}

#generate the variance-covariance matrix
sigma<-matrix(rho, ncol = nvars, nrow = nvars)

#set the values of the leading diagonal to be 1 for variances
diag(sigma)<-1
sigma[2,3] = 0; sigma[3,2] = 0
mat<-mvrnorm(n = nobs, mu = c(0, 0, 0), Sigma = sigma, tol = 0.1) #noted idea;replace y=z and then y=z+a+b*trt

#assign x1 as the 2nd column and z as col 1
x1 = mat[, 2]; z = mat[, 1]; x2 = mat[, 3]; x = (x1/x2)
x1dev=(x1-mean(x1))^2; x2dev=(x2-mean(x2))^2

#define new y conditional on values of trt and x
y = z + a + b*trt + theta1*x + theta2*x1dev +theta3*x2dev
mat.new<-data.frame(y, trt, x, x1dev, x2dev)
#run the lm on these data 
lm.mod<-lm(y ~ factor(trt) + x + x1dev + x2dev, data = mat.new)
#print(lm.mod)
return(lm.mod)
  }
  
#replicate n=1000 times
sims<-replicate(nreps, helper.fun(), simplify = FALSE )

#average intercept
intercept=sims %>% map_df(tidy) %>% dplyr::filter(term == "(Intercept)") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#average treatment
trt2=sims %>% map_df(tidy) %>% dplyr::filter(term == "factor(trt)1") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#average x=x1/x2
x=sims %>% map_df(tidy) %>% dplyr::filter(term == "x") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
#average x1dev
x1dev=sims %>% map_df(tidy) %>% dplyr::filter(term == "x1dev") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
#average x2dev
x2dev=sims %>% map_df(tidy) %>% dplyr::filter(term == "x2dev") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
rslts=data.frame(rbind(intercept, trt2, x, x1dev, x2dev))
row.names(rslts)<-c("Intercept", "trt2", "x1/x2", "(x1-x1bar)^2", "(x2-x2bar)^2")
return(rslts)
}
## Case 1: adjusting for non-prognostic baseline covariates: when correlation between x and y is 0
model.four()%>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting four: adjusting for non-prognostic baseline covariates", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE) %>% 
  add_footnote("These results are based on an incorrect model specification. The model is under-specified since the main effects of the prognostic baseline variables are not included.")
## Case 2: adjusting for a fairly prognostic variable (r=0.3) in a mis-specified model
model.fourb = model.four(rho = 0.3)
model.fourb %>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting four: adjusting for a fairly prognostic variable", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)%>% 
  add_footnote("These results are based on an incorrect model specification. The model is under-specified since the main effects of the prognostic baseline variables are not included.")
## Case 3: adjusting for a prognostic variables (r=0.8)
model.fourc = model.four(rho = 0.8)
model.fourc%>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting four: adjusting for a prognostic variable", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)%>% 
  add_footnote("These results are based on an incorrect model specification. The model is under-specified since the main effects of the prognostic baseline variables are not included.")
# Scenario 5: Adjust for x1/x2, (x1-x1bar)^2, (x2-x2bar)^2, x1, x2
model.five<-function(nobs = 500, nvars = 3, nreps = 1000, rho = 0, a = 0, b = 1, theta1 = 1, theta2 = 1, theta3 = 1, theta4 =1, theta5 = 1){
helper.fun<-function(){
#create the treatment variable where n=250 is trt = 0 and 1 elsewhere
trt<-c()
for (i in 1:nobs){
  if (i<=nobs/2){
    trt[i]<-0
  }else{
    trt[i]<-1
  }
}

#generate the variance-covariance matrix
sigma<-matrix(rho, ncol = nvars, nrow = nvars)

#set the values of the leading diagonal to be 1 for variances
diag(sigma)<-1
sigma[2,3] = 0; sigma[3,2] = 0
mat<-mvrnorm(n = nobs, mu = c(0, 0, 0), Sigma = sigma, tol = 0.1) 
#assign x1 as the 2nd column and z as col 1
x1 = mat[, 2]; z = mat[, 1]; x2 = mat[, 3]; x = (x1/x2)
x1dev=(x1-mean(x1))^2; x2dev=(x2-mean(x2))^2

#define new y conditional on values of trt and x
y = z + a + b*trt + theta1*x1 + theta2*x2 + theta3*x + theta4*x1dev +theta5*x2dev
mat.new<-data.frame(y, trt, x1, x2, x, x1dev, x2dev)
#run the lm on these data 
lm.mod<-lm(y ~ factor(trt) + x1 + x2 + x + x1dev + x2dev, data = mat.new)
#print(lm.mod)
return(lm.mod)
  }
  
#replicate n=1000 times
sims<-replicate(nreps, helper.fun(), simplify = FALSE )

#average intercept
intercept=sims %>% map_df(tidy) %>% dplyr::filter(term == "(Intercept)") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#average treatment
trt2=sims %>% map_df(tidy) %>% dplyr::filter(term == "factor(trt)1") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#average x1
x1=sims %>% map_df(tidy) %>% dplyr::filter(term == "x1") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#average x2
x2=sims %>% map_df(tidy) %>% dplyr::filter(term == "x2") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))

#average x=x1/x2
x=sims %>% map_df(tidy) %>% dplyr::filter(term == "x") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
#average x1dev
x1dev=sims %>% map_df(tidy) %>% dplyr::filter(term == "x1dev") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
#average x2dev
x2dev=sims %>% map_df(tidy) %>% dplyr::filter(term == "x2dev") %>% dplyr::summarise(estimate = mean(estimate), stderr=mean(std.error))
rslts=data.frame(rbind(intercept, trt2, x1, x2, x, x1dev, x2dev))
row.names(rslts)<-c("Intercept", "trt2", "x1", "x2", "x1/x2", "(x1-x1bar)^2", "(x2-x2bar)^2")
return(rslts)
}
## Case 1: adjusting for non-prognostic baseline covariates: when correlation between x and y is 0
model.five()%>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting five: adjusting for non-prognostic baseline covariates", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE) %>% 
  add_footnote("These results are based on an incorrect model specification. The model is specified to include the main effects of the prognostic baseline variables.")
## Case 2: adjusting for a fairly prognostic variable (r=0.3) in a mis-specified model
model.fiveb = model.five(rho = 0.3)
model.fiveb %>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting five: adjusting for a fairly prognostic variable", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)%>% 
  add_footnote("These results are based on an incorrect model specification. The model is specified to include the main effects of the prognostic baseline variables.")
## Case 3: adjusting for a prognostic variables (r=0.8)
model.fivec = model.five(rho = 0.8)
model.fivec %>% kable(format = "html", 
                      caption = "Regression estimates and standard errors from the adjusted regression model in simulation setting five: adjusting for a prognostic variable", 
                      digits = 4, 
                      col.names = c("Estimate", "Std. Error")) %>% 
  kable_styling(latex_options=c("stripped", "HOLD_position"), full_width=FALSE)%>% 
  add_footnote("These results are based on an incorrect model specification. The model is specified to include the main effects of the prognostic baseline variables.")
```








