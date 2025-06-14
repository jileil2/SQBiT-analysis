# Smoothed Quantile Regression for Spatial Data

In this work, we develop a procedure for fitting the partially linear quantile spatial model.


## Quantile Spatial Model (QSM) (Why it is useful?)
Spatial data often exhibit both spatial heterogeneity and tail heterogeneity. In this work, we develop an estimation procedure for a model we call the quantile spatial model. In this model, the conditional quantile function is expressed as a linear combination of covariates with both constant and spatially varying effects. The varying effects are allowed to change smoothly across space. This framework offers considerable flexibility and interpretability, as it captures how the effects of covariates vary across both spatial locations and quantile levels.

## Methodology (How flexible it is?)
We propose a Smoothed Quantile Bivariate Triangulation method (SQBiT), which combines penalized splines on triangulated spatial domains with convolution smoothing of the quantile loss function. This approach offers two key advantages. First, it is well-suited for irregularly shaped spatial domains, where traditional tensor product spline methods can suffer from significant bias. Second, it is computationally efficient. Quantile regression is often computationally demanding in high dimensions, a particular challenge in spatial settings, where a large number of basis functions may be needed for adequate accuracy. By using a convolution-smoothed quantile check loss, which is twice differentiable and locally strictly convex, our method enables faster convergence and more efficient optimization compared to the standard (non-smooth) quantile loss.

## Inference

We discuss the inference of the QSM for constant and varying coefficients separately. 

1. For varying coefficients, we develop an interval estimation procedure based on wild bootstrap. Its implementation is based on `smqsvcm_admm.wb.R`, see Section Main Functions.
2. For constant coefficients, we develop interval estimation procedures based on asymptotic normality and wild bootstrap. The former is based on `SQBiT.R` or `tune.lambda.R` while the latter is based on `smqsvcm_admm.wb.R`. See Section Main Functions for details.

## Instructions for use

### Dependencies

Note that our codes heavily depend on the packages Rcpp, RcppArmadillo, and BPST. In particular, the linear algebra is built on Rcpp. The packages Rcpp and RcppArmadillo can be downloaded using `install.packages()` in R. The package BPST can be downloaded from Github using `devtools::install_github(FIRST-Data-Lab/BPST)`.

### Main Functions

There are several main functions for implementing SQBiT. 

1. `SQBiT.R` contains the codes for SQBiT for a given penalization parameter.
2. `tune.lambda.R` contains the codes for SQBiT with penalization parameters selected by GACV.
3. `smqsvcm_admm.wb.R` contains the codes for wild bootstrap inference of SQBiT
4. `cv.pred.SQBiT.R` contains the codes for computing k-fold cross-validation prediction error of SQBiT for given bandwidth h and a penalization parameter.

### Simulation Examples

In the Simulation Examples folder, it contains the necessary codes and data for simulation examples. Note that

1. `plm_sim_t.R` contains the codes for the simulation study with homoskedastic t(2) error.

2. `plm_sim_t_heter.R` contains the codes for the simulation study with heteroskedastic t(2) error.

3. `plm_sim_mn_heter.R` contains the codes for the simulation study with heteroskedastic mixture normal error.


### Application Examples

In the Application Examples folder, it contains the necessary codes and data for US mortality data analysis. Note that

1. `main_tau50.R` contains the codes for estimating the quantile spatial model for the median case.

2. `main_tau75.R` contains the codes for estimating the quantile spatial model for the third quartile case.

## Contact
If you have any questions about the codes, please contact jileil2@gwmail.gwu.edu.
