covsim
==========

This package contains  implementations of three simulation procedures
that simulate from a multivariate distributions with a pre-specified
covariance matrix. In addition, the user may pre-specify marginal information. The main function is vita, which implements the very general VIne-To-Anything (VITA) algorithm of Grønneberg and Foldnes (Psychometrika, 2017).


How to install
--------------


You can install:

-   the stable release on CRAN:

    ``` r
    install.packages("covsim")
    ```

-   the latest development version:

    ``` r
    devtools::install_github("njaalf/covsim")
    ```

------------------------------------------------------------------------

Package overview
----------------



### VITA

The user must completely specify the marginals from the distributions available in the stats package. Also the user must specify the population covariance matrix. Note that the variances in this matrix must equal the population variances of the marginal distributions. 
A regular vine may also be provided, that is, a list of bicopulas and a sequence for their calibrations, i.e. a sequence of trees. If no vine
is provided, a default D-vine is assumed. Note that in many cases the vine, the marginals and the covariance matrix are not compatible. That is, no vine distribution exists whose covariance equals the target covariance. Also note that calibration of the vine is computationally burdensome. To speed up calibration, one could reduce the numpoints and/or the Nmax parameters. This would reduce the precision of the adherence of the vita vine to the target population covariance matrix.

### IG

In addition to vita, a simpler and much more restricted simulation algorithm is provided in the rIG function, which implements the independent generator approach of Foldnes and Olsson (MBR, 2016)
The IG approach is comparable to the Vale-Maurelli approach, in that
it is fast and precise. Only skewness and kurtosis of the marginals is provided. The marginal distributions are therefore not fully controlled using IG. 

### PLSIM 

has a similar interface as that of IG. Skewness and kurtosis is pre-specified. PLSIM is flexible beyond the default values. Numsegments can be increased, or the breakpoints manually set by providing gammalist.  


References
----------
Grønneberg, S., Foldnes, N. Covariance Model Simulation Using Regular Vines. Psychometrika 82, 1035–1051 (2017). https://doi.org/10.1007/s11336-017-9569-6

Njål Foldnes & Ulf Henning Olsson (2016) A Simple Simulation Technique for Nonnormal Data with Prespecified Skewness, Kurtosis, and Covariance Matrix, Multivariate Behavioral Research, 51:2-3, 207-219, DOI: 10.1080/00273171.2015.1133274

Njål Foldnes & Steffen Grønneberg (2021) Non-normal data simulation using piecewise linear transforms. In review.


