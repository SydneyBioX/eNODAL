# NEMoE 

eNODAL is an R package to analysis relationships between complex experimental conditions and high dimensional omics features. It support two types of experimental conditions as input, which can be both (mutivariate) numerical or categorical variable. eNODAL first use an ANOVA like test procedure to categorize high diemensional omics features into interpretable groups, which are significantly affected by 1) condition1, 2) condition2, 3) additive effect 4) interaction effect. Then eNODAL do consensus clustering for each interpretable group. eNODAL also provides variate of visualization of the result including nutrition geometry framework(NGF).

## Installation

To use this package you can :
 
install it through
 
``` r
remotes::install_github("SydneyBioX/eNODAL")
```
## Vignette

You can find the vignette at our website: https://sydneybiox.github.io/eNODAL/.

## Contact us
If you have any enquiries, especially about running eNODAL model on your own data, then please contact bioinformatics@maths.usyd.edu.au. You can also open an issue on GitHub.
