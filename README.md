# brgeegor
The SAS codes which implement simulation study and fit the mothodology described in the paper "A bias-reduced generalized estimating equation approach for proportional odds models with small sample longitudinal ordinal data", by Yukio Tada and Tosiya Sato (20xx, BMC Medical Research Methodology). 

## brgeegor.sas
The SAS macro to apply the bias-reduced generalized estmating equation (BR-GEE) to estimate the regression parameter of the proportional odds model (POM) for longitudinal ordinal data.  

## gen_ordinal_gor3.sas
The SAS macro to generate longitudinal ordinal data with the within-subject association specified by global odds ratio. I referred to a SAS code described at an appendix in "Generating correlated discrete ordinal data using R and SAS IML", by Noor Akma Ibrahim, Suliadi Suliadi (2011) to create this program. 

## sim_genord_gor3_d4_h1.sas
The SAS macro to implement the simulation study to generate longitudinal ordinal data under alternative hypothesis and to apply BR-GEE to estimate the regression parameter of the POM using two macros.  

## sim_genord_gor3_d4_h0.sas
The SAS macro to implement the simulation study to generate longitudinal ordinal data under null hypothesis and to apply BR-GEE to estimate the regression parameter of the POM using two macros.  
