# LMGPEpi
This package contains all the functions required for the paper, Latent Mechanistic Gaussian Processes for Partially Observed Stochastic Epidemics.

To install:

```
install.packages("devtools")
devtools::install_github("Suchismita-Statistics/LMGPEpi")
```

# About the package

LMGP (Latent Mechanistic Gaussian Processes) proposes a computationally efficient method about the parameters -  infection rate ($\beta$) and recovery rate ($\gamma$). It further estimates the limiting proportion of susceptible and infected ($\rho$). The DSA method does not necessarily require information on the number of initial susceptible individuals. It can also be estimated.  

This package contains all the R codes required to implement the paper. 
new_data_ct_lkd() is the main function to apply count likelihood for a new daily count data, which simulates HMC posterior samples via stan. We used rstan::version 2.21.8.  
Please look at the Example.R file to see the usage of the important functions. To do so, write the following code in the terminal:  

```
git clone https://github.com/Suchismita-Statistics/DSA.CountData.git
cd DSA.CountData
Rscript Example.R
```
