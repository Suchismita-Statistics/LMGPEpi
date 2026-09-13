# LMGPEpi
This package contains all the functions required for the paper, Latent Mechanistic Gaussian Processes for Partially Observed Stochastic Epidemics.

To install:

```
install.packages("devtools")
devtools::install_github("Suchismita-Statistics/LMGPEpi")
```

# About the package


The package includes two models: one for accurate count data and one for under-reported count data. Inference uses LMGP (Latent Mechanistic Gaussian Processes), which proposes a computationally efficient method for estimating the parameters: infection rate ($\beta$), recovery rate ($\gamma$) and the limiting proportion of susceptible and infected ($\rho$). For the under-reported model, it also estimates the under-reporting probability.

This package contains all the R functions required to implement the method. 
LMGP() is the main function to apply the method to new daily count data, which simulates HMC posterior samples via stan. We used rstan::version 2.32.7.  
Please see the Example.R file for usage of the key functions. To do so, write the following code in the terminal:  

```
git clone https://github.com/Suchismita-Statistics/LMGPEpi.git
cd LMGPEpi
Rscript Example.R
```

DRC_Ebola_Application.R applies the method to the Ebola 2026 epidemic in the DRC. The stan files are in the "inst/stan" folder. 
