# The roles  of pre-season immunity, age, viral shedding, and community exposures in shaping influenza household transmission dynamics
## Molly Sauter, Jackie Kleynhans, Jocelyn Moyes, Meredith L McMorrow, Florette K Treurnicht, Orienka Hellferscee, Anne von Gottberg, Nicole Wolter, Amelia Buys, Lorens Maake, Neil A Martinson, Kathleen Kahn, Limakatso Lebina, Katlego Mothlaoleng, Floidy Wafawanaka, Francesc Xavier Gómez-Olivé, Stefano Tempia, Bryan Grenfell, Cecile Viboud, Kaiyuan Sun, Cheryl Cohen
These files provide a mock data set (the data is not yet public) and the code used to reproduce the viral shedding model, household transmission model, and figures from this manuscript. This mock data set only runs the model and creates figures for one subtype/lineage (A(H3N2)); all models and figures in the paper were run seperately for each of the four influenza virus subtypes/lineages observed in the real data set. 

# System Requirements

### Operating System: 
Windows, Linux (Ubuntu 20.04 or later), macOS 10.15 or later

### Programming Language:
R (version 4.1.0 or later)
Required R Packages:
tidyverse, dpylr, ggplot2, rstan, purrr, shinystan, bayesplot, lme4, sjPlot, sjmisc, sjlabelled, lmtest, nlme, broom.mixed, lmerTest, lares, gridExtra, grid, jtools, mertools, zoo 

### Tested Environments: 
Windows Intel Core i5 processor, 8GB RAM 
R version: 4.1.0

# Installation
### Install R
Ensure that R (version 4.1.0 or later) is installed on your system. You can download it from the CRAN website (https://cran.r-project.org/) for the appropriate operating system. 

### Install RStan to allow R to interface with Stan
First, R needs to be configured to be able to compile C++ code. Instructions to do so depending on the operating system can be found here: 

Windows: https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Windows

Mac: https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Mac

Linux: https://github.com/stan-dev/rstan/wiki/Configuring-C-Toolchain-for-Linux

Then users can install the latest version of rstan using: 
```
install.packages("rstan", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
```

### Clone the repository and install the required R packages
To clone the entire repository, users can run the following in a command-line interface: 
```
git clone https://github.com/msauter19/phirst-flu-shed-trans.git  
cd phirst-flu-shed-trans  
```
Users can install all required packages in the R terminal using: 
```
install.packages(c("tidyverse", "dpylr", "ggplot2", "rstan", "purrr", "shinystan", "bayesplot", "lme4", "sjPlot", "sjmisc", "sjlabelled", "lmtest", "nlme", "broom.mixed", "lmerTest", "lares", "gridExtra", "grid", "jtools", "mertools", "zoo")) 
```
And ensure that the working directory is set to the cloned repository folder using:
```
setwd("path_to/phirst-flu-shed-trans")
```
Installation time: approximately 10-15 minutes to install all necessary components on a standard desktop computer. 

# Demo
### Viral Shedding Model: 
After ensuring the successful load and access of the [![mock_dataset.csv](https://github.com/msauter19/phirst-flu-shed-trans/blob/main/mock_dataset.csv)] and [![flu_ct_model.stan](https://github.com/msauter19/phirst-flu-shed-trans/blob/main/flu_ct_model.stan)] files, a demo of the viral shedding model can be run with [![shed_model_run.R](https://github.com/msauter19/phirst-flu-shed-trans/blob/main/shed_model_run.R)]. The regressions to demo the results of the viral shedding model are found in [![shed_regression.R](https://github.com/msauter19/phirst-flu-shed-trans/blob/main/shed_regression.R)].

This model should take 1-5 minutes to run. 

### Household Transmission Model: 
After running the viral shedding model, the results are read into the household transmission model within [![transmission_model.R](https://github.com/msauter19/phirst-flu-shed-trans/blob/main/transmission_model.R)]. 

This model should take 0-1 minutes to run.

### Visualizations 
Demos for all figures in the main text of the manuscript are produced in [![figures.R](https://github.com/msauter19/phirst-flu-shed-trans/blob/main/figures.R)]. Demos for the figures in the supporting information are produced in [![supporting_figs.R](https://github.com/msauter19/phirst-flu-shed-trans/blob/main/supporting_figs.R)]. 


