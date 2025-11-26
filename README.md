# On the use of the concentration function to compare predictive distributions in ARMA models
## Author: Cristian Cruz-Torres [aut, cre], Fabrizio Ruggeri [aut]

The study explores the use of the **concentration function** a tool based on predictive density ratios to compare the predictive performance of competing time series models. The methodology is illustrated using Bayesian ARMA models and comparison score.

The goal of this project is to demonstrate how the concentration function can:

- Quantify differences between predictive distributions
- Provide a clear, interpretable, probabilistic comparison between two competing models
- Complement traditional metrics with the concentration function, Pietra and Gini index

---

## Structure

- `R/functions.R`  
  - Core density function  
  - General analytic estimator  
  - Model-specific wrappers  
  - RNG for each special case  
  - Numeric MAP estimator  
  - MLE for the Gamma model
