# On the use of the concentration function to compare predictive distributions in ARMA models
## Author: Cristian Cruz-Torres [aut, cre], Fabrizio Ruggeri [aut]

The study explores the use of the **concentration function** a tool based on predictive density ratios to compare the predictive performance of competing time series models. The methodology is illustrated using Bayesian ARMA models and comparison score.

The goal of this project is to demonstrate how the concentration function can:

- Quantify differences between predictive distributions
- Provide a clear, interpretable, probabilistic comparison between two competing models
- Complement traditional metrics with the concentration function, Pietra and Gini index

**Methodology**
- Bayesian ARMA/SARIMA via bayesforecast::stan_sarima (HMC in Stan, with posterior draws).
- Predictive densities computed by Monte Carlo averaging over posterior draws.
- Concentration functions: w.r.t. uniform reference (Lorenz-style). Relative concentration of one model vs another, ordering outcomes by the predictive likelihood ratio.
- Inequality measures: Pietra index (max vertical gap). Gini concentration coefficient (2*area).
- Predictive scores: NLPD (negative log predictive density). CRPS (Monte Carlo from posterior predictive).
- Marginal log-likelihood: shifted-gamma estimator (Raftery et al.).
