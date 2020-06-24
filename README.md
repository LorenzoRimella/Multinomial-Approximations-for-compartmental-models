# Multinomial-Approximations-for-compartmental-models
Explanation on how to use the Multinomial Approximations for compartmental models proposed in "Inference in Stochastic Epidemic Models via Multinomial Approximations" by N. Whiteley and L. Rimella. The repository contains two tutorials:
- application to the Ebola outbreak in the Democratic Republic of Congo in 1995 (in Python);
- application to the beginning of the COVID-19 epidemic in Wuhan (in R).

## Requirements
EBOLA:
- The codes are written in Python and they can be used to run our experiments.
  Once you open the folder "EBOLA" you can find the jupyter notebook "Tutorial.ipynb"
  which explains how to use the main functions in our codes.
- You require the following libraries:
  * numpy
  * pandas
  * scipy
  * matplotlib

COVID-19:
- The codes are written in R and they follow the implementation available at
  https://github.com/ adamkucharski/2020-ncov/. Once you open the folder COVID-19
  there are two other folders: "Kucharski-modified", which shows how to run Kucharski's
  code; "MultinomialApproximation", which explains how to use our code. Both folder 
  contains a script "Tutorial.R" which simply reproduce the plots about our experiment
  on COVID-19. 
- You require the packages:
  * foreach
  * lubridate
  * magrittr
  * coda
  * tidyverse
  * rootSolve
  * mgcv
  * doMC: this might give some problems to Windows' users, you could simply avoid using it
          and not run the code in parallel
