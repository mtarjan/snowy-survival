##snovery plover survival analysis
##M tarjan
##September 8, 2017

##install required package
#install.packages("devtools")
#library(devtools) ##required to install packages from github

#devtools::install_github('kfdnm','NMML')
library(kfdnm)

# num_kf Target numer of known-fate individuals in the group
# num_years Number of years for the study
# num_surveys Number of surveys within each year
# recruit_rate Average recruitment rate
# init_rate Initial average group size in year 1 survey 1.
# survival_dnm Survival rate for non-kf individuals.
# survival_kf Survival rate for KF individuals.
# perfect_survey_rate The probability that all non-kf individuals are observed in a survey.
# detection The detection probability for abundance surveys.
sim.dat<-sim_data(num_kf = 3, num_years = 10, num_surveys = 3, recruit_rate = .5, init_rate = 30, survival_dnm = .25, survival_kf = .25, perfect_survey_rate = .2, detection = .6)

##make likelihood function
kfdnm_lik<-make_ikfdnm_likelihood(kf_data = sim.dat$kf_data, dnm_data = sim.dat$dnm_data)

##
kfdnm.samp<-make_ikfdnm_sampler()

