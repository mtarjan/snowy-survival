##snovery plover survival analysis
##M tarjan
##September 8, 2017

##install required package
#install.packages("devtools")
#library(devtools) ##required to install packages from github

#devtools::install_github('kfdnm','NMML')
library(kfdnm)

##SIMULATE DATA
# num_kf Target numer of known-fate individuals in the group
# num_years Number of years for the study
# num_surveys Number of surveys within each year
# recruit_rate Average recruitment rate
# init_rate Initial average group size in year 1 survey 1.
# survival_dnm Survival rate for non-kf individuals.
# survival_kf Survival rate for KF individuals.
# perfect_survey_rate The probability that all non-kf individuals are observed in a survey.
# detection The detection probability for abundance surveys.
sim.dat<-sim_data(num_kf = 5, num_years = 10, num_surveys = 8, recruit_rate = 5, init_rate = 15, survival_dnm = .90, survival_kf = .90, perfect_survey_rate = .3, detection = .7)

kf_data<-sim.dat$kf_data
dnm_data<-sim.dat$dnm_data

##LOAD SNPL DATA
##rows must be first ordered by group, then time within group
##create kf_data (known fate) with colnames id, year, survey, CH (1=alive, 0=dead)
##create dnm_data (n-mix) with colnames "year" "survey" "perfect_survey" "Y" "R" "S"  "G" "N" "n"  "group"
  ##perfect_survey = were all non-tagged indivs observed during the survey (1,0)
  ##Y --set to 0 in example
  ##R = the number of new known-fate indivs sampled from population at each time (set to 0 in example)
  ##S -- rbinom(1, N[r-1]-out$R[r-1], survival_dnm) in example ##the number of surviving animals from last time?? yes!
  ##G -- rpois(1,recruit_rate) in example ##This is the number of recruits??
  ##N -- S+G in example ##true number of animals in the pop??
  ##n = observed abundance count at each survey time (N - R)
  ##group = a factor variable indicating the survey group

  ##only actually requires group, n, R for analysis. other variables were used to generate these values in the example

###load SNPL data from database



###
### Create list of fixed parameters
###

fixed_list = list(
  recruit=ifelse(sim.dat$dnm_data$year!=1 & sim.dat$dnm_data$survey==1, NA, 0),
  detection=ifelse(sim.dat$dnm_data$perfect==1, 1, NA)
)

###
### Make likelihood function
### 
lik.func <- make_ikfdnm_likelihood(survival=~1, recruit=~1, detection=~1, kf_survival_effects=NULL, fixed_list=fixed_list, kf_data=kf_data, dnm_data=dnm_data, N_max=50)

###
### Optimize and obtain estimates and variance-covariance matrix
### 

par_start=c(qlogis(0.95), log(4), qlogis(0.5))
#par_start=rep(0,3)
lik.func(par_start)


mle=optim(par_start, lik.func, method="BFGS", control=list(REPORT=1, trace=1), hessian=TRUE)
par=mle$par
se = sqrt(diag(2*solve(mle$hessian)))


###
### Estimates and 95% CI
###

# omega
cat("Survival: \n")
cat(plogis(par[1]), "(", plogis(par[1]-2*se[1]), ",",plogis(par[1]+2*se[1]), ")\n")

# gamma
cat("Recruitment rate:")
cat(exp(par[2]), "(", exp(par[2]-2*se[2]), ",",exp(par[2]+2*se[2]), ")\n")

# p
cat("Detection prob.:\n")
cat(plogis(par[3]), "(", plogis(par[3]-2*se[3]), ",",plogis(par[3]+2*se[3]), ")\n")
