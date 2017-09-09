##Code authored by https://github.com/NMML/kfdnm/blob/master/R/sim_group.R
##Apologies for replicating the code here; I tried to load the package and received a message that it is not available for my version of R, thus the copy-paste

# Simulate group data for analysis with known-fate dynamic N-mixture models
#
# This function generates simulated data for a single group in the correct form for use in the kfdnm package. In addition to the required
# elements for estimation, additional latent information is also provided. The true abundance, recruits, and survivors are retained for user reference. 
# 
# num_kf Target numer of known-fate individuals in the group
# num_years Number of years for the study
# num_surveys Number of surveys within each year
# recruit_rate Average recruitment rate
# init_rate Initial average group size in year 1 survey 1.
# survival_dnm Survival rate for non-kf individuals.
# survival_kf Survival rate for KF individuals.
# perfect_survey_rate The probability that all non-kf individuals are observed in a survey.
# detection The detection probability for abundance surveys.
#
#
#
sim_data = function(num_kf, num_years, num_surveys, recruit_rate, init_rate, 
                    survival_dnm, survival_kf, perfect_survey_rate, detection){
  out = data.frame(year=rep(1:num_years, each=num_surveys))
  out$survey=rep(1:num_surveys, num_years)
  out$perfect_survey=rbinom(num_years*num_surveys, 1, perfect_survey_rate)
  out$Y = 0
  out$R = 0
  N = rep(NA,num_years*num_surveys) 
  N[1] = rpois(1,init_rate)+num_kf
  S = rep(NA,num_years*num_surveys)
  G = rep(0,num_years*num_surveys)
  for(r in 2:(num_years*num_surveys)){
    S[r]=rbinom(1, N[r-1]-out$R[r-1], survival_dnm)
    if(out$survey[r]==1) {
      G[r]=rpois(1,recruit_rate)
    }
    N[r] = S[r]+ G[r]
    #     if(out$survey[r]==1){
    #       out$R[r]= min(num_kf-out$Y[r], N[r])
    #     }
  }
  det=rep(detection, num_years*num_surveys)
  det=ifelse(out$perfect_survey==1, 1, det)
  out$S=S
  out$G=G
  out$N=N
  out$n = rbinom(num_years*num_surveys, N-out$R, det)
  
  out_kf=NULL
  for(i in 1:num_years){
    Y = matrix(nrow=num_surveys, ncol=num_kf)
    Y[1,] = 1
    for(j in 2:num_surveys){
      Y[j,] = rbinom(n=num_kf, size=Y[j-1,], prob=survival_kf)
    }
    CH=as.vector(Y)
    id = rep(paste(i, 1:num_kf, sep="-"), each=num_surveys)
    survey=rep(1:num_surveys, num_kf)
    out_kf = rbind(out_kf, data.frame(id, year=i, survey, CH))
  }
  
  return(list(dnm_data=out, kf_data=out_kf))
}
