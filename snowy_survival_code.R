##snovery plover survival analysis
##M tarjan
##September 8, 2017
##

##install required package
#install.packages("devtools")
#library(devtools) ##required to install packages from github

##Please download and install Rtools 3.4 from http://cran.r-project.org/bin/windows/Rtools/
#devtools::install_github('kfdnm','NMML')
library(kfdnm)

library(dplyr) ##required for groupby

library(tidyr)

library(ggplot2)

library(stringr) ##required for string processing

library(RODBC) ##required to pull data from access databases

library(lubridate)

##SIMULATE DATA
# num_kf Target number of known-fate individuals in the group
# num_years Number of years for the study
# num_surveys Number of surveys within each year
# recruit_rate Average recruitment rate = additions of individuals to group
# init_rate Initial average group size in year 1 survey 1.
# survival_dnm Survival rate for non-kf individuals.
# survival_kf Survival rate for KF individuals.
# perfect_survey_rate The probability that all non-kf individuals are observed in a survey.
# detection The detection probability for abundance surveys.
sim.dat<-sim_data(num_kf = 10, num_years = 5, num_surveys = 7, recruit_rate = 5, init_rate = 15, survival_dnm = .70, survival_kf = .70, perfect_survey_rate = .3, detection = .65)

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

##known fate individuals
wb<-"S:/Science/Banding Database/SFBBO Banding Database.accdb" ##file path for use when connected to SFBBO server
#wb<- "C:/Users/max/Desktop/Tarjan/Science/Plover/SFBBO Banding Database copy 22Sep2017.accdb" #filepath for work locally on Birds25
con<-odbcConnectAccess2007(wb) ##open connection to database

##get resight data
qry<-"SELECT rrBandNumber AS id, 'NA' AS Age, b.Sex, rrDate AS [Date], year(rrDate) AS year, rrLocation AS Location, 'resight' AS Type, IIF(rrStatus = 'Alive', 1, 0) AS CH FROM ResightRecords as r
  LEFT OUTER JOIN BandingRecords AS b ON r.rrBandNumber = b.BandNumber
  WHERE rrSpeciesCode = 'SNPL' AND rrStatus IN ('Alive', 'Dead') "

kf.dat1<-sqlQuery(con, qry); head(kf.dat1) ##import the queried table

##get capture data
##assumes all birds alive when banded
qry<-"SELECT BandNumber AS id, Age, Sex, CaptureDate AS [Date], year(CaptureDate) AS year, SNO.snoNestPondID AS Location, 'capture' AS Type, 1 AS CH FROM BandingRecords AS BR
  LEFT OUTER JOIN SNPLnestOffspring AS SNO ON BR.BRautoID = SNO.snoBandingRecordID
  WHERE SpeciesCode = 'SNPL' "

kf.dat2<-sqlQuery(con, qry); head(kf.dat2) ##import the queried table

##get capture location for parents
qry<-"SELECT IIF(snpMaleBandNum IS NULL, snpFemaleBandNum, snpMaleBandNum) AS id, snpNestYear AS year, snpNestPondID AS Location
FROM snplNestParents"

parent.locs<-sqlQuery(con, qry); head(parent.locs)

##when finished with db, close the connection
odbcCloseAll()

##remove records with unknown id
kf.dat1<-kf.dat1[complete.cases(kf.dat1$id),]

##remove NA from parent locs
parent.locs<-parent.locs[complete.cases(parent.locs),]
##add parent locs to kf.dat2
kf.dat3<-dplyr::left_join(x=kf.dat2, y=subset(parent.locs, select=c(id, year, Location)), by = c("id","year"))
kf.dat2$Location[which(is.na(kf.dat2$Location))]<-kf.dat3$Location.y[which(is.na(kf.dat2$Location))] ##for NA locations in original dataset, replace them with parent locations in kf.dat3

##missing locations provided by Ben
mlocs<-read.csv("birds without ponds.csv")
mlocs$Location<-mlocs$Pond.Association
mlocs<-subset(mlocs, select=c(id, Age, Sex, Date, year, Location, Type, CH))
mlocs<-mlocs[complete.cases(mlocs),]
kf.dat4<-dplyr::left_join(x=kf.dat2, y=subset(mlocs, select=c(id, year, Location)), by = c("id","year"))
kf.dat2$Location[which(is.na(kf.dat2$Location))]<-kf.dat4$Location.y[which(is.na(kf.dat2$Location))] ##for NA locations in original dataset, replace them with  locations in mlocs (missing locations)

##combine resight and capture data
kf.dat<-unique(rbind(kf.dat2, kf.dat1))

##fix spelling of crittenden->don
#kf.dat$Location[which(kf.dat$Location=="Crittenden Marsh East")]<-"Crittendon Marsh East"


##add pond group to locations
pond.groups<-read.csv("pond names.csv")
kf.dat<-dplyr::left_join(x=kf.dat, y=subset(pond.groups, select=c(Pond, Pond.Group)), by = c("Location" = "Pond"))


##ADD CHICKS PRESUMED DEAD ###
##change to add assumption that any chick not seen alive is dead by fledging. but when do we assume they died? randomly assign them a time that is within the distribution of deaths for known-dead individuals. or is it a better assumption that they were dead the day after they were last observed?? that's maybe a better assumption. one week after they were last observed
broods<-read.csv("Brood Groups.csv")
broods<-broods[,1:3]

##add brood to kf.dat
kf.dat<-dplyr::left_join(x=kf.dat, y=subset(broods, select=c(Nest.Group, Band.Number)), by = c("id" = "Band.Number"))


for (j in 1:length(unique(kf.dat$id))) { ##for each bird
  dat.temp<-subset(kf.dat, id==unique(kf.dat$id)[j])
  if (length(which(dat.temp$CH==0))>0 ) {next} ##if we already know the death date or the brood group is unknown, then move to the next bird
  if (is.na(dat.temp$Nest.Group[1]) | dat.temp$Nest.Group[1]=="NA" | dat.temp$Nest.Group[1]=="") {next}##if the brood group is unknown, then move to the next bird
  brood.temp<-subset(kf.dat, Nest.Group==dat.temp$Nest.Group[1])
  b.day<-dat.temp$Date[which.min(dat.temp$Date)] ##focal bird bday
  last.seen<-dat.temp$Date[which.max(dat.temp$Date)] ##date focal bird last observed alive
  ##next date when brood-mates were observed after "last seen" date
  brood.dates<-subset(brood.temp, Date > last.seen)$Date
  if (length(brood.dates)==0) {next} ##if the brood is never observed after the last seen date of the focal bird, move to next bird
  brood.resight<-min(brood.dates)
  if(brood.resight > b.day+ 30*24*3600) {next} ##if the sighting is after the fledge date, move to the next bird
  kf.dat<-rbind(kf.dat, data.frame(id=dat.temp$id[1], Age = dat.temp$Age[1], Sex= dat.temp$Sex[1], Date=brood.resight, year=format(brood.resight, format="%Y"), Location="", Type="missing", CH=0, Pond.Group = "", Nest.Group=dat.temp$Nest.Group[1]))
}

#write.csv(subset(kf.dat, Type=="missing"), "presumed_dead.csv", row.names=F)

##NEED TO ASSIGN LOCATION WHERE BIRD FIRST APPEARED?? OTHERWISE AN INDIVIDUAL WILL BE ANALYZED IN TWO DIFFERENT GROUPS
##their group is their first CAPTURE location
kf.dat$group<-rep(NA, nrow(kf.dat))
for (j in 1:length(unique(kf.dat$id))) {
  kf.dat$group[which(kf.dat$id==unique(kf.dat$id)[j])]<-as.character(kf.dat$Pond.Group[which(kf.dat$id == unique(kf.dat$id)[j] & kf.dat$Type=="capture")][1])
}

##create kf_data (known fate) with colnames id, year, survey, CH (1=alive, 0=dead)
#kf.dat<-subset(kf.dat, select=c(id, year, Date, CH, Location))
#kf.dat<-subset(kf.dat, select=c(id, year, Date, CH, Location), subset= str_detect(string = kf.dat$Location, pattern = "^E") ) ##select only sites in eden landing

##remove birds with no group
kf.dat<-subset(kf.dat, group!="")

##order data by group and then time
kf.dat<-kf.dat[order(kf.dat$group, kf.dat$id, kf.dat$Date),]; head(kf.dat)

##Ages
##https://www.pwrc.usgs.gov/bbl/manual/age.cfm
##0=UNK
##1 = AHY after hatching year
##2=HY
##4=LOCAL young bird incapable of substantial flight
##5=SY*
##6=ASY*
##7=TY*
##8=ATY*
##*NOTE NEEDED

##unknown fate individuals
wb<-"S:/Science/Waterbird/Databases - enter data here!/SNPL/SNPL.accdb"  ##file path for use when connected to SFBBO server
#wb<- "C:/Users/max/Desktop/Tarjan/Science/Plover/SNPL copy 22Sep2017.accdb" #filepath for work locally on Birds25
con<-odbcConnectAccess2007(wb) ##open connection to database

qry<-"SELECT year(i.Date) AS year, i.Date, s.SurveyID AS survey, 0 AS perfect_survey, 0 AS R, s.TotalObserved AS n, p.Complex AS [complex], p.PondNumber, IIF(IsNull(NumMales), 0, NumMales) + IIF(IsNull(NumFemales), 0, NumFemales) + IIF(IsNull(NumUnknown), 0, NumUnknown) + IIF(IsNull(NumJuveniles), 0, NumJuveniles) AS adults, AgeSexStatus
  FROM ([SNPL SURVEY] AS s
  LEFT OUTER JOIN [SNPL Observers Info] AS i ON s.SurveyID = i.SurveyID)
  LEFT OUTER JOIN LookupPond AS p ON s.PondNumber = p.PondNumber
  WHERE p.Complex <> 'Error With Record'"

uf.dat<-sqlQuery(con, qry); head(uf.dat) ##import the queried table

##when finished with db, close the connection
odbcCloseAll()


##EXCLUDE CHICK COUNTS; ONLY ADULTS ##
uf.dat<-subset(uf.dat, AgeSexStatus %in% c("F", "J", "M", "None", "U", "X")) ##exclude counts of chicks
uf.dat$adults2<-ifelse(uf.dat$AgeSexStatus =="X", uf.dat$adults, uf.dat$n)
uf.dat$n<-uf.dat$adults2
uf.dat<-subset(uf.dat, select=-c(adults,AgeSexStatus, adults2))


##clean up data
##remove nonsensical year
uf.dat<-subset(uf.dat, year > 1980 & year < 2018)

##restrict surveys to "window surveys" when only reproductive individuals are present (ie exclude counts of migratory birds)
#window.dates<-read.csv("SNPL Breeding Window Dates.csv")
#window.dates$start.date<-as.Date(window.dates$start.date, "%m/%d/%Y")
#window.dates$end.date<-as.Date(window.dates$end.date, "%m/%d/%Y")

#uf.dat.win<-dim(0) ##counts of untagged indivs during the "window surveys" for each year
#for (j in 1:nrow(window.dates)) {
#  year.temp<-window.dates$Year[j]
#  dat.temp<-subset(uf.dat, year==year.temp)
#  win.temp<-window.dates[j,]
#  out.temp<- subset(dat.temp, as.Date(dat.temp$Date) >= win.temp$start.date & as.Date(dat.temp$Date) <= win.temp$end.date)
  
#  if (nrow(out.temp)>0) {
    ##for survey dates that apply to all locations
#    if (win.temp$locations=="all") {
#      uf.dat.win<-rbind(uf.dat.win, out.temp)
#    }
    
    ##for survey windows that apply to select locations
#    if (win.temp$locations=="Hayward") {
#      uf.dat.win<-rbind(uf.dat.win, subset(out.temp, complex=="Hayward"))
#    }
    
#    if (win.temp$locations=="Hayward, Napa") {
#      uf.dat.win<-rbind(uf.dat.win, subset(out.temp, complex %in% c("Hayward", "Napa")))
#    }
    
#    if (win.temp$locations=="all-Hayward") {
#      uf.dat.win<-rbind(uf.dat.win, subset(out.temp, complex !=  "Hayward"))
#    }
    
#    if (win.temp$locations=="all-Hayward, Napa") {
#      uf.dat.win<-rbind(uf.dat.win, subset(out.temp, complex != "Hayward" & complex != "Napa"))
#    }
#  }
#}

#uf.dat<-uf.dat.win

###ALTERNATIVE: RESTRICT TO MAY AND JUNE
uf.dat$month<-as.numeric(format(uf.dat$Date, "%m"))
uf.dat$week<-as.numeric(week(uf.dat$Date))
uf.dat$wday<-as.numeric(wday(uf.dat$Date))
#uf.dat$wday<-wday(uf.dat$Date, label=T)
#uf.dat<-subset(uf.dat, month %in% c(5,6) & PondNumber %in% read.csv("Ponds Surveyed Weekly.csv", header=F)[,1] & year>=2004) ##use data in May and June from ponds that are surveyed weekly; also choose years when groups were tracked closely

dat.sub<-subset(uf.dat, year==3000)
for (j in 2004:2017) {
  ##find the first monday in May
  year.temp<-j
  for (i in 1:7) {
    date.temp<-as.Date(str_c(year.temp, "-05-0", i), format="%Y-%m-%d")
    if (wday(date.temp)==2) {first.mon<-date.temp; break}
  }
  ##find the last saturday in June
  for (i in 30:22) {
    date.temp<-as.Date(str_c(year.temp, "-06-", i), format="%Y-%m-%d")
    if (wday(date.temp)==6) {last.sat<-date.temp; break}
  }
  
  dat.sub.temp<-subset(uf.dat, Date >=first.mon & Date <=last.sat)
  dat.sub<-rbind(dat.sub, dat.sub.temp)
}

uf.dat<-dat.sub

##FIGURE OUT WHAT LOCATIONS WE SHOULD MAKE THE "GROUP" VARIABLE
uf.dat$group<-uf.dat$complex

##use Ben's groups
#uf.dat<-dplyr::left_join(x=uf.dat, y=subset(pond.groups, select=c(Pond, Pond.Group)), by = c("PondNumber" = "Pond"))
#uf.dat$group<-uf.dat$Pond.Group ##pond group as the location
#uf.dat$group<-uf.dat$PondNumber ##ponds as their own location

## GROUP BY POND ##

##get summary n for a given date and location
#uf.dat.sum<-uf.dat %>% group_by(Date, PondNumber) %>% mutate(group.n=sum(n)) %>% data.frame()
##remove n and remove duplicates; rename survey.n to n
#uf.dat.sum<-subset(uf.dat.sum, select = -c(n)); uf.dat.sum<-unique(uf.dat.sum)
#names(uf.dat.sum)[length(uf.dat.sum)]<-"n"
#head(uf.dat.sum)

## ADD SURVEYS BASED ON 2-WEEK INTERVALS ##
#uf.dat$survey<-ifelse(chron::days(uf.dat$Date)< 16, str_c(uf.dat$year,".", uf.dat$month, ".1"), str_c(uf.dat$year,".", uf.dat$month, ".2"))

uf.dat$survey<-str_c(uf.dat$year, ".", uf.dat$week)

##ALTERNATIVE TO GROUP BY WEEK
uf.dat.wk<-dim(0)
for (j in 1:length(unique(uf.dat$group))) {
  dat.temp<-subset(uf.dat, group==unique(uf.dat$group)[j]) ##get data for a given group
  
  ## SUM BY POND FOR EACH DATE
  
  dat.temp<-dat.temp %>% group_by(Date, PondNumber) %>% mutate(pond.date.n=sum(n)) %>% data.frame() #sum n by pond
  dat.temp<-subset(dat.temp, select=-n) ##remove old n before summarized
  dat.temp<-unique(dat.temp) #remove duplicates
  
  ## TAKE MAX COUNT FOR PONDS WITH MULTIPLE DATES WITHIN THE SAME SURVEY/WEEK
  
  dat.temp<-dat.temp %>% group_by(survey, PondNumber) %>% mutate(pond.survey.n=max(pond.date.n, na.rm=T)) %>% data.frame()
  dat.temp<-subset(dat.temp, select=-c(pond.date.n, Date, wday)) ##remove old n before summarized
  dat.temp<-unique(dat.temp) #remove duplicates
  
  
  #pond.count<-length(unique(dat.temp$PondNumber)) ##get unique number of ponds
  #dat.temp<-spread(dat.temp, PondNumber, pond.survey.n)
  
  ## SUM BY GROUP
  
  ##get the sum of the counts for all ponds in the pond group. note that not all ponds were necessarily counted on each date
  ##spread function above allows us to change NA.RM to F if want to know where we are missing counts for some ponds
  #if(pond.count>1) {
  #  dat.temp$n<-apply(dat.temp[,(ncol(dat.temp)-pond.count+1):ncol(dat.temp)], FUN=sum, na.rm=T, MARGIN = 1)
  #} else {
  #  dat.temp$n<-last(dat.temp)
  #}
  
  dat.temp<-dat.temp %>% group_by(survey, group) %>% mutate(n=sum(pond.survey.n, na.rm=T)) %>% data.frame()
  dat.temp<-subset(dat.temp, select=-c(pond.survey.n, PondNumber)) ##remove old n before summarized
  dat.temp<-unique(dat.temp) #remove duplicates
  
  dat.temp<-subset(dat.temp, select=c(year, month, survey, perfect_survey, R, group, n))
  uf.dat.wk<-rbind(uf.dat.wk, dat.temp)
}
##take the max if more than one count for each date
#uf.dat.wk<-uf.dat.wk %>% group_by(survey, Pond.Group) %>% mutate(n.max=max(n)) %>% data.frame()
#uf.dat.wk$n<-uf.dat.wk$n.max; uf.dat.wk<-subset(uf.dat.wk, select= -n.max) ##replace n with n.max and remove duplicates
#uf.dat.wk<-uf.dat.wk[which(duplicated(uf.dat.wk)==F),]

##remove na counts
uf.dat.wk<-subset(uf.dat.wk, is.na(n)==F)

uf.dat.sum<-uf.dat.wk

#uf.dat.sum$group<-uf.dat.sum$Pond.Group

##order by group, then time
uf.dat.sum<-uf.dat.sum[order(uf.dat.sum$group, uf.dat.sum$survey),]; head(uf.dat.sum)

##OVERWRITE WITH SEQUENTIAL SURVEY NUMBER; requires that the data are ordered by group and then by date
##ADD SURVEY NUMBER
##rename survey to be all 1, since we will represent each yearly count during the window survye as the one and only survey
uf.dat.sum$survey<-rep(1, nrow(uf.dat.sum))

for (j in 1:nrow(uf.dat.sum)) {
  if (duplicated(subset(uf.dat.sum, select=c(year, group)))[j]) { ##if this group and year already appeared
    uf.dat.sum$survey[j]<-uf.dat.sum$survey[j-1]+1
  }
}

#dnm_data<-subset(uf.dat.sum, group == "Eden Landing")
dnm_data<-uf.dat.sum
kf_data<-kf.dat

##visualize the data
vis <- ggplot(data = dnm_data, aes(x = year, y = n))
vis <- vis + geom_point()
vis <- vis + facet_grid(facets = dnm_data$group~., scales = "free_y")
vis

##plot of counts across years for particular ponds
j<-0
j<-j+1
fig <- ggplot(data = subset(dnm_data, group == unique(dnm_data$group)[j]), aes(x = survey, y = n))
fig <- fig + geom_point()
fig <- fig + facet_wrap(~year, scales = "free_x")
fig <- fig + xlab(label = str_c(unique(dnm_data$group)[j], " survey number"))
fig

##visual data for known-fate birds
j<-0
j<-j+1
fplot <- ggplot(data = subset(kf_data, subset = id == unique(kf_data$id)[j]), aes(x = Date, y = CH))
fplot <- fplot + geom_point()
fplot

###
### Create list of fixed parameters
###

#fixed_list = list(
#  recruit=ifelse(dnm_data$year!=1 & dnm_data$survey==1, NA, 0),
#  detection=ifelse(dnm_data$perfect==1, 1, NA)
#)

###
### Make likelihood function
### 
#lik.func <- make_ikfdnm_likelihood(survival=~1, recruit=~1, detection=~1, kf_survival_effects=NULL, fixed_list=fixed_list, kf_data=kf_data, dnm_data=dnm_data, N_max=400)

###
### Optimize and obtain estimates and variance-covariance matrix
### 

#par_start=c(qlogis(0.95), log(4), qlogis(0.5))
#par_start=rep(0,3)
#lik.func(par_start)


#mle=optim(par_start, lik.func, method="BFGS", control=list(REPORT=1, trace=1), hessian=TRUE)
#par=mle$par
#se = sqrt(diag(2*solve(mle$hessian)))


###
### Estimates and 95% CI
###

# omega
#cat("Survival: \n")
#cat(plogis(par[1]), "(", plogis(par[1]-2*se[1]), ",",plogis(par[1]+2*se[1]), ")\n")

# gamma
#cat("Recruitment rate:")
#cat(exp(par[2]), "(", exp(par[2]-2*se[2]), ",",exp(par[2]+2*se[2]), ")\n")

# p
#cat("Detection prob.:\n")
#cat(plogis(par[3]), "(", plogis(par[3]-2*se[3]), ",",plogis(par[3]+2*se[3]), ")\n")


### SURVIVAL ANALYSIS USING SURVIVAL PACKAGE ###

#install.packages("survival")
library(survival)

##get time elapsed from first to last date as numeric (capture to last resight) and fate (in this case 1 is dead and 0 is alive)
surv.dat<-dim(0)
for (j in 1:length(unique(kf.dat$id))) {
  id.temp<-unique(kf.dat$id)[j]
  bday.temp<-subset(kf.dat, id==id.temp & Type=="capture" & Age==4)$Date
  last.temp<-max(kf.dat$Date[which(kf.dat$id==id.temp)])
  if (length(bday.temp)==0) {next} ##if the bird was not captured as a chick or doesn't have a capture date then move to next bird
  dur.temp<-as.numeric(last.temp-bday.temp)
  fate.temp<-subset(kf.dat, id==id.temp & Date==last.temp)$CH
  fate.surv.temp<-ifelse(fate.temp==1, 0, 1) ##if fate is 0, then make it 1, else make it 0
  
  #if (dur.temp==0) {next} ##if the bird was only observed at tagging, drop it
  #if (dur.temp==0) {dur.temp<-7; fate.surv.temp<-1} ##if bird was never observed after tagging, assume it died within a week
  ##default assumption is that bird was observed alive on day 0
  year.temp<-as.numeric(format(bday.temp, "%Y"))
  brood.temp<-subset(kf.dat, id==id.temp)$Nest.Group[1]
  pond.temp<-subset(kf.dat, id==id.temp & Date==bday.temp & Type=="capture")$Location
  
  if (dur.temp>31) {dur.temp<-31; fate.surv.temp<-0} ##if the bird was observed after fledging, assume it was alive at fledging and make its last observation at 31 days
  

  surv.dat<-rbind(surv.dat, data.frame(id=id.temp, year=year.temp, duration=dur.temp, fate=fate.surv.temp, pond=pond.temp, brood=brood.temp))
}

##SUBSET DATA

surv.dat<-subset(surv.dat, duration == 31 | fate == 1) ##use data where the bird either died or was tracked a full 31 days
#surv.dat<-subset(surv.dat, is.na(brood)==F) ## take only birds with known broods
surv.dat<-subset(surv.dat, pond %in% c("E14","E6B","E8", "E16B")) ##take only certain locations
#surv.dat<-subset(surv.dat, year %in% c(2016))

##SURVIVAL MODEL
surv.object <- Surv(time = surv.dat$duration, event = surv.dat$fate)

#surv.fitted.Motulsky <- survfit(surv.object ~ 1, conf.type = "log-log")
#plot(surv.fitted.Motulsky)

surv.fitted.default <- survfit(surv.object ~ 1)
plot(surv.fitted.default, xlab="Days", ylab="Proportion surviving")

##plot censored individuals
#plot(surv.fitted.default, mark.time = T, mark=16)

print(str_c("The proportion of chicks that survive to day ", last(surv.fitted.default$time)," is ", round(last(surv.fitted.default$surv),2), " (", round(last(surv.fitted.default$upper),2), ", ", round(last(surv.fitted.default$lower),2),")"))

par(mfrow=c(2,2), mar=c(4,4,3,1))
for (j in 1:length(unique(surv.dat$pond))) {
  pond.temp<-unique(surv.dat$pond)[j]
  dat.temp<-subset(surv.dat, pond==pond.temp)
  
  surv.object <- Surv(time = dat.temp$duration, event = dat.temp$fate)
  surv.fitted.default <- survfit(surv.object ~ 1)
  plot(surv.fitted.default, xlab="Days", ylab="Proportion surviving", main=pond.temp)
  text(x=5, y=0.2, labels = str_c("n alive = ", length(which(dat.temp$fate==0))))
  text(x=5, y=0.1, labels = str_c("n dead = ", length(which(dat.temp$fate==1))))
  mtext(text = print(str_c("Prop survive to day ", last(surv.fitted.default$time)," is ", round(last(surv.fitted.default$surv),2), " (", round(last(surv.fitted.default$upper),2), ", ", round(last(surv.fitted.default$lower),2),")")), side=3)
}
par(mfrow=c(1,1))

##CONS: assumes perfect detection rate. very sensitive to not detecting deaths. Curve only drops if death is detected, so need to make these plots for birds with confident fates

##can compare survival curves across sites and years: https://rpubs.com/brouwern/MotulskyCh5

##this approach is good for E14; good tracking of broods over time

###  N-MIXTURE MODEL  ###
head(uf.dat.sum)

#install.packages("jointNmix")
library(jointNmix)

##SIMULATION
## simulating data with negative binomial latent abundances
R <- 10 # sites
T <- 10 # time occasions
lambda <- 5 # abundance parameter
p <- .3 # probability of detection
phi <- 1 # dispersion parameter
set.seed(1234); Ni <- rnbinom(R, mu=lambda, size=phi) # latent abundances
y <- matrix(0, ncol=T, nrow=R)
set.seed(1234); for(i in 1:R) y[,i] <- rbinom(T, Ni, p) # observed abundances

##y is a matrix. rows are sites. columns are time. values are counts
##REAL DATA
data<-subset(uf.dat.sum, year==2017, select=c(group, survey, n))
y<-spread(data, survey, n)
##excuse missing data in last column
y<-y[,1:(ncol(y)-1)]
##remove sites that don't have all data
y<-y[complete.cases(y),]
groups.temp<-y$group
y<-as.matrix(subset(y, select=-group))
R<-nrow(y)
T<-ncol(y)

## fitting the Poisson N-mixture model
#fitp <- Nmix(y, Xp=cbind(rep(1, R*T)), Xl=cbind(rep(1, R)), mixture="P", K=25)
fitp <- Nmix(y, Xp=cbind(rep(1, R*T)), Xl=cbind(rep(1, R)), mixture="P", K=max(y,na.rm=T)*5)

## fitting the negative binomial N-mixture model
#fitnb <- Nmix(y, Xp=cbind(rep(1, R*T)), Xl=cbind(rep(1, R)), mixture="NB", K=25)
fitnb <- Nmix(y, Xp=cbind(rep(1, R*T)), Xl=cbind(rep(1, R)), mixture="NB", K=max(y,na.rm=T)*5)

## likelihood-ratio test between Poisson and negbin models
anova(fitp, fitnb)

## comparing using AIC
lapply(list(fitp, fitnb), AIC)

## conditional posterior probability functions for abundances
plot(fitnb, posterior = TRUE)

##p is probability of detection
##lambda is abundance parameter
##theta is dispersion parameter
fitnb
##https://rdrr.io/cran/jointNmix/man/Nmix.html

cbind(as.character(groups.temp),1:length(groups.temp))

##shows the numbers relative to the maximum number that a site can sustain

## GET POP SIZE ESTIMATES FOR EACH GROUP AND YEAR ##
fits<-list()
h<-0
fit.est<-dim(0)
#for (j in 10:length(unique(uf.dat.sum$group))) {
  for (i in 1:length(unique(uf.dat.sum$year))) {
    #group.temp<-unique(uf.dat.sum$group)[j]
    year.temp<-unique(uf.dat.sum$year)[i]
    #data.temp<-subset(uf.dat.sum, group==group.temp & year==year.temp, select=c(group, survey, n))
    data.temp<-subset(uf.dat.sum, year==year.temp, select=c(group, survey, n))
    if (nrow(data.temp)==0) {next}
    
    y<-spread(data.temp, survey, n)
    if(ncol(y)==2) {next} ##move to next if only one count date
    ##excuse missing data in last column
    y<-y[,1:(ncol(y)-1)]
    ##remove sites that don't have all data
    y<-y[complete.cases(y),]
    if (nrow(y)==1) {next} ##move to next if only one site
    groups.temp<-y$group
    y<-as.matrix(subset(y, select=-group))
    R<-nrow(y)
    T<-ncol(y)
    
    fitnb <- Nmix(y, Xp=cbind(rep(1, R*T)), Xl=cbind(rep(1, R)), mixture="NB", K=max(y,na.rm=T)*2)
    h<-h+1
    fits[[h]]<-fitnb
    #print(fitnb)
    
    fit.est<-rbind(fit.est, data.frame(year=rep(year.temp, length(groups.temp)), group=as.character(groups.temp), group.num=c(1:length(groups.temp)), list.item=rep(h,length(groups.temp)), Ni=getranef.uniNmix(fitnb)$Ni1))
  }
#} 

plot(fits[[1]], posterior=T)

getranef.uniNmix(fits[[1]])
coef(fits[[1]])

##plot the estimates
fig <- ggplot(data = fit.est, aes(x = year, y = Ni))
fig <- fig + geom_point(color="blue")
fig <- fig + geom_line(aes(color="model"))
fig <- fig + geom_line(data=subset(uf.dat.sum, group %in% c("Eden Landing", "Warm Springs", "Ravenswood", "Alviso")) %>% group_by(year, group) %>% mutate(n.max=max(n, na.rm=T)) %>% data.frame(), aes(x=year, y=n.max, color="N max", linetype="N max"))
fig <- fig + geom_line(data=subset(uf.dat.sum, group %in% c("Eden Landing", "Warm Springs", "Ravenswood", "Alviso")) %>% group_by(year, group) %>% mutate(n.mean=mean(n, na.rm=T)) %>% data.frame(), aes(x=year, y=n.mean, color="N mean"))
fig <- fig + geom_point(data = subset(uf.dat.sum, group %in% c("Eden Landing", "Warm Springs", "Ravenswood", "Alviso")), aes(x = year, y = n), shape=1, color="black")

fig <- fig + facet_wrap(~group, scales = "free_y")
fig <- fig + labs(color="") + scale_color_manual(values = c("blue", "black", "black")) 
fig <- fig + scale_linetype_manual(values = c("dashed", "dashed", "solid"))
#fig <- fig + scale_shape_manual(values = c(2))
fig <- fig + theme_bw()
fig <- fig + scale_x_continuous(breaks=seq(min(uf.dat.sum$year),max(uf.dat.sum$year), 2), expand=c(0.05,0))
fig

##get detection probability by dividing mean count by estimated pop size??

## ANOTHER N MIXTURE MODEL ##
##https://cran.r-project.org/web/packages/unmarked/unmarked.pdf pcount function

##read unmarked paper