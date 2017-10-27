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

library(ggplot2)

library(stringr) ##required for string processing

library(RODBC) ##required to pull data from access databases

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
sim.dat<-sim_data(num_kf = 5, num_years = 27, num_surveys = 10, recruit_rate = 5, init_rate = 15, survival_dnm = .90, survival_kf = .90, perfect_survey_rate = .3, detection = .7)

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


##ADD CHICKS PRESUMED DEAD
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

qry<-"SELECT year(i.Date) AS year, i.Date, s.SurveyID AS survey, 0 AS perfect_survey, 0 AS R, s.TotalObserved AS n, p.Complex AS [complex], p.PondNumber
  FROM ([SNPL SURVEY] AS s
  LEFT OUTER JOIN [SNPL Observers Info] AS i ON s.SurveyID = i.SurveyID)
  LEFT OUTER JOIN LookupPond AS p ON s.PondNumber = p.PondNumber
  WHERE p.Complex <> 'Error With Record' "

uf.dat<-sqlQuery(con, qry); head(uf.dat) ##import the queried table

##when finished with db, close the connection
odbcCloseAll()


##clean up data
##remove nonsensical year
#uf.dat.sum<-subset(uf.dat.sum, year > 1980)

##restrict surveys to "window surveys" when only reproductive individuals are present (ie exclude counts of migratory birds)
window.dates<-read.csv("SNPL Breeding Window Dates.csv")
window.dates$start.date<-as.Date(window.dates$start.date, "%m/%d/%Y")
window.dates$end.date<-as.Date(window.dates$end.date, "%m/%d/%Y")

uf.dat.win<-dim(0) ##counts of untagged indivs during the "window surveys" for each year
for (j in 1:nrow(window.dates)) {
  year.temp<-window.dates$Year[j]
  dat.temp<-subset(uf.dat, year==year.temp)
  win.temp<-window.dates[j,]
  out.temp<- subset(dat.temp, as.Date(dat.temp$Date) >= win.temp$start.date & as.Date(dat.temp$Date) <= win.temp$end.date)
  
  if (nrow(out.temp)>0) {
    ##for survey dates that apply to all locations
    if (win.temp$locations=="all") {
      uf.dat.win<-rbind(uf.dat.win, out.temp)
    }
    
    ##for survey windows that apply to select locations
    if (win.temp$locations=="Hayward") {
      uf.dat.win<-rbind(uf.dat.win, subset(out.temp, complex=="Hayward"))
    }
    
    if (win.temp$locations=="Hayward, Napa") {
      uf.dat.win<-rbind(uf.dat.win, subset(out.temp, complex %in% c("Hayward", "Napa")))
    }
    
    if (win.temp$locations=="all-Hayward") {
      uf.dat.win<-rbind(uf.dat.win, subset(out.temp, complex !=  "Hayward"))
    }
    
    if (win.temp$locations=="all-Hayward, Napa") {
      uf.dat.win<-rbind(uf.dat.win, subset(out.temp, complex != "Hayward" & complex != "Napa"))
    }
  }
}

##FIGURE OUT WHAT LOCATIONS WE SHOULD MAKE THE "GROUP" VARIABLE
#uf.dat$group<-uf.dat$complex

##use Ben's groups
uf.dat.win<-dplyr::left_join(x=uf.dat.win, y=subset(pond.groups, select=c(Pond, Pond.Group)), by = c("PondNumber" = "Pond"))
uf.dat.win$group<-uf.dat.win$Pond.Group

##get summary n for a given date and location
uf.dat.sum<-uf.dat.win %>% group_by(year, group) %>% mutate(year.n=sum(n)) %>% data.frame()
##order by group, then time
uf.dat.sum<-uf.dat.sum[order(uf.dat.sum$group, uf.dat.sum$year),]; head(uf.dat.sum)
##rename survey to be all 1, since we will represent each yearly count during the window survye as the one and only survey
uf.dat.sum$survey<-rep(1, nrow(uf.dat.sum))
##remove n and remove duplicates; rename survey.n to n
uf.dat.sum<-subset(uf.dat.sum, select = -c(n,Date)); uf.dat.sum<-unique(uf.dat.sum)
names(uf.dat.sum)[length(uf.dat.sum)]<-"n"
head(uf.dat.sum)


#dnm_data<-subset(uf.dat.sum, group == "Eden Landing")
dnm_data<-uf.dat.sum
kf_data<-kf.dat

##visualize the data
vis <- ggplot(data = dnm_data, aes(x = year, y = n))
vis <- vis + geom_point()
vis <- vis + facet_grid(facets = dnm_data$group~., scales = "free_y")
vis

##visual data for known-fate birds
j<-0
j<-j+1
fplot <- ggplot(data = subset(kf_data, subset = id == unique(kf_data$id)[j]), aes(x = Date, y = CH))
fplot <- fplot + geom_point()
fplot

###
### Create list of fixed parameters
###

fixed_list = list(
  recruit=ifelse(dnm_data$year!=1 & dnm_data$survey==1, NA, 0),
  detection=ifelse(dnm_data$perfect==1, 1, NA)
)

###
### Make likelihood function
### 
lik.func <- make_ikfdnm_likelihood(survival=~1, recruit=~1, detection=~1, kf_survival_effects=NULL, fixed_list=fixed_list, kf_data=kf_data, dnm_data=dnm_data, N_max=400)

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


### SURVIVAL ANALYSIS USING SURVIVAL PACKAGE ###

#install.packages("survival")
library(survival)

##get time elapsed from first to last date as numeric (capture to last resight) and fate (in this case 1 is dead and 0 is alive)
surv.dat<-dim(0)
for (j in 1:length(unique(kf.dat$id))) {
  id.temp<-unique(kf.dat$id)[j]
  bday.temp<-subset(kf.dat, id==id.temp & Type=="capture" & Age==4)$Date
  last.temp<-max(kf.dat$Date[which(kf.dat$id==id.temp)])
  if (length(bday.temp)==0) {next}
  dur.temp<-as.numeric(last.temp-bday.temp)
  if (dur.temp==0) {next}
  fate.temp<-subset(kf.dat, id==id.temp & Date==last.temp)$CH
  fate.surv.temp<-ifelse(fate.temp==1, 0, 1) ##if fate is 0, then make it 1, else make it 0
  
  if (is.na(subset(kf.dat, id==id.temp)$Nest.Group[1])) {next} ##fit the survival curve for known broods only
  if (dur.temp>31) {next} ##fit the survival curve only until the age of fledging
  
  surv.dat<-rbind(surv.dat, data.frame(duration=dur.temp, fate=fate.surv.temp))
}

surv.object <- Surv(time = surv.dat$duration, event = surv.dat$fate)

surv.fitted.Motulsky <- survfit(surv.object ~ 1, conf.type = "log-log")
plot(surv.fitted.Motulsky)

surv.fitted.default <- survfit(surv.object ~ 1)
plot(surv.fitted.default, xlab="Days", ylab="Proportion surviving")

##plot censored individuals
#plot(surv.fitted.default, mark.time = T, mark=16)

print(str_c("The proportion of chicks that survive to day ", last(surv.fitted.default$time)," is ", round(last(surv.fitted.default$surv),2), " (", round(last(surv.fitted.default$upper),2), ", ", round(last(surv.fitted.default$lower),2),")"))

