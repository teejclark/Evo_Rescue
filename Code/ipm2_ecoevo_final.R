library(tidyverse)
library(nimble)
library(ggplot2)
library(MCMCvis)
library(lubridate)
library(IPMbook)

# BUILD ECO-EVO MODEL FOR PENGUINS

##################################################################################################################

#### LOAD FUNCTIONS ####

##nimblefunction
ddataBern <- nimbleFunction(
  run = function(x = double(1), prob = double(1), length = integer(), log = double()) {
    ll <- 0
    for(i in 1:length) {
      ll <- ll + dbinom(x[i], prob = prob[i], size = 1, log=TRUE)
    }
    returnType(double())
    if(log) return(ll) else return(exp(ll))
  }
)
rdataBern <- nimbleFunction(
  run = function(n = integer(), prob = double(1), length = integer()) {
    print('this should never run')
    ##x <- numeric(length)
    declare(x, double(1, length))
    returnType(double(1))
    return(x)
  }
)
ddataPois <- nimbleFunction(
  run = function(x = double(1), lambda = double(1), length = integer(), log = double()) {
    ll <- 0
    for(i in 1:length) {
      ll <- ll + dpois(x[i], lambda[i], log=TRUE)
    }
    returnType(double())
    if(log) return(ll) else return(exp(ll))
  }
)
rdataPois <- nimbleFunction(
  run = function(n = integer(), lambda = double(1), length = integer()) {
    print('this should never run')
    ##x <- numeric(length)
    declare(x, double(1, length))
    returnType(double(1))
    return(x)
  }
)

ddataNorm <- nimbleFunction(
  run = function(x = double(1), mean = double(1), tau = double(), length = integer(), log = double()) {
    sd <- 1/sqrt(tau)
    ll <- 0
    for(i in 1:length) {
      ll <- ll + dnorm(x[i], mean[i], sd=sd, log=TRUE)
    }
    returnType(double())
    if(log) return(ll) else return(exp(ll))
  }
)
rdataNorm <- nimbleFunction(
  run = function(n = integer(), mean = double(1), tau = double(), length = integer()) {
    print('this should never run')
    ##x <- numeric(length)
    declare(x, double(1, length))
    returnType(double(1))
    return(x)
  }
)

djdistri <- nimbleFunction(
  run = function(x = double(1), rate = double(1), length = integer(), log = integer()){
    rate1<-rate[1:length]/sum(rate[1:length])
    x1<-x[1:length]/sum(x[1:length])
    prob<-0
    for(i in 1:length) {
      prob<-prob+min(x1[i],rate1[i])
    }
    returnType(double(0))
    logProb <- log(prob)
    if( log ) return(logProb)
    else return(prob)
  }
)
rjdistri <- nimbleFunction(
  run = function(n = integer(), rate = double(1), length = integer()) {
    print('this should never run')
    ##x <- numeric(length)
    declare(x, double(1, length))
    returnType(double(1))
    return(x)
  }
)

dcmr <- nimbleFunction(
  run = function(x = double(1), P = double(1), nind = integer(),  log  = integer()){
    # Likelihood 
    LL<-sum(x[1:nind]*log(P[1:nind]))
    returnType(double(0))
    if( log ) return(LL)
    else return(exp(LL))
  }
)
rcmr <- nimbleFunction(
  run = function(n = integer(), P = double(1), nind = integer()) {
    print('this should never run')
    ##x <- numeric(length)
    declare(x, double(1, nind))
    returnType(double(1))
    return(x)
  }
)


registerDistributions(list(
  djdistri = list(
    BUGSdist = 'djdistri(rate, length)',
    types    = c('value = double(1)', 'rate = double(1)', 'length = integer()')
  ),
  dcmr = list(
    BUGSdist = 'dcmr(P, nind)',
    types    = c('value = double(1)', 'P = double(1)', 'nind = integer()')
  ),
  
  ddataBern = list(
    BUGSdist = 'ddataBern(prob, length)',
    types    = c('value = double(1)', 'prob = double(1)', 'length = integer()')
  ),
  ddataPois = list(
    BUGSdist = 'ddataPois(lambda, length)',
    types    = c('value = double(1)', 'lambda = double(1)', 'length = integer()')
  ),
  ddataNorm = list(
    BUGSdist = 'ddataNorm(mean, tau, length)',
    types    = c('value = double(1)', 'mean = double(1)', 'tau = double()', 'length = integer()')
  )
))

#OTHER FUNCTIONS
eq.vec <- function(x,a) {
  dec=T
  k=0
  while(k<length(x)&dec==T){
    k=k+1
    if(x[k]!=a[k]){dec=F}
  }
  return(dec)
}

get.first <- function(x) min(which(x==1))
get.nobs<- function(x) length(which(x==1))
get.last <- function(x) max(which(x==1))

##################################################################################################################

#### LOAD DATA ####

# LOAD POPULATION DATA #

# REPRODUCTION: NUMNESTS = number of nests, numfledged = number fledged!
repro <- read.csv("ReproSuccess_9.21.21.csv", header = F)
repro$V1[1] <- "1984" # fixed!
colnames(repro) <- c("bookyear", "numnests", "numeggs", "clutchsize", "numhatch", "meanhatch", "numfledged", "rs")
repro$bookyear <- as.Date(as.character(repro$bookyear), "%Y")
repro <- repro %>% arrange(bookyear)

# number of surveyed nests
R <- as.numeric(c(repro$numnests[1:28], NA, repro$numnests[29:36]))

# number of fledges
J <- as.numeric(c(repro$numfledged[1:28], NA, repro$numfledged[29:36]))

rz <- c(R[1:28],140,R[30:37])

# POPULATION DATA (using Ginger's normal adjustment, not sliding scale...)
pop <- read.csv("StakeSurveyTOMACT22_9.21.21.csv", header = F)
pop[1,1] <- 1987 
colnames(pop) <- c("bookyear","00N01E", "00N02E", "00N03E", "00N04E", "00N05E", "00N06E", "00N07E", "00N08E",
                   "00N09E", "00N12E", "00N14E", "00N15E", "00N16E", "00N17E", "00N18E", "01N14E", "01S14E",
                   "02N14E", "02S14E", "03S14E", "04S14E", "05S14E", "total")
pop$density <- pop$total/22 # number of 100 meter survey areas = 22

# convert by 1.36
pop$blah <- c(pop$density[1:28]/1.36, pop$density[29:33])
pop$popest2 <- pop$blah*35240 # total habitat reported in the paper...(Table 1, 2012 survey)

y <- c(rep(NA, 4), pop$popest2[1:24], NA, pop$popest2[25:32])

# LOAD SURVIVAL DATA #
combo_size2 <- read.csv("juvenile_survival_by_ind_12422.csv")

# survival model functions
known.state.cjs <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}

cjs.init.z <- function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2 <- max(which(ch[i,]==1))
    ch[i,f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}

combo_size3 <- combo_size2

# capture histories
CH <- combo_size3 %>%
  dplyr::select(-sex, -bd, -bl, -f, -fl, -w,
                -pc1overall, -pc2overall, -recruit, -total) %>%
  data.matrix(.)
CH2 <- CH[,-c(1:2)]
#CH2 <- cbind(rep(0,nrow(CH)),CH) # add in 1982! (none marked)

# create vector with first occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH2, 1, get.first)

# last occasion = NOTE: this might not matter
last <- which(f==ncol(CH2))
f <- f[-last] # remove individuals marked on the last occasion, since they provide no data! (removes 77)
CH2 <- CH2[-last,]
combo_size3 <- combo_size3[-last,]

# LOAD INHERITANCE DATA #
family <- read.csv("HereditaryInfo_1.25.22.csv") # mother/father - offspring info
colnames(family) <- c("pengid","fatherid", "motherid")

juvsize <- read.csv("juvenile_survival_by_ind_12422.csv") # size of each fledge
juvsize$pengid <- as.character(juvsize$pengid)
juvsize2 <- juvsize[,-c(9:46)] # clean up code and extraneous variables

adultsize <- read.csv("adult_survival_by_ind_12522.csv") # size of each adult
adultsize$pengid <- as.character(adultsize$pengid)
adultsize <- adultsize[,-c(10:47)] # clean up code and extraneous variables

herit <- inner_join(juvsize2, family, by = "pengid") # join family table to juveniles = 6419 matches!

heritfather <- inner_join(herit, juvsize2, by = c("fatherid" = "pengid")) # join family table to fathers = 5188 matches

heritmother <- inner_join(herit, juvsize2, by = c("motherid" = "pengid"))

heritall <- full_join(heritfather, heritmother, by = "pengid") %>%
  rowwise() %>%
  mutate(averagePC1 = mean(c(pc1overall.y.x,pc1overall.y.y), na.rm = T),
         averagePC2 = mean(c(pc2overall.y.x,pc2overall.y.y), na.rm = T),
         averagebd = mean(c(bd.y.x,bd.y.y), na.rm = T),
         averagebl = mean(c(bl.y.x,bl.y.y), na.rm = T),
         averagefl = mean(c(fl.y.x,fl.y.y), na.rm = T),
         averagef = mean(c(f.y.x,f.y.y), na.rm = T),
         averagew = mean(c(w.y.x,w.y.y), na.rm = T),
         pc1overall.x = mean(c(pc1overall.x.x, pc1overall.x.y), na.rm = T))


# LOAD ENVIRONMENTAL DATA #
# OCEAN - MIGRATION - SSTA
ssta_migration <- read.csv("ssta.PC1.migration.11.12.21.csv") %>%
  mutate(bookyear = year)
ssta_migration <- ssta_migration[-1,]


# OCEAN - BREEDING - SSTA
ssta_breeding <- read.csv("ssta_summarized.12.2.21.csv")
ssta_breeding <- ssta_breeding[-1,]

# load ppt
weather <- read.csv("Site Weather.csv")

# replace NAs in precipitation
weather1 <- weather %>% 
  mutate(Precipitation = replace(Precipitation, which(is.na(Precipitation)), 0.0))

# convert character dates to actual dates
weather1$WeathDate <- as.Date(weather1$WeathDate)

# calculate total precipitation between Oct 15 and Dec 15
rain_60 <-  weather1 %>%
  filter(WeathDate >= as.Date(paste(year(WeathDate), 10, 15, sep = "-")),
         WeathDate <= as.Date(paste(year(WeathDate), 12, 15, sep = "-"))) %>%
  group_by(BookYear) %>%
  filter(!is.na(Precipitation)) %>%
  summarize(amt_rain = sum(Precipitation)) %>%
  rename(bookyear = BookYear)

rain_60 <- rbind(rain_60[1:28,], c(2011, 36.6), rain_60[29:36,])

# load temp
num_days_25 <- weather1 %>% group_by(BookYear) %>%
  filter(!is.na(MaxTemp)) %>%
  summarize(n_days = n(),
            n_gt25 = sum(MaxTemp > 25),
            p_gt25 = n_gt25/n_days) %>%
  rename(bookyear = BookYear)

num_days_25 <- rbind(num_days_25[1:28,], c(2011, 0.383), num_days_25[29:36,])

##################################################################################################################

#### MANIPULATE DATA ####

# TRAIT DATA #
# set up breakpoints, meshpoints, etc!
nn <- 50 # edit to run model quickly
minsize <- -5 # minimum possible size 
maxsize <- 5 # maximum possible size 
L <- 1.1*minsize; U <- 1.1*maxsize
b <- L+c(0:nn)*(U-L)/nn # boundary points (edges of boxes)
meshpoints <- 0.5*(b[1:nn]+b[2:(nn+1)]) # meshpoints ()


# create initial population data by trait (inipopN1 and inipopNad)
sdj <- sd(combo_size3$pc1overall[combo_size3$X1983==1]) # standard deviation on trait for juveniles (2529 individuals)
sda <- sd(combo_size3$pc1overall[combo_size3$recruit == 1]) # standard deviation on trait for surviving adults
ini1 <- rnorm(1000, mean(-1*combo_size3$pc1overall[combo_size3$X1983==1]), sdj)
ini2 <- rnorm(1000, mean(-1*combo_size3$pc1overall[combo_size3$recruit == 1]), sda) 
ini3 <- rnorm(1000, mean(-1*combo_size3$pc1overall[combo_size3$recruit == 1]), sda)
iniP <- rnorm(1000, mean(-1*combo_size3$pc1overall[combo_size3$recruit == 1]), sda)
iniB <- rnorm(1000, mean(-1*combo_size3$pc1overall[combo_size3$recruit == 1]), sda)
iniN <- rnorm(1000, mean(-1*combo_size3$pc1overall[combo_size3$recruit == 1]), sda)

inipopN1 <- t(as.matrix(table(cut(ini1,breaks=b))))
inipopN2 <- t(as.matrix(table(cut(ini2,breaks=b))))
inipopN3 <- t(as.matrix(table(cut(ini3,breaks=b))))
inipopNP <- t(as.matrix(table(cut(iniP,breaks=b))))
inipopNB <- t(as.matrix(table(cut(iniB,breaks=b))))
inipopNN <- t(as.matrix(table(cut(iniN,breaks=b))))

inipopN1 <- inipopN1/sum(inipopN1)*25000 
inipopN2 <- inipopN2/sum(inipopN2)*25000 
inipopN3 <- inipopN3/sum(inipopN3)*25000 
inipopNP <- inipopNP/sum(inipopNP)*25000
inipopNB <- inipopNB/sum(inipopNB)*250000
inipopNN <- inipopNN/sum(inipopNN)*50000

# create ini matrix
inivalsN1 <- matrix(rep(inipopN1/sum(inipopN1)*25000, length(y)), ncol = length(y))
inivalsN2 <- matrix(rep(inipopN2/sum(inipopN2)*25000, length(y)), ncol = length(y))
inivalsN3 <- matrix(rep(inipopN3/sum(inipopN3)*25000, length(y)), ncol = length(y))
inivalsNP <- matrix(rep(inipopNP/sum(inipopNP)*25000, length(y)), ncol = length(y))
inivalsNB <- matrix(rep(inipopNB/sum(inipopNB)*250000, length(y)), ncol = length(y))
inivalsNN <- matrix(rep(inipopNN/sum(inipopNN)*50000, length(y)), ncol = length(y))

# ENVIRONMENTAL DATA FOR INHERITANCE AND SURVIVAL #
# get occasion of first marking and add it to 1982 to get year of marking
ssta_breeding_lag <- ssta_breeding
ssta_breeding_lag$bookyear <- ssta_breeding_lag$bookyear+1

ssta_migration_lag <- ssta_migration
ssta_migration_lag$bookyear <- ssta_migration_lag$bookyear-1

juvsize_env <- juvsize %>%
  mutate(tag = apply(juvsize[,9:45], 1, get.first),
         bookyear = 1983+tag) %>%
  inner_join(ssta_breeding, by = "bookyear") %>%
  inner_join(ssta_migration, by = "bookyear") %>% 
  semi_join(heritall, by = "pengid") # now join with heritability data

# combine
herit_env <- inner_join(juvsize_env, heritall, by = "pengid")

# Survival (41,887 individuals)
surv_env <- inner_join(ssta_breeding, ssta_migration, by = "bookyear") %>%
  mutate(p_gt25 = scale(num_days_25$p_gt25), amt_rain = scale(rain_60$amt_rain)) %>%
  pivot_wider(names_from = bookyear, values_from = c(av_s, mean.2, p_gt25, amt_rain))

surv_env2 <- matrix(rep(c(ssta_breeding$av_s, ssta_migration$mean.2, scale(num_days_25$p_gt25), scale(rain_60$amt_rain)), nrow(CH2)),
                    nrow = nrow(CH2),
                    ncol = length(c(ssta_breeding$av_s, ssta_migration$mean.2, num_days_25$p_gt25, rain_60$amt_rain)),
                    byrow = T)
colnames(surv_env2) <- c(colnames(surv_env[,12:85]),"mean.2_2020", colnames(surv_env[,86:159]))


# traits
surv_traits <- matrix(rep(combo_size3$pc1overall, length(ssta_breeding$av_s)),
                      nrow = length(combo_size3$pc1overall),
                      ncol = length(ssta_breeding$av_s),
                      byrow = F)

for (i in which(f>1)){
  surv_traits[i,1:(f[i]-1)]=0
}


##################################################################################################################

#### IPM2 ####

# number of individuals per unique capture history
surv_num <- 100 
b_num <- L+c(0:surv_num)*(U-L)/surv_num # boundary points
meshpoints_num <- 0.5*(b_num[1:surv_num]+b_num[2:(surv_num+1)]) # midpoints
ints_num <- findInterval(-1*combo_size3$pc1overall, b_num, all.inside=T)
trait_mids_num <- (b_num[ints_num] + b_num[ints_num + 1])/2

# calculate new survival traits
surv_traits_num <- matrix(rep(trait_mids_num, length(ssta_breeding$av_s)),
                          nrow = length(combo_size3$pc1overall),
                          ncol = length(ssta_breeding$av_s),
                          byrow = F)
for (i in which(f>1)){
  surv_traits_num[i,1:(f[i]-1)]=0
}

RCH_traits <- cbind(CH2, surv_traits_num)
uCH_traits <- unique(RCH_traits) 

nbru <- readRDS("nbru_100traits.rds")

CH4 <- uCH_traits[,1:37]
surv_traits_4 <- uCH_traits[,38:74]

f <- apply(CH4, 1, get.first)
last <- apply(CH4, 1, get.last)
nobsu <- apply(CH4, 1, get.nobs)

max_age <- 3
x <- createAge(f, nyears=ncol(CH4)+1, mAge = max_age+2, age = rep(1,length(f)))

# now convert to three age classes, corresponding to three survival probabilities: juvies, prebreeders, adults
#x <- x-1 # convert to 0s, 1s, etc.
x[x==3] <- 2; x[x==4] <- 2 
x[x==5] <- 3 

# put data in correct order for the survival model
nobso <- nobsu
nobso[nobsu>1]=0
CH4_alt <- CH4[order(nobso),]
f_alt <- f[order(nobso)]
nbru_alt <- nbru[order(nobso)]
last_alt <- last[order(nobso)]
x_alt <- x[order(nobso),]
surv_traits_alt <- surv_traits_4[order(nobso),]
nobsu2 <- nobsu[order(nobso)]

lasto <- rep(1,length(last))
lasto[last_alt==ncol(CH4)]=0
CH4_alt <- CH4_alt[order(lasto),]
f_alt <- f_alt[order(lasto)]
nbru_alt <- nbru_alt[order(lasto)]
last_alt <- last_alt[order(lasto)]
nobsu2 <- nobsu2[order(lasto)]
x_alt <- x_alt[order(lasto),]
surv_traits_alt <- surv_traits_alt[order(lasto),]
surv_env3 <- surv_env2[1:nrow(surv_traits_alt),]
CH4_alt[CH4_alt==0]=2 # convert 0s to 2s


# Bayesian code
IPM2code_peng1 <- nimbleCode({
  
  ####
  # 1. Define priors for the parameters
  ####
  
  # initial population sizes by trait
  N1[1:nn, 1] <- inipopN1[1:nn]
  N2[1:nn, 1] <- inipopN2[1:nn]
  N3[1:nn, 1] <- inipopN3[1:nn]
  NP[1:nn, 1] <- inipopNP[1:nn]
  NB[1:nn, 1] <- inipopNB[1:nn]
  NN[1:nn, 1] <- inipopNN[1:nn]
  
  # mean demographic parameters
  l.mfec ~ dnorm(0, 0.34)
  l.mim ~ dnorm(0, 0.34)
  l.minh ~ dnorm(0, 0.34)
  
  # covariates
  for (i in 1:5){alpha[i] ~ dnorm(0, 0.34)} # inheritance 
  
  for (a in 1:max_age){
    l.mphi[a] ~ dnorm(0, 0.34)
    l.p[a] ~ dnorm(0, 0.34)
    for (i in 1:6){beta[i,a] ~ dnorm(0, 0.34)}} # survival
  
  # transition parameters
  gamma ~ dunif(0, 1) # probability of breeding for the first time
  delta ~ dunif(0, 1) # probability of staying a breeder
  epsilon ~ dunif(0, 1) # probability of staying a nonbreeder
  
  # variance fxns
  inh_var ~ dunif(0, 100) # variance of the inheritance function
  tau_I <- 1 / (inh_var^2)
  
  # observation error
  tauy <- pow(sigma.y, -2) # precision
  sigma.y ~ dunif(0, 200000) # sd (50000)
  
  ####
  # 2. Constrain parameters
  ####
  
  for (t in 1:(nyears-1)){
    log(fec[t]) <- l.mfec 
    log(omega[t]) <- l.mim
  }
  
  ####
  # 3. Derived parameters
  ####

  # # population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1]/Ntot[t]
    l.lambda[t] <- log(lambda[t])
  }

  # # geometric mean
  geomean.lambda <- exp((1/(nyears-1))*sum(l.lambda[1:(nyears-1)]))
  
  ####
  # 4. Likelihoods
  ####
  
  # likelihood for population count data
  # system process
  # inheritance by traits
  for (t in 1:(nyears-1)){ 
    for(trait in (1:nn)){ # nn = # of breakpoints in traits
      I[1:nn,trait,t] <- exp(-((X[1:nn] - (l.minh + alpha[1]*X[trait] + alpha[2]*SSTA[t] + alpha[3]*PLUME[t] + alpha[4]*X[trait]*SSTA[t] + alpha[5]*X[trait]*PLUME[t]))^2)/(2*inh_var*inh_var))/(sqrt(2*pi)*inh_var)
      Inorm[1:nn,trait,t] <- I[1:nn,trait,t]/(sum(I[1:nn,trait,t]) + 0.000001)
    }
    
    # survival by traits
    for (a in 1:max_age){
      logit(s[1:nn,t,a]) <- l.mphi[a] + beta[1,a]*X[1:nn] + beta[2,a]*SSTA[t] + beta[3,a]*PLUME[t] + beta[4,a]*X[1:nn]*SSTA[t] + beta[5,a]*X[1:nn]*PLUME[t] + beta[6,a]*TEMP[t]
    }
  }
  
  # population distribution
  # NOTE: adding in more pre-breeding stages and a nonbreeding stage
  for (t in 2:nyears){
    N1[1:nn,t] <- Inorm[1:nn,1:nn,t-1] %*% (0.5 * fec[t-1] * s[1:nn,t-1,1] * NB[1:nn,t-1]) # juveniles
    N2[1:nn,t] <- s[1:nn,t-1,2] * N1[1:nn,t-1] # 2nd years
    N3[1:nn,t] <- s[1:nn,t-1,2] * N2[1:nn,t-1] # 3rd years
    NP[1:nn,t] <- (s[1:nn,t-1,2] * N3[1:nn,t-1]) + (s[1:nn,t-1,2] * (1-gamma) * NP[1:nn,t-1])
    NB[1:nn,t] <- (s[1:nn,t-1,2] * gamma * NP[1:nn,t-1]) + (s[1:nn,t-1,3] * delta * NB[1:nn,t-1]) + (s[1:nn,t-1,3] * (1-epsilon) * NN[1:nn,t-1]) + (s[1:nn,t-1,3] * NB[1:nn,t-1] * omega[t-1]) # breeders + immigrants
    NN[1:nn,t] <- (s[1:nn,t-1,3] * (1-delta) * NB[1:nn,t-1]) + (s[1:nn,t-1,3] * epsilon * NN[1:nn,t-1])
  }
  
  # only observed breeding individuals
  for (t in 1:nyears){
    Ntot[t] <- sum(NB[1:nn,t])
  }
  
  # observation process
  for (t in 1:nyears){
    y[t] ~ T(dnorm(Ntot[t], tauy),0,)
  }
  
  # likelihood for individual data
  # inheritance by individual
  muI[1:Ndatainh] <- l.minh + alpha[1]*IX[1:Ndatainh] + alpha[2]*ISSTA[1:Ndatainh] + alpha[3]*IPLUME[1:Ndatainh] + alpha[4]*IX[1:Ndatainh]*ISSTA[1:Ndatainh] + alpha[5]*IX[1:Ndatainh]*IPLUME[1:Ndatainh]
  Idata[1:Ndatainh] ~ ddataNorm(muI[1:Ndatainh], tau_I, Ndatainh)
  
  # survival by individual
  for (i in 1:nind){
    for (a in f[i]:(n.occasions)){
      
    logit(phi[i,a]) <- l.mphi[age[i,a]] +
      beta[1,age[i,a]]*SX[i,a] +
      beta[2,age[i,a]]*SSSTA[i,a] +
      beta[3,age[i,a]]*SPLUME[i,a] +
      beta[4,age[i,a]]*SX[i,a]*SSSTA[i,a] +
      beta[5,age[i,a]]*SX[i,a]*SPLUME[i,a] +
      beta[6,age[i,a]]*STEMP[i,a]

    logit(p[i,a]) <- l.p[age[i,a]]
  }}
  
  #m-array survival
  # annual probabilities of individual capture histories
  for (i in 1:nind2){ # individuals that have been recaptured at least once 
    for (t in f[i]:(last[i]-1)) { # probabilities from first to last observation
      MO[i,t,1]<-phi[i,t]*p[i,t]
      MO[i,t,2]<-phi[i,t]*(1-p[i,t])
      p1[i,t]<-MO[i,t,ys[i,(t+1)]]
    }
  }
  for (i in (nindobs1+1):nind){ # individuals that have been recaptured before last observation
    # probability that they died the year after last observation
    p2[i,last[i]] <- (1-phi[i,last[i]])
    # probability that they died between last observation and the end of the study period
    for (t in (last[i]+1):(n.occasions)){
      p2[i,t] <- (1-phi[i,t])*prod(phi[i,last[i]:(t-1)]*(1-p[i,last[i]:(t-1)]))
    }
    # probability that they never die
    p3[i] <- prod(phi[i,last[i]:(n.occasions-1)]*(1-p[i,last[i]:(n.occasions-1)]))
    
  }
  # total capture-history probabilities
  for (i in 1:nindobs1){# individuals that have been recaptured in the last year
    P[i]<-prod(p1[i,f[i]:(last[i]-1)])}
  for (i in (nindobs1+1):nind2){ # individuals that have been recaptured at least once and not in the last year
    P[i]<-prod(p1[i,f[i]:(last[i]-1)])*(sum(p2[i,(last[i]):(n.occasions-1)])+p3[i])}
  for (i in (nind2+1):nind){ # individuals that have never been recaptured
    P[i]<-sum(p2[i,(last[i]):(n.occasions-1)])+p3[i]}
  # likelihood
  M[1:nind]~dcmr(P[1:nind],nind)
  
  # reproduction
  for (t in 1:(nyears-1)){
    J[t] ~ dpois(rho[t])
    rho[t] <- R[t]*fec[t]
  }
  
})

my.data <- list(y = y,
                X = meshpoints, 
                inipopN1 = as.vector(inipopN1), 
                inipopN2 = as.vector(inipopN2),
                inipopN3 = as.vector(inipopN3),
                inipopNP = as.vector(inipopNP),
                inipopNB = as.vector(inipopNB),
                inipopNN = as.vector(inipopNN),
                Idata = -1*herit_env$pc1overall.x, 
                M = nbru_alt, 
                J = J, 
                R = rz)

my.constants <- list(nyears = length(y),
                     ys = as.matrix(CH4_alt), 
                     max_age = max_age, 
                     f = f_alt, 
                     last = last_alt, 
                     nind = dim(CH4_alt)[1], 
                     nindobs1 = length(which(last_alt==ncol(CH4_alt))), 
                     nind2 = length(which(nobsu2>1)), 
                     n.occasions = dim(CH4_alt)[2], 
                     pi = pi, 
                     nn = nn, 
                     Ndatainh = nrow(herit_env), 
                     age = x_alt, 
                     #RAIN = as.vector(scale(rain_60$amt_rain)), 
                     TEMP = as.vector(scale(num_days_25$p_gt25)), 
                     SSTA = as.vector(ssta_breeding$av_s[1:37]), 
                     PLUME = as.vector(c(ssta_migration$mean.2[2:38]))) 

# create initial values for Ntot
blah <- y
blah[1:5] <- seq(280000,250000,length.out=5)
blah[29] <- 180000

# inits
initial.values <- function(){list(
  l.mphi = rnorm(max_age, 0, 0.5),
  l.p = rnorm(max_age, 0, 1),
  l.mfec = rnorm(1, 0, 0.5),
  l.mim = rnorm(1, 0, 0.5),
  l.minh = rnorm(1, 0, 0.5),
  inh_var = runif(1, 0, 100), 
  alpha = rnorm(5, 0, 1), 
  beta = array(rnorm(18,0,1), dim = c(6,max_age)), 
  gamma = runif(1, 0, 1),
  delta = runif(1, 0, 1),
  epsilon = runif(1, 0, 1),
  sigma.y = runif(1,0,50000),
  N1 = matrix(rpois(length(inivalsN1), inivalsN1), ncol = length(y)), 
  N2 = matrix(rpois(length(inivalsN2), inivalsN2), ncol = length(y)),
  N3 = matrix(rpois(length(inivalsN3), inivalsN3), ncol = length(y)),
  NP = matrix(rpois(length(inivalsNP), inivalsNP), ncol = length(y)), 
  NB = matrix(rpois(length(inivalsNB), inivalsNB), ncol = length(y)),
  NN = matrix(rpois(length(inivalsNN), inivalsNN), ncol = length(y)), 
  Ntot = rpois(length(blah),blah),
  SX = surv_traits_alt,
  SSSTA = surv_env3[,1:37],
  SPLUME = surv_env3[,39:75], 
  #SRAIN = surv_env3[,76:112], 
  STEMP = surv_env3[,113:149], 
  ISSTA = herit_env$av_s,
  IPLUME = herit_env$mean.2, 
  IX = herit_env$averagePC1*-1)} 

# run nimble code!
Rmodel <- nimbleModel(code = IPM2code_peng1,
                      constants = my.constants,
                      data = my.data,
                      inits = initial.values(), check=F)
Cmodel <- compileNimble(Rmodel,showCompilerOutput = T)
conf<- configureMCMC(Rmodel, useConjugacy = F, monitors=c( "l.p", "l.mphi", "l.mfec", "l.mim", "l.minh", "inh_var",
                                                           "sigma.y","alpha", "beta", "gamma", "delta", "epsilon", 
                                                           "Ntot", "N1", "NB", "y", "lambda", "geomean.lambda"),thin=10)

# model selection
conf$removeSamplers(c('alpha[4]','alpha[5]',
                      'beta[2,2]','beta[4,2]',
                      'beta[5,1]',
                      'beta[6,1]','beta[6,2]','beta[6,3]'))
#Cmodel$alpha[3:5] <- 0; Cmodel$beta[2,2] <- 0; Cmodel$beta[4,2] <- 0; Cmodel$beta[5,1] <- 0; Cmodel$beta[6,1] <- 0; Cmodel$beta[6,2] <- 0; Cmodel$beta[6,3] <- 0
#Rmodel$alpha[3:5] <- 0; Rmodel$beta[2,2] <- 0; Rmodel$beta[4,2] <- 0; Rmodel$beta[5,1] <- 0; Rmodel$beta[6,1] <- 0; Rmodel$beta[6,2] <- 0; Rmodel$beta[6,3] <- 0
Cmodel$alpha[4:5] <- 0; Cmodel$beta[2,2] <- 0; Cmodel$beta[4,2] <- 0; Cmodel$beta[5,1] <- 0; Cmodel$beta[6,1] <- 0; Cmodel$beta[6,2] <- 0; Cmodel$beta[6,3] <- 0
Rmodel$alpha[4:5] <- 0; Rmodel$beta[2,2] <- 0; Rmodel$beta[4,2] <- 0; Rmodel$beta[5,1] <- 0; Rmodel$beta[6,1] <- 0; Rmodel$beta[6,2] <- 0; Rmodel$beta[6,3] <- 0
#conf$addSampler(target = 'alpha[1:4]', type = "RW")

# run the rest
Rmcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel, resetFunctions = T)
chain_output <- runMCMC(Cmcmc, niter=10000, nburnin = 5000, nchains=3, progressBar=T)
