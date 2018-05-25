############################################################################################################
# Discrete generation, simple model of evolution of anadromy. 
# The evolvable trait here is the "underlying threshold for residency", i.e. the condition value above which 
# residency occurs, below which anadromy occurs.
# "Condition" here refers to a trait that varies across individuals randomly (variation is purely environmentally 
# driven in this simple model, i.e. there is no genetic basis to condition).  

############################################################################################################
rm(list=ls())        # clear all 

## PARAMETERS:

nyears <- 10        # number of years for each simulation run
nreps<- 1           # number of simulations to run

mu_thresh <- 0       # initial mean genetic value of the treshold trait
h2 <- 0.5            # heritability of the threshold trait 
Vp <- 1              # phenotypic variance of the threshold trait
Va <- Vp*h2          # additive genetic variance of the threshold trait
Ve <- Vp - Va        # environmental variance of the threshold trait (Vp = Va + Ve)

mu_cond  <- 0        # mean value of the condition trait
V_cond   <- 1        # variance of the condition trait (in this case, all the phenotypic variance is environmental)

N_females <- 500                     # Number of females 
N_males <- 500                       # Number of males 
N_init <- N_females + N_males       # initial population size (number of fry) is the sum of males and females 

# NB**** If any of the following 4 parameters are changed, this will drive evolution towards one or other tactic 
# NB**** (e.g. if S_anadromous is increased, then the pop will evolve towards a higher fraction of anadromous fish)

S_anadromousF <- 0.35 #  fraction of female anadromous fish that survive to breeding age
S_anadromousM <- 0.35 #  fraction of male anadromous fish that survive to breeding age
S_residentF <-   0.5  #  fraction of female resident fish that survive to breeding age
S_residentM <-   0.5  #  fraction of male resident fish that survive to breeding age
F_anadromous <- 10     #  fecundity of anadromous fish (number of offspring) 
F_resident <-   5      #  fecundity of resident fish (number of offspring)

K <-  100            # Carrying capacity (this is a simple ceiling population size, above which individuals 
# are randomly culled, to avoid the population sky-rocketing to huge sizes

# The following bits simply create empty vectors, in which the results of each model run can be stored: 
sim <- c()             # Store simulation number
year <- c()            # Store year number
pop.size1 <- c()       # Store population size of before survival and fecundity selection
pop.size1_F <- c()     # Store population size of females before survival and fecundity selection
pop.size1_M <- c()     # Store population size of males before survival and fecundity selection
pop.size2 <- c()       # Store population size of after survival and fecundity selection
pop.size2_F <- c()     # Store population size of females after survival and fecundity selection
pop.size2_M <- c()     # Store population size of males after survival and fecundity selection
mean.thresh <- c()   # Store realised mean genetic threshold across all individuals
va.realiz<-  c()     # Store realised variance in genetic thresholds across all individuals
frac.anadF <- c()     # Store realised fraction of anadromous female fish (number from 0 to 1)
frac.anadM <- c()     # Store realised fraction of anadromous male fish (number from 0 to 1)
real.surv.anadF <- c()# Store realised survival of anadromous female fish
real.surv.anadM <- c()# Store realised survival of anadromous male fish
real.surv.resF <- c() # Store realised survival of resident female fish
real.surv.resM <- c() # Store realised survival of resident male fish

###########################################################
### SIMULATION STARTS HERE!!! 
###########################################################

### cycle over simulations
for (Sim in 1:nreps) {
  
  ## Initiate the population with N_init number of breeding values for the threshold trait. These are the individual-specific,
  ## genetic values for the evolvable trait
  a_threshF <- rnorm(N_females,  mu_thresh, (sqrt(Va)))
  a_threshM <- rnorm(N_males,  mu_thresh, (sqrt(Va)))
  threshSize<-length(a_threshF) + length(a_threshM)
  
  ### cycle over years
  for (Y in 1:nyears)
  {
    
    ### additional random mortality if K exceeded, to keep population size in check
    if (length(a_threshF) + length(a_threshM)>K)
    {
      surv <- ifelse(runif(threshSize,0,1)<(K/threshSize),1,0) # sum(surv==1)
      # runif draws a random number from a uniform 
      #distribution, here between 0 and 1. If this number 
      #is less than the desired fraction of survivors (K/N), 
      #that individual survives, otherwise it dies
      a_threshF <- a_threshF[surv==1] # subset the breeding values for only those "individuals" that survive
      a_threshM <- a_threshM[surv==1] # subset the breeding values for only those "individuals" that survive
      a_threshF <- na.omit(a_threshF)
      a_threshM <- na.omit(a_threshM)
    }
    pop.size1 <- c(pop.size1,length(a_threshF) + length(a_threshM))  # calculate and store pop size ( = length of a_thresh vector)
    
    IDF <- 1:length(a_threshF)    # allocate arbitraty identities to each "individual"
    IDM <- 1:length(a_threshM) + 500   # allocate arbitraty identities to each "individual"
    
    e_threshF <- rnorm(length(a_threshF),  0, (sqrt(Ve))) # draw environmental deviations for the threshold trait
    e_threshM <- rnorm(length(a_threshM),  0, (sqrt(Ve))) # draw environmental deviations for the threshold trait
    
    z_threshF <- a_threshF + e_threshF     # the phenotypic value (z) for each individual is the sum of the genetic 
    # (breeding) value and the environmental deviation#
    z_threshM <- a_threshM + e_threshM
    
    condF <-  rnorm(length(a_threshF), mu_cond,  sqrt(V_cond)) # draw random condition values for each individual
    condM <-  rnorm(length(a_threshM), mu_cond,  sqrt(V_cond)) # draw random condition values for each individual
    
    anadromousF <- ifelse(condF > z_threshF, 0, 1)   # define migration tactics based on threshold model (if an individual's 
    # condition is > it's threshold value, it becomes resident (i.e. anadromous==0).
    # otherwise it becomes anadromous (anadromous == 1))
    
    anadromousM <- ifelse(condF > z_threshF, 0, 1)  
    
    
    frac.anadF <- c(frac.anadF, mean(anadromousF))   # calculate and store the fraction of anadromous fish in the pop at this timepoint
    frac.anadM <- c(frac.anadM, mean(anadromousM))   # calculate and store the fraction of anadromous fish in the pop at this timepoint
    
    anad_fishF <- IDF[anadromousF==1]   # pull out the IDs for the anadromous fish
    res_fishF  <- IDF[anadromousF==0]   # pull out the IDs for the resident fish
    anad_fishM <- IDM[anadromousM==1]   # pull out the IDs for the anadromous fish
    res_fishM  <- IDM[anadromousM==0]   # pull out the IDs for the resident fish
    
    # Calculate survival of anadromous fish by drawing a random number from uniform distribution between 0 and 1. If this number
    # is < S_anadromous (i.e. the inputted expected survival of anadromous fish), then that individual survives, otherwise it dies 
    surv_anadF <- runif(length(anad_fishF)) < rep(S_anadromousF, length(anad_fishF))
    surv_anadM <- runif(length(anad_fishM)) < rep(S_anadromousM, length(anad_fishM))
    
    # Do same as above for the resident fish:
    surv_resF <- runif(length(res_fishF)) < rep(S_residentF, length(res_fishF))
    surv_resM <- runif(length(res_fishM)) < rep(S_residentM, length(res_fishM))
    
    real.surv.anadF <- c(real.surv.anadF, mean(surv_anadF))  # Calculate and store the realised mean survival of anadromous fish, as a check
    real.surv.anadM <- c(real.surv.anadM, mean(surv_anadM))  # Calculate and store the realised mean survival of anadromous fish, as a check
    
    real.surv.resF <- c(real.surv.resF, mean(surv_resF))     # Calculate and store the realised mean survival of resident fish, as a check
    real.surv.resM <- c(real.surv.resM, mean(surv_resM))     # Calculate and store the realised mean survival of resident fish, as a check
     
    anad_fishF <- anad_fishF[surv_anadF]  # Extract the IDs of the surviving anadromous fish
    anad_fishM <- anad_fishM[surv_anadM]  # Extract the IDs of the surviving anadromous fish
    res_fishF <-  res_fishF[surv_resF]    # Extract the IDs of the resident anadromous fish
    res_fishM <-  res_fishM[surv_resM]    # Extract the IDs of the resident anadromous fish
    
    parents     <- c(anad_fishF,anad_fishM,res_fishF,res_fishM) # Make a new vector of parent IDs by stringing together the anadromous and resident IDs
    parents     <- data.frame(parents)
    parents$gender <- ifelse(parents > 500,"male", "female")
    parents     <- data.frame(parents)
    head(parents)
    # Make a vector of family sizes (number of offspring) for each parent, which will be used below in the mating and reproduction loop
    # This vector is ordered the same way as the vector of parents' IDs (anadromous fish first, then residents). 
    family_size <- c( rep(F_anadromous, (length(anad_fishF)) + length(anad_fishM) )
                                            , rep(F_resident, (length(res_fishF) + length(res_fishM))))
    # check with this  length(family_size) == length(parents$parents)
    # Mating and reproduction:
    #     sum(parents$parents<0) & sum(parents$parents>500)>=2                             # FALSE condition to test                         
    if (sum(parents$parents<500,na.rm = T) & sum(parents$parents>500,na.rm=T)>=2)                            # if more than 1 individual, mating happens
    {
      # sum(parents$parents<500) + sum(parents$parents>500)
      ### mating
      # npairs <- length(parents)%/%2                    # number of mating pairs = no. of parents ÷ 2, and rounded down
      motherPotential <- sample(parents$parents[parents$parents<500], replace=F)     # randomly select n=npairs "mothers" out of the list of parents IDs
      fatherPotential <- sample(parents$parents[parents$parents>500], replace=F)     # randomly select n=npairs "fathers" out of the list of parents IDs
     
    # we have to correct for the fact that the numbers of mothers and fathers won't match, this ifelse statement 
    # does that by sampling from the longer parent vector by the length of the smaller one 
    if(length(motherPotential) > length(fatherPotential)){
                         father<-fatherPotential;
                         mother<-sample(motherPotential,length(fatherPotential), replace = FALSE)
                          } else {
                        mother<-motherPotential;
                        father<-sample(fatherPotential,length(motherPotential), replace = FALSE)
                         }
      
      a_thresh_fath <- a_threshM[match(father,IDM)]      # extract the corresponding breeding values for these fathers from the vector a_thresh
      a_thresh_moth <- a_threshF[match(mother,IDF)]      # extract the corresponding breeding values for these mothers from the vector a_thresh
      mid <- (a_thresh_fath + a_thresh_moth)/2         # calculate the mid-parental genetic value (mean of breeding values of each parent)  
      
      ### breeding
      BV <- c()       # create an empty vector (BV= "breeding values") to store the new genetic values of the offspring
      
      for (i in 1:length(mother))  # cycle over the n mothers
      {
        # Generate new genetic values for the offpring that are centred on the mid-parental values of their parents,
        # plus a random deviation drawn from 0.5Va.  This essentially generates genetic variation among siblings,
        # with the expected genetic variation among siblings equal to half the population-wide additive genetic variance. 
        # This comes directly from quantitative genetic theory, see Chapter 4 of Roff 2009 book on "Modelling Evolution"
        BVo <- rep(mid[i] , family_size[parents==mother[i]]) + rnorm(family_size[parents==mother[i]], 0, sqrt(0.5*Va)) 
        # Add these offspring genetic values for each family to the vector BV:
        BV <- c(BV, BVo)
      }
      idx <- sample.int(length(BV),size=length(BV)/2,replace=FALSE) 
      a_threshF <- BV[idx]
      a_threshM <- BV[-idx]
      #  Now we replace the parental genetic values with a new list of offspring genetic values. 
      #  This is because we here assume that all parents die immediately after reproduction, along with their
      #  genetic values!
     
      
      
      ### store results
      sim <- c(sim,Sim)
      year <- c(year,Y)
      pop.size2 <- c(pop.size2,length(a_threshF) + length(a_threshM))   # store population size after survival and mating
      mean.thresh <- c(mean.thresh,mean(c(a_threshF,a_threshM, na.rm=T))) # store realised mean genetic value. 
      va.realiz<-  c(va.realiz,var(c(a_threshF,a_threshM, na.rm=T)))     # store realised variance in genetic values. 
      
    } else
    {
      sim <- c(sim,Sim)
      year <- c(year,Y)
      pop.size2 <- c(pop.size2,0)   # If the number of parents is <2, then simply store a 0 for pop size 
      mean.thresh <- c(mean.thresh,NA) # And store an NA here
      va.realiz<-  c(va.realiz, NA)    # And store an NA here
    }
  } # close the year loop
} # close the simulation replicate loop


# Create a data frame called r1 ("results 1") to store all the results:
r1 <- data.frame(sim,year,pop.size1,pop.size2,mean.thresh,va.realiz,frac.anadF,frac.anadM,real.surv.anadF,real.surv.anadM
                 ,real.surv.resF, real.surv.resM)

# Calculate the mean population size per year across all replicate simulations:
mN <- tapply(r1$pop.size1, r1$year, quantile, 0.5, na.rm=T)
lciN <- tapply(r1$pop.size1, r1$year, quantile, 0.05, na.rm=T) # calculate the lower confidence interval for this variable
uciN <- tapply(r1$pop.size1, r1$year, quantile, 0.95, na.rm=T) # calculate the upper confidence interval for this variable

# Calculate the mean fraction of anadromous fish per year across all replicate simulations:
mA <- tapply(r1$frac.anad, r1$year, quantile, 0.5, na.rm=T)
lcA <- tapply(r1$frac.anad, r1$year, quantile, 0.05, na.rm=T)
ucA <- tapply(r1$frac.anad, r1$year, quantile, 0.95, na.rm=T)

# Calculate the mean genetic threshold value per year across all replicate simulations:
mT <- tapply(r1$mean.thresh, r1$year, quantile, 0.5, na.rm=T)
lcT <- tapply(r1$mean.thresh, r1$year, quantile, 0.05, na.rm=T)
ucT <- tapply(r1$mean.thresh, r1$year, quantile, 0.95, na.rm=T)

# Calculate the mean realised survival of anadromous fish per year across all replicate simulations:
mSA <- tapply(r1$real.surv.anad, r1$year, quantile, 0.5, na.rm=T)
lcSA <- tapply(r1$real.surv.anad, r1$year, quantile, 0.05, na.rm=T)
ucSA <- tapply(r1$real.surv.anad, r1$year, quantile, 0.95, na.rm=T)

# Calculate the mean realised survival of resident fish per year across all replicate simulations:
mSR <- tapply(r1$real.surv.res, r1$year, quantile, 0.5, na.rm=T)
lcSR <- tapply(r1$real.surv.res, r1$year, quantile, 0.05, na.rm=T)
ucSR <- tapply(r1$real.surv.res, r1$year, quantile, 0.95, na.rm=T)

yr<- 1:nyears  # create a vector of year IDs, for plotting purposes below

par(mfrow=c(2,2))  # open a 2 x 2 tiled plotting window

# Plot 1 = Population size against year:
plot(yr, mN, xlab='Year', ylab='Population size', type='l', lwd=1, ylim=c(0,250), cex.lab=1.2)
points(yr, lciN, type='l', lwd=1, lty=2, ylim=c(0,250))
points(yr, uciN, type='l', lwd=1, lty=2, ylim=c(0,250))

# Plot 2 = Survival or anadromous and resident fish against year:
plot(yr, mSA, xlab='Year', ylab='Survival anadromous/resident', type='l', lwd=1, ylim=c(0,1), cex.lab=1.2, col="blue")
points(yr, mSR,  type='l', lwd=1, ylim=c(0,1), cex.lab=1.2, col="red")
legend("topright",legend=c("andadromous","resident"), col=c("blue","red"), lty=c(1,1))

# Plot 3 = Fraction of anadromous fish against year:
# ***NB This fraction will change as the population evolves!!
plot(yr, mA, xlab='Year', ylab='Fraction anadromous', type='l', lwd=1, ylim=c(0,1), cex.lab=1.2)
points(yr, lcA, type='l', lwd=1, lty=2, ylim=c(0,1))
points(yr, ucA, type='l', lwd=1, lty=2, ylim=c(0,1))

# Plot 4 = Mean genetic threshold against year:
# ***NB The mean genetic threshold will change as the population evolves!!
plot(yr, mT, xlab='Year', ylab='Mean genetic threshold', type='l', lwd=1, ylim=c(-5,5), cex.lab=1.2)
points(yr, lcT, type='l', lwd=1, lty=2, ylim=c(0,1))
points(yr, ucT, type='l', lwd=1, lty=2, ylim=c(0,1))



