##########################################################
#### ANALYSIS OF THE PUBLICATION ENTITLED #### 
#    Defining the in vivo clinical phenotype of artemisinin resistant falciparum malaria : 
#    A modelling approach 
##########################################################

setwd('d://FOLDERNAME/FILENAME') # set working directory

# THE PACKAGES USED IN THIS ANALYSIS
library(boot)
library(MASS)
library(segmented)
library(sm)
library(mixtools)
library(vioplot)


M<-5 #maximum number of components
pval<-0.1 #significance level for comparing models
nboot<-100 # number of iterations for bootstrap
nsim<-1000 # number of iterations for creating probaility of resistance vs HL graphs
P<-2 # use P or more samples to get geometric means and discard all other samples for permutation analysis
T<-100 # number of permutations for permutation analysis
smax=5000 # number of simulations for the day 3 postive model (must alter vector precres if smax is small)

nboot<-10 # number of iterations for bootstrap
nsim<-100 # number of iterations for creating probaility of resistance vs HL graphs
P<-2 # use P or more samples to get geometric means and discard all other samples for permutation analysis
T<-10 # number of permutations for permutation analysis
smax=100 # number of simulations for the day 3 postive model (must alter vector precres if smax is small)


#### ANALYSIS OF THE DATA FROM THE NORTHWESTERN THAI-MYANMAR BORDER STRATIFIED BY YEAR  ####

mixdat <-read.csv('simulated_cloneHLdata_SMRUbyyear.csv')  # import dataset
N<-ncol(mixdat)-1

# create output matrices
output.mu <- matrix(NA,nrow=M,ncol=N)
output.sigma <- matrix(NA,nrow=M,ncol=N)
output.lambda <- matrix(NA,nrow=M,ncol=N)
output.loglik <- matrix(NA,nrow=M,ncol=N)
output.mu.se <- matrix(NA,nrow=M,ncol=N)
output.sigma.se <- matrix(NA,nrow=M,ncol=N)
output.lambda.se <- matrix(NA,nrow=M,ncol=N)
AIC<-matrix(0,nrow=M,ncol=N)
AICdelta<-matrix(0,nrow=M,ncol=N)

nb<-na.omit(mixdat[,N+1])

# fit single component model
for (i in 1:N){
	# 1 COMPONENT LOG NORMAL
	nmixdat<-na.omit(mixdat[,i])
	lmixdat<- log(nmixdat)
	xll<-fitdistr(lmixdat,"normal")
	output.loglik[1,i]<- xll$loglik
	output.mu[1,i]<-xll$estimate[1]
	output.lambda[1,i]<-1
	output.sigma[1,i]<-xll$estimate[2]
	output.mu.se[1,i]<-xll$sd[1]
	output.sigma.se[1,i]<-xll$sd[2]
	output.lambda.se[1,i]<-0
	AIC[1,i]<-2*(3*1-1)-2*output.loglik[1,i]
	AICdelta[1,i]<-0
	}

# fit multiple component models sequntially
for (i in 1:N){
	nmixdat<-na.omit(mixdat[,i])
	lmixdat<- log(nmixdat)
	# >=2 COMPONENTS LOG NORMAL
	j<-1
  # stop if j-component model is more parsimonious than (j-1)-compnent model
	while((j<=M-1) && AICdelta[j,i]<=pval){
		j<-j+1
		res <- normalmixEM(lmixdat, lambda = matrix((1/j),nrow=1,ncol=j), mu = 2*(1:j)/j, sigma = 0.3*matrix(1,nrow=1,ncol=j))
		resboot <- boot.se(res, B = nboot)
		resboot[c("lambda.se", "mu.se", "sigma.se","loglik.se")]	
		output.loglik[j,i]<-res$loglik
		AIC[j,i]<-2*(3*j-1)-2*output.loglik[j,i]
    AICdelta[j,i]<-exp(-(AIC[j-1,i]-AIC[j,i])/2)
		if(AICdelta[j,i]<=pval){
		  output.mu[1:j,i]<-res$mu
      output.sigma[1:j,i]<-res$sigma
      output.lambda[1:j,i]<-res$lambda
      output.mu.se[1:j,i]<-resboot$mu.se
      output.sigma.se[1:j,i]<-resboot$sigma.se
      output.lambda.se[1:j,i]<-resboot$lambda.se			
		  }
		}
	}

##########################################################
# PLOTTING

# GRAPHS FOR SUPPORTING INFORMATION FILE 3
Sys.sleep(0.02)
for (ds in 1:N){
  Sys.sleep(0.02)
  nmixdat<-na.omit(mixdat[,ds])
  plam<-na.omit(output.lambda[,ds])
  pmu<-na.omit(output.mu[,ds])
  psig<-na.omit(output.sigma[,ds])
  hist(nmixdat,freq=FALSE,main = paste("Northwestern Thai-Myanmar border",2000+ds),xlab = "Clearance half-life (hours)",ylim=c(0,0.6),col="grey",lwd=2,ps=20,breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12))
  x <- seq(0.1, max(nmixdat), length=1000)

    hx<-plam[1]*dlnorm(x,meanlog=(pmu[1]),sdlog=psig[1])
    if(length(plam)>1){
      for(k in 2:length(plam)){
        hx<-hx+plam[k]*dlnorm(x,meanlog=(pmu[k]),sdlog=psig[k])
      }
    }
    lines(x,hx,col="red", lwd=5)
  Sys.sleep(0.02)
}

# TOP GRAPH OF FIGURE 2 
Sys.sleep(0.02)

vioplot(na.omit(mixdat[,1]),na.omit(mixdat[,2]), na.omit(mixdat[,3]), na.omit(mixdat[,4]), na.omit(mixdat[,5]), na.omit(mixdat[,6]), na.omit(mixdat[,7]),na.omit(mixdat[,8]), na.omit(mixdat[,9]), na.omit(mixdat[,10]), na.omit(mixdat[,11]), na.omit(mixdat[,12]),col="grey",pchMed=30)


mev1<-matrix(1,nrow=N*M,ncol=1)
mev2<-matrix(1,nrow=N*M,ncol=1)
mev3<-matrix(1,nrow=N*M,ncol=1)
for (a in 1:N){
  for (b in 1:M){
    mev1[M*(a-1)+b]<-2000+a+b/1000
    mev2[M*(a-1)+b]<-exp(output.mu[b,a])
    mev3[M*(a-1)+b]<-(output.lambda[b,a])^0.5
  }
}
dfx = data.frame(ev1=mev1, ev2=mev2, ev3=mev3)
symbols(x=dfx$ev1-2000, y=dfx$ev2, circles=dfx$ev3, inches=1/7, ann=F, bg="red", fg=NULL,xaxp = c(2001, 2015, 14),ylim=c(2,8),add=TRUE)
title(xlab = "Time (years)",ylab="Clearance half-life (hours)",ps=20)


# BOTTOM GRAPH OF FIGURE 2 (BLACK LINE) 

plot(seq(2001,2012),1-output.lambda[1,],type="l",ylim=c(0,1),xlab = "Time (years)",ylab="predicted proportion resistant")
par(new=T)
pres.low<-(1-output.lambda[1,])-output.lambda.se[1,]
pres.high<-(1-output.lambda[1,])+output.lambda.se[1,]
#polygon(c(seq(2001,2012), rev(seq(2001,2012))), c(pres.high, rev(pres.low)),col = "grey80", border = "black", lty="dashed")
polygon(c(seq(2001,2012), rev(seq(2001,2012))), c(pres.high, rev(pres.low)),col = "grey80", border = NA)
par(new=T)
plot(seq(2001,2012),1-output.lambda[1,],type="l",ylim=c(0,1),xlab = "Time (years)",ylab="predicted proportion resistant",lwd=5,ps=20)
pres.med<-1-output.lambda[1,]
##########################################################

#### ANALYSIS OF THE COMBINED DATASET STRATIFIED BY COUNTRY  ####

mixdat <-read.csv('simulated_cloneHLdata_ALLbycountry.csv')  # only samples with clones from SMRU

N<-ncol(mixdat)-1


output.mu <- matrix(NA,nrow=M,ncol=N)
output.sigma <- matrix(NA,nrow=M,ncol=N)
output.lambda <- matrix(NA,nrow=M,ncol=N)
output.loglik <- matrix(NA,nrow=M,ncol=N)
output.mu.se <- matrix(NA,nrow=M,ncol=N)
output.sigma.se <- matrix(NA,nrow=M,ncol=N)
output.lambda.se <- matrix(NA,nrow=M,ncol=N)
AIC<-matrix(0,nrow=M,ncol=N)
AICdelta<-matrix(0,nrow=M,ncol=N)

nb<-na.omit(mixdat[,N+1])

for (i in 1:N){
  # 1 COMPONENT LOG NORMAL
  nmixdat<-na.omit(mixdat[,i])
  lmixdat<- log(nmixdat)
  xll<-fitdistr(lmixdat,"normal")
  output.loglik[1,i]<- xll$loglik
  output.mu[1,i]<-xll$estimate[1]
  output.lambda[1,i]<-1
  output.sigma[1,i]<-xll$estimate[2]
  output.mu.se[1,i]<-xll$sd[1]
  output.sigma.se[1,i]<-xll$sd[2]
  output.lambda.se[1,i]<-0
  AIC[1,i]<-2*(3*1-1)-2*output.loglik[1,i]
  AICdelta[1,i]<-0
}
for (i in 1:N){
  nmixdat<-na.omit(mixdat[,i])
  lmixdat<- log(nmixdat)
  # >=2 COMPONENTS LOG NORMAL
  j<-1
  while((j<=M-1) && AICdelta[j,i]<=pval){
    j<-j+1
    res <- normalmixEM(lmixdat, lambda = matrix((1/j),nrow=1,ncol=j), mu = 2*(1:j)/j, sigma = 0.3*matrix(1,nrow=1,ncol=j))
    resboot <- boot.se(res, B = nboot)
    resboot[c("lambda.se", "mu.se", "sigma.se","loglik.se")]	
    output.loglik[j,i]<-res$loglik
    AIC[j,i]<-2*(3*j-1)-2*output.loglik[j,i]
    AICdelta[j,i]<-exp(-(AIC[j-1,i]-AIC[j,i])/2)
    if(AICdelta[j,i]<=pval){
      output.mu[1:j,i]<-res$mu
      output.sigma[1:j,i]<-res$sigma
      output.lambda[1:j,i]<-res$lambda
      output.mu.se[1:j,i]<-resboot$mu.se
      output.sigma.se[1:j,i]<-resboot$sigma.se
      output.lambda.se[1:j,i]<-resboot$lambda.se			
    }
  }
}


##########################################################
# PLOTTING

# A PLOT OF THE DATA VERSUS MODEL PREDICTIONS (NOT INCLUDED IN PUBLICATION) 
Sys.sleep(0.02)
for (ds in 1:N){
  Sys.sleep(0.02)
  nmixdat<-na.omit(mixdat[,ds])
  plam<-na.omit(output.lambda[,ds])
  pmu<-na.omit(output.mu[,ds])
  psig<-na.omit(output.sigma[,ds])
  hist(nmixdat,freq=FALSE,main = paste("dataset",ds,",",length(plam),"components"),xlab = "Clearance half-life (hours)",ylim=c(0,0.6),col="grey",lwd=2,ps=20,breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12))
  x <- seq(0.1, max(nmixdat), length=1000)
  
  hx<-plam[1]*dlnorm(x,meanlog=(pmu[1]),sdlog=psig[1])
  if(length(plam)>1){
    for(k in 2:length(plam)){
      hx<-hx+plam[k]*dlnorm(x,meanlog=(pmu[k]),sdlog=psig[k])
    }
  }
  lines(x,hx,col="red", lwd=5)
  Sys.sleep(0.02)
}

# GRAPH OF FIGURE 3
Sys.sleep(0.02)
mev1<-matrix(1,nrow=N*M,ncol=1)
mev2<-matrix(1,nrow=N*M,ncol=1)
mev3<-matrix(1,nrow=N*M,ncol=1)
for (a in 1:N){
  for (b in 1:M){
    mev1[M*(a-1)+b]<-2000+a+b/1000
    mev2[M*(a-1)+b]<-exp(output.mu[b,a])
    mev3[M*(a-1)+b]<-(output.lambda[b,a])^0.5
  }
}
dfx = data.frame(ev1=mev1, ev2=mev2, ev3=mev3)
vioplot(na.omit(mixdat[,1]),na.omit(mixdat[,2]), na.omit(mixdat[,3]),col="grey",pchMed=30)
symbols(x=dfx$ev1-2000, y=dfx$ev2, circles=dfx$ev3, inches=1/5, ann=F, bg="red", fg=NULL,xaxp = c(2001, 2015, 14),ylim=c(2,8),add=TRUE,ps=20)
title(xlab = "Site",ylab="Clearance half-life (hours)",ps=20)

sensmu<-output.mu[1,3]
sensmu.se<-output.mu.se[1,3]
senssigma<-output.sigma[1,3]
senssigma.se<-output.sigma.se[1,3]
resmu<-output.mu[2,3]
resmu.se<-output.mu.se[2,3]
ressigma<-output.sigma[2,3]
ressigma.se<-output.sigma.se[2,3]


#### APPLICAITON OF THE MODEL TO PREDICT THE LIKELIHOOD OF AN INDIVIDUAL INFECTION BEING RESISTANT ####

hlx<-seq(0,10,by=0.1)

# ASSUME AN UNDERLYING PROPORTION RESISISTANT OF 0.1
presistant<-0.1

sensmean<-rnorm(nsim,mean=(output.mu[1,3]),sd=(output.mu.se[1,3]))
resmean<-rnorm(nsim,mean=(output.mu[2,3]),sd=(output.mu.se[2,3]))
senssd<-rnorm(nsim,mean=(output.sigma[1,3]),sd=(output.sigma.se[1,3]))
ressd<-rnorm(nsim,mean=(output.sigma[2,3]),sd=(output.sigma.se[2,3]))

y<-matrix(NA,nrow=nsim,ncol=length(hlx))

for (i in 1:nsim){
  y[i,]<-presistant*dlnorm(hlx,meanlog=resmean[i],sdlog=ressd[i])/((1-presistant)*dlnorm(hlx,meanlog=sensmean[i],sdlog=senssd[i])+presistant*dlnorm(hlx,meanlog=resmean[i],sdlog=ressd[i]))    
}

propres.2.5<-matrix(NA,nrow=1,ncol=length(hlx))
propres.97.5<-matrix(NA,nrow=1,ncol=length(hlx))
propres.5<-matrix(NA,nrow=1,ncol=length(hlx))
propres.95<-matrix(NA,nrow=1,ncol=length(hlx))
propres.10<-matrix(NA,nrow=1,ncol=length(hlx))
propres.90<-matrix(NA,nrow=1,ncol=length(hlx))
propres.25<-matrix(NA,nrow=1,ncol=length(hlx))
propres.75<-matrix(NA,nrow=1,ncol=length(hlx))
propres.50<-matrix(NA,nrow=1,ncol=length(hlx))

for (i in 1:length(hlx)){
  propres.2.5[i]<-quantile(y[,i], probs = 0.025, na.rm = TRUE)
  propres.97.5[i]<-quantile(y[,i], probs = 0.975, na.rm = TRUE)  
  propres.5[i]<-quantile(y[,i], probs = 0.05, na.rm = TRUE)
  propres.95[i]<-quantile(y[,i], probs = 0.95, na.rm = TRUE)  
  propres.10[i]<-quantile(y[,i], probs = 0.1, na.rm = TRUE)
  propres.90[i]<-quantile(y[,i], probs = 0.9, na.rm = TRUE)  
  propres.25[i]<-quantile(y[,i], probs = 0.25, na.rm = TRUE)
  propres.75[i]<-quantile(y[,i], probs = 0.75, na.rm = TRUE)  
  propres.50[i]<-quantile(y[,i], probs = 0.50, na.rm = TRUE)  
}

# PLOT GREEN LINE IN FIGURE 4
plot(hlx,propres.50,type="l",ylim=c(0,1),xlab = "Clearance half-life (hours)",ylab="predicted proportion resistant",ps=20,col="darkgreen")
par(new=T)
polygon(c(hlx, rev(hlx)), c(propres.97.5, rev(propres.2.5)),col = "aquamarine1", border = NA)
par(new=T)
polygon(c(hlx, rev(hlx)), c(propres.95, rev(propres.5)),col = "aquamarine2", border = NA)
par(new=T)
polygon(c(hlx, rev(hlx)), c(propres.90, rev(propres.10)),col = "aquamarine3", border = NA)
par(new=T)
polygon(c(hlx, rev(hlx)), c(propres.75, rev(propres.25)),col = "aquamarine4", border = NA)
par(new=T)
plot(hlx,propres.50,type="l", lwd=2,ylim=c(0,1),xlab = "Clearance half-life (hours)",ylab="predicted proportion resistant",ps=20,col="darkgreen")

# ASSUME AN UNDERLYING PROPORTION RESISISTANT OF 0.5
presistant<-0.5

sensmean<-rnorm(nsim,mean=(output.mu[1,4]),sd=(output.mu.se[1,4]))
resmean<-rnorm(nsim,mean=(output.mu[2,4]),sd=(output.mu.se[2,4]))
senssd<-rnorm(nsim,mean=(output.sigma[1,4]),sd=(output.sigma.se[1,4]))
ressd<-rnorm(nsim,mean=(output.sigma[2,4]),sd=(output.sigma.se[2,4]))

y<-matrix(NA,nrow=nsim,ncol=length(hlx))

for (i in 1:nsim){
  y[i,]<-presistant*dlnorm(hlx,meanlog=resmean[i],sdlog=ressd[i])/((1-presistant)*dlnorm(hlx,meanlog=sensmean[i],sdlog=senssd[i])+presistant*dlnorm(hlx,meanlog=resmean[i],sdlog=ressd[i]))    
}

propres.2.5<-matrix(NA,nrow=1,ncol=length(hlx))
propres.97.5<-matrix(NA,nrow=1,ncol=length(hlx))
propres.5<-matrix(NA,nrow=1,ncol=length(hlx))
propres.95<-matrix(NA,nrow=1,ncol=length(hlx))
propres.10<-matrix(NA,nrow=1,ncol=length(hlx))
propres.90<-matrix(NA,nrow=1,ncol=length(hlx))
propres.25<-matrix(NA,nrow=1,ncol=length(hlx))
propres.75<-matrix(NA,nrow=1,ncol=length(hlx))
propres.50<-matrix(NA,nrow=1,ncol=length(hlx))

for (i in 1:length(hlx)){
  propres.2.5[i]<-quantile(y[,i], probs = 0.025, na.rm = TRUE)
  propres.97.5[i]<-quantile(y[,i], probs = 0.975, na.rm = TRUE)  
  propres.5[i]<-quantile(y[,i], probs = 0.05, na.rm = TRUE)
  propres.95[i]<-quantile(y[,i], probs = 0.95, na.rm = TRUE)  
  propres.10[i]<-quantile(y[,i], probs = 0.1, na.rm = TRUE)
  propres.90[i]<-quantile(y[,i], probs = 0.9, na.rm = TRUE)  
  propres.25[i]<-quantile(y[,i], probs = 0.25, na.rm = TRUE)
  propres.75[i]<-quantile(y[,i], probs = 0.75, na.rm = TRUE)  
  propres.50[i]<-quantile(y[,i], probs = 0.50, na.rm = TRUE)  
}

# PLOT BLUE LINE IN FIGURE 4
par(new=T)
plot(hlx,propres.50,type="l",ylim=c(0,1),xlab = "Clearance half-life (hours)",ylab="predicted proportion resistant",ps=20,col="blue4")
par(new=T)
polygon(c(hlx, rev(hlx)), c(propres.97.5, rev(propres.2.5)),col = "steelblue1", border = NA)
par(new=T)
polygon(c(hlx, rev(hlx)), c(propres.95, rev(propres.5)),col = "dodgerblue2", border = NA)
par(new=T)
polygon(c(hlx, rev(hlx)), c(propres.90, rev(propres.10)),col = "dodgerblue3", border = NA)
par(new=T)
polygon(c(hlx, rev(hlx)), c(propres.75, rev(propres.25)),col = "dodgerblue4", border = NA)
par(new=T)
plot(hlx,propres.50,type="l", lwd=2,ylim=c(0,1),xlab = "Clearance half-life (hours)",ylab="predicted proportion resistant",ps=20,col="blue4")

# ASSUME AN UNDERLYING PROPORTION RESISISTANT OF 0.5
presistant<-0.9

sensmean<-rnorm(nsim,mean=(output.mu[1,4]),sd=(output.mu.se[1,4]))
resmean<-rnorm(nsim,mean=(output.mu[2,4]),sd=(output.mu.se[2,4]))
senssd<-rnorm(nsim,mean=(output.sigma[1,4]),sd=(output.sigma.se[1,4]))
ressd<-rnorm(nsim,mean=(output.sigma[2,4]),sd=(output.sigma.se[2,4]))

y<-matrix(NA,nrow=nsim,ncol=length(hlx))

for (i in 1:nsim){
  y[i,]<-presistant*dlnorm(hlx,meanlog=resmean[i],sdlog=ressd[i])/((1-presistant)*dlnorm(hlx,meanlog=sensmean[i],sdlog=senssd[i])+presistant*dlnorm(hlx,meanlog=resmean[i],sdlog=ressd[i]))    
}

propres.2.5<-matrix(NA,nrow=1,ncol=length(hlx))
propres.97.5<-matrix(NA,nrow=1,ncol=length(hlx))
propres.5<-matrix(NA,nrow=1,ncol=length(hlx))
propres.95<-matrix(NA,nrow=1,ncol=length(hlx))
propres.10<-matrix(NA,nrow=1,ncol=length(hlx))
propres.90<-matrix(NA,nrow=1,ncol=length(hlx))
propres.25<-matrix(NA,nrow=1,ncol=length(hlx))
propres.75<-matrix(NA,nrow=1,ncol=length(hlx))
propres.50<-matrix(NA,nrow=1,ncol=length(hlx))

for (i in 1:length(hlx)){
  propres.2.5[i]<-quantile(y[,i], probs = 0.025, na.rm = TRUE)
  propres.97.5[i]<-quantile(y[,i], probs = 0.975, na.rm = TRUE)  
  propres.5[i]<-quantile(y[,i], probs = 0.05, na.rm = TRUE)
  propres.95[i]<-quantile(y[,i], probs = 0.95, na.rm = TRUE)  
  propres.10[i]<-quantile(y[,i], probs = 0.1, na.rm = TRUE)
  propres.90[i]<-quantile(y[,i], probs = 0.9, na.rm = TRUE)  
  propres.25[i]<-quantile(y[,i], probs = 0.25, na.rm = TRUE)
  propres.75[i]<-quantile(y[,i], probs = 0.75, na.rm = TRUE)  
  propres.50[i]<-quantile(y[,i], probs = 0.50, na.rm = TRUE)  
}

# PLOT PURPLE LINE IN FIGURE 4
par(new=T)
plot(hlx,propres.50,type="l",ylim=c(0,1),xlab = "Clearance half-life (hours)",ylab="predicted proportion resistant",ps=20,col="purple4")
par(new=T)
polygon(c(hlx, rev(hlx)), c(propres.97.5, rev(propres.2.5)),col = "plum3", border = NA)
par(new=T)
polygon(c(hlx, rev(hlx)), c(propres.95, rev(propres.5)),col = "mediumpurple2", border = NA)
par(new=T)
polygon(c(hlx, rev(hlx)), c(propres.90, rev(propres.10)),col = "mediumpurple3", border = NA)
par(new=T)
polygon(c(hlx, rev(hlx)), c(propres.75, rev(propres.25)),col = "mediumpurple4", border = NA)
par(new=T)
plot(hlx,propres.50,type="l", lwd=2,ylim=c(0,1),xlab = "Clearance half-life (hours)",ylab="predicted proportion resistant",ps=20,col="purple4")


##########################################################
#### PERMUTATION ANALYSIS ####
##########################################################
# FOR RESULTS SECTION "COMPARISON WITH DISTRIBUTIONS BY GENETIC RELATEDNESS"

# create a new mixdat file of all geometric means of true and false clone labels
clone <-read.csv('simulated_Allclones.csv')

clonemixdat<-function(data){
  C<-data
  C<-C[order(C[,2]),]
  H<-matrix(NA,nrow=length(C[,1]),ncol=5)
  H[,1]<-C[,1]
  H[,2]<-C[,2]
  L<-length(C[,1])
  
  # number of samples
  H[1,3]<-1
  for (k in 2:L){
    if(H[k,2]==H[k-1,2]){
      H[k,3]<-H[k-1,3]+1
    }
    if (H[k,2]!=H[k-1,2]){
      H[k,3]<-1
    }
  }
  # total number of samples
  H[L,4]<-H[L,3]
  for (k in 1:(L-1)){
    if (H[L-k,3]<H[L-k+1,3]){
      H[L-k,4]<-H[L-k+1,4]
    }
    if (H[L-k,3]>=H[L-k+1,3]){
      H[L-k,4]<-H[L-k,3]
    }
  }
  #remove unwanted samples
  H2<-H[H[,4]>=P,]
  L2<-length(H2[,1])  
  # get geomeans
  for (k in 1:L2){
    if (H2[k,3]==H2[k,4]){
      H2[k,5]<-exp((1/H2[k,4])*sum(log(H2[(k-H2[k,4]+1):k])))
    }
  }
  H2
}

H2<-clonemixdat(clone)
L2<-length(H2[,1])
MIX<-matrix(NA,nrow=L2,ncol=T+2)

MIX[,1]<-H2[,5]

# shuffle clone numbers
for (t in 2:(T+1)){
  H2[,2]<-H2[sample(L2),2] 
  H2[,3:5]<-NA
  H2<-clonemixdat(H2)
  MIX[,t]<-H2[,5]
}

for (i in 1:T+1){
  MIX[i+1,T+2]=length(MIX[,i])-sum(is.na(MIX[,i]))  
}

##########################################################
mixdat <-MIX  # random permutations of clone lables

N<-ncol(mixdat)-1


output.mu <- matrix(NA,nrow=M,ncol=N)
output.sigma <- matrix(NA,nrow=M,ncol=N)
output.lambda <- matrix(NA,nrow=M,ncol=N)
output.loglik <- matrix(NA,nrow=M,ncol=N)
output.mu.se <- matrix(NA,nrow=M,ncol=N)
output.sigma.se <- matrix(NA,nrow=M,ncol=N)
output.lambda.se <- matrix(NA,nrow=M,ncol=N)
AIC<-matrix(0,nrow=M,ncol=N)
AICdelta<-matrix(0,nrow=M,ncol=N)

nb<-na.omit(mixdat[,N+1])

for (i in 1:N){
  # 1 COMPONENT LOG NORMAL
  nmixdat<-na.omit(mixdat[,i])
  lmixdat<- log(nmixdat)
  xll<-fitdistr(lmixdat,"normal")
  output.loglik[1,i]<- xll$loglik
  output.mu[1,i]<-xll$estimate[1]
  output.lambda[1,i]<-1
  output.sigma[1,i]<-xll$estimate[2]
  output.mu.se[1,i]<-xll$sd[1]
  output.sigma.se[1,i]<-xll$sd[2]
  output.lambda.se[1,i]<-0
  AIC[1,i]<-2*(3*1-1)-2*output.loglik[1,i]
  AICdelta[1,i]<-0
}
for (i in 1:N){
  nmixdat<-na.omit(mixdat[,i])
  lmixdat<- log(nmixdat)
  # >=2 COMPONENTS LOG NORMAL
  j<-1
  while((j<=M-1) && AICdelta[j,i]<=pval){
    j<-j+1
    res <- normalmixEM(lmixdat, lambda = matrix((1/j),nrow=1,ncol=j), mu = 2*(1:j)/j, sigma = 0.3*matrix(1,nrow=1,ncol=j))
    resboot <- boot.se(res, B = nboot)
    resboot[c("lambda.se", "mu.se", "sigma.se","loglik.se")]	
    output.loglik[j,i]<-res$loglik
    AIC[j,i]<-2*(3*j-1)-2*output.loglik[j,i]
    AICdelta[j,i]<-exp(-(AIC[j-1,i]-AIC[j,i])/2)
    if(AICdelta[j,i]<=pval){
      output.mu[1:j,i]<-res$mu
      output.sigma[1:j,i]<-res$sigma
      output.lambda[1:j,i]<-res$lambda
      output.mu.se[1:j,i]<-resboot$mu.se
      output.sigma.se[1:j,i]<-resboot$sigma.se
      output.lambda.se[1:j,i]<-resboot$lambda.se			
    }
  }
}


##########################################################
# PLOTTING
ncomponent<-matrix(NA,nrow=1,ncol=T+1)

# A PLOT OF THE RESULTS OF THE PERMUATION ANALYSIS (NOT INCLUDED IN PUBLICATION)
for (i in 1: length(output.lambda[1,])){
  ncomponent[i]<-length(na.omit(output.lambda[,i]))
}
plot(seq(1,(length(MIX[1,])-1)),ncomponent,type="l",xlab = "simulation number",ylab="predictedmean half-life (hours) and proportion resistant",lwd=1,ps=20,ylim=c(0,10))
par(new=T)
plot(seq(1,(length(MIX[1,])-1)),exp(output.mu[1,]),type="p", pch=16,xlab = "simulation number",ylab="predictedmean half-life (hours) and proportion resistant",lwd=1,ps=20, col="green",ylim=c(0,10))
par(new=T)
plot(seq(1,(length(MIX[1,])-1)),exp(output.mu[2,]),type="p", pch=16,xlab = "simulation number",ylab="predictedmean half-life (hours) and proportion resistant",lwd=1,ps=20, col="red",ylim=c(0,10))


percbimodal=100*((length(output.lambda[1,])-sum(output.lambda[1,]==1)-1)/(length(output.lambda[1,])-1))

# OUTPUT THE PERCENTAGE OF ALL SHUFFLED DATASETS WHICH PREDICT A BIMODAL DISTRIBUTION
percbimodal


##########################################################

##### DAY 3 POSITIVE MODEL SIMULATION #####
# THE MODEL IS USED TO SIMULATE THE PROPORTIONS OF PATIENTS STILL PARASITAEMIC ON DAYS 1, 2 AND 3 OF TREATMENT
# A RANGE OF UNDERLYING PROPORTIONS OF RESISTANT INFECTIONS IS ASSUMED
# A RANGE OF GEOMETRIC MEAN ADMISSION PARASITAEMIA IS ASSUMED
# A RANGE OF SAMPLE SIZES IS ASSUMED

ipdat = read.csv('simulated3Dposdata.csv', header = TRUE)


for (nn in 1:dim(ipdat)[1]){
  
  rm(list=c("hres","hsen","lag","TC","nTC1","nTC2","nTC3","percres","perc1day","perc2day","perc3day"))
  
  arrayN=ipdat[nn,2]
  arraymuip=log10(exp(ipdat[nn,3]))
  sigip=log10(exp(ipdat[nn,4]))
  #percres=c(1,2,3,4,5,10,20,40,60,80,99)
  percres=c(2,4,10,20,50,99)  
  #percres=c(1, 2, 4, 6, 8, 10, 15, 20, 30, 40, 60, 80, 99)
  
  cil=2.5
  ciu=97.5
  
  for (px in 1:length(arrayN)){
    for (py in 1:length(arraymuip)){
      N=arrayN[px] # number of patients in simulation
      muip=arraymuip[py] # mean log initial para
      perc1day <- perc2day<- perc3day <- array(NA, c(length(percres), smax))
      p101<-p102<-p103<-array(NA, nn)
      for (s in 1:smax){
        100*(((py+(px-1)*length(arraymuip))-1)*smax+s)/(length(arrayN)*length(arraymuip)*smax)
        for (rr in 1:length(percres)){
          pres=percres[rr] # precentage resistant [y axis]
          nres=round(pres*N/100)
          nsen=N-nres
          dlim=50 # detection limit  
          # generate a ditribution of initial parasitaemias
          sigip=sigip   #max(0,random('norm',0.75,0.2)) # sd log initial para
          ip = rnorm(N, mean = muip*rep(1,N), sd = sigip)
          
          # generate distribution of sensitive halflives
          muhs_log = max(c(0,rnorm(1, mean = sensmu, sd = sensmu.se))) # mean sensitive halflife in hours
          sighs_log=max(c(0,rnorm(1, mean = senssigma, sd = senssigma.se))) # sd sensitive halflife in hours
          hsen=exp(rnorm(nsen, mean = muhs_log*rep(1,nsen), sd = sighs_log)) # log-normal distribution
          
          # generate distribution of resistant halflives
          muhr_log= max(c(0.01,rnorm(1, mean = resmu, sd = resmu.se))) # mean resistant halflife in hours
          sighr_log=max(c(0.01,rnorm(1, mean = ressigma, sd = ressigma.se))) # sd resistant halflife in hours
          hres=exp(rnorm(nres, mean = muhr_log*rep(1,nres), sd = sighr_log)) # log-normal distribution
          
          # lag times
          mulag=5
          siglag=3
          lag = apply(cbind(rep(0,N), rlnorm(N, meanlog = log(mulag)*rep(1,N), sdlog = siglag/mulag)), 1, max)
          
          # calculate times to clear
          TC = array(NA, nsen+nres)
          for (i in 1:nsen){
            TC[i]=(hsen[i]*(log(10^ip[i])-log(dlim))/(log(2)))+lag[i]
          }
          for (i in 1:nres){
            TC[nsen+i]=(hres[i]*(log(10^ip[nsen+i])-log(dlim))/(log(2)))+lag[i]
          }
          
          # DAY 1	
          nTC1 = array(NA, N)
          tob1 = rlnorm(1, meanlog = log(24), sdlog = 3.5/24)
          for (i in 1:N){
            if (TC[i]>=rnorm(1, mean = tob1, sd = 6)){ 
              nTC1[i]=1
            }else{	
              nTC1[i]=0	
            }
          }
          perc1day[rr,s]=sum(nTC1)
          
          # DAY 2
          nTC2 = array(NA, N)
          tob2=rlnorm(1, meanlog = log(48), sdlog = 3.5/48)
          for (i in 1:N){
            if (TC[i]>=rnorm(1, mean = tob2, sd = 6)){
              nTC2[i]=1
            }else{
              nTC2[i]=0
            }
          }
          perc2day[rr,s]=sum(nTC2) 
          
          # DAY 3
          nTC3 = array(NA, N)
          tob3=rlnorm(1, meanlog = log(72), sdlog = 3.5/72)
          for (i in 1:N){
            if (TC[i]>=rnorm(1, mean = tob3, sd = 6)){
              nTC3[i]=1
            }else{
              nTC3[i]=0
            }
          }
          perc3day[rr,s]=sum(nTC3)
        }
      }
      # GRAPHS FOR SUPPORTING INFORMATION FILE 3
      #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
      plot(100*apply(perc1day, 1, quantile, prob = cil/10)/N,percres,col = "black", lty = "dashed", lwd = 5, typ = "l", xlim = c(0,100),ylim = c(0,100),main = sprintf("N=%s, P0 = %s, Day 1", ipdat[nn,2], round(exp(ipdat[nn,3]))), xlab = "", ylab = "")
      lines(100*apply(perc1day, 1, quantile, prob = 0.5)/N,percres,col = "black", lty = "solid", lwd = 5, typ = "l")
      lines(100*apply(perc1day, 1, quantile, prob = ciu/100)/N,percres,col = "black", lty = "dashed", lwd = 5, typ = "l")
      lines(c(ipdat[nn,5],ipdat[nn,5]),c(0,100),col = "red", lty = "dashed", lwd = 5, typ = "l")
      #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
      plot(100*apply(perc2day, 1, quantile, prob = cil/100)/N,percres,col = "black", lty = "dashed", lwd = 5, typ = "l", xlim = c(0,100),ylim = c(0,100),main = sprintf("N=%s, P0 = %s, Day 2", ipdat[nn,2], round(exp(ipdat[nn,3]))),xlab = "", ylab = "")
      lines(100*apply(perc2day, 1, quantile, prob = 0.5)/N,percres,col = "black", lty = "solid", lwd = 5, typ = "l")
      lines(100*apply(perc2day, 1, quantile, prob = ciu/100)/N,percres,col = "black", lty = "dashed", lwd = 5, typ = "l")
      lines(c(ipdat[nn,6],ipdat[nn,6]),c(0,100),col = "red", lty = "dashed", lwd = 5, typ = "l")
      #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
      # FOR FIGURE 5 
      plot(100*apply(perc3day, 1, quantile, prob = cil/100)/N,percres,col = "black", lty = "dashed", lwd = 5, typ = "l", xlim = c(0,100),ylim = c(0,100),main = sprintf("N=%s, P0 = %s, Day 3", ipdat[nn,2], round(exp(ipdat[nn,3]))), xlab = "", ylab = "")
      lines(100*apply(perc3day, 1, quantile, prob = 0.5)/N,percres,col = "black", lty = "solid", lwd = 5, typ = "l")
      lines(100*apply(perc3day, 1, quantile, prob = ciu/100)/N,percres,col = "black", lty = "dashed", lwd = 5, typ = "l")
      lines(c(ipdat[nn,7],ipdat[nn,7]),c(0,100),col = "red", lty = "dashed", lwd = 5, typ = "l")
      #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
    }
  }
}


##########################################################

##### DAY 3 POSITIVE MODEL CASE STUDY #####
# A SIMULATION MODEL IS USED TO PROVIDE AN ALTERNATIVE ESTIMATE OF THE PROPORTION RESISTANT ON THE THAI-MYANMAR BORDER
# STARTING WITH DAY-3 POSITIVITY DATA THE SIMULATION MODEL IS USED TO INFER THE PROPORTION RESISTANT EACH YEAR

ipdat = read.csv('simulated_SMRU3Dposdata.csv', header = TRUE)

pres3D<-matrix(NA,ncol=3,nrow=length(ipdat[,1]))


for (nn in 1:dim(ipdat)[1]){
  
  rm(list=c("hres","hsen","lag","TC","nTC1","nTC2","nTC3","percres","perc1day","perc2day","perc3day"))
  
  arrayN=ipdat[nn,2]
  arraymuip=log10(exp(ipdat[nn,3]))
  sigip=log10(exp(ipdat[nn,4]))
  #percres=c(1,2,3,4,5,10,20,40,60,80,99)
  percres=c(2,4,10,20,50,99)
  #percres=c(1, 2, 4, 6, 8, 10, 15, 20, 30, 40, 60, 80, 99)
  
  cil=2.5
  ciu=97.5
  
  for (px in 1:length(arrayN)){
    for (py in 1:length(arraymuip)){
      N=arrayN[px] # number of patients in simulation
      muip=arraymuip[py] # mean log initial para
      perc1day <- perc2day<- perc3day <- array(NA, c(length(percres), smax))
      p101<-p102<-p103<-array(NA, nn)
      for (s in 1:smax){
        100*(((py+(px-1)*length(arraymuip))-1)*smax+s)/(length(arrayN)*length(arraymuip)*smax)
        for (rr in 1:length(percres)){
          pres=percres[rr] # precentage resistant [y axis]
          nres=round(pres*N/100)
          nsen=N-nres
          dlim=50 # detection limit  
          # generate a ditribution of initial parasitaemias
          sigip=sigip   #max(0,random('norm',0.75,0.2)) # sd log initial para
          ip = rnorm(N, mean = muip*rep(1,N), sd = sigip)
          
          # generate distribution of sensitive halflives
          muhs_log = max(c(0,rnorm(1, mean = sensmu, sd = sensmu.se))) # mean sensitive halflife in hours
          sighs_log=max(c(0,rnorm(1, mean = senssigma, sd = senssigma.se))) # sd sensitive halflife in hours
          hsen=exp(rnorm(nsen, mean = muhs_log*rep(1,nsen), sd = sighs_log)) # log-normal distribution
          
          # generate distribution of resistant halflives
          muhr_log= max(c(0.01,rnorm(1, mean = resmu, sd = resmu.se))) # mean resistant halflife in hours
          sighr_log=max(c(0.01,rnorm(1, mean = ressigma, sd = ressigma.se))) # sd resistant halflife in hours
          hres=exp(rnorm(nres, mean = muhr_log*rep(1,nres), sd = sighr_log)) # log-normal distribution
          
          # lag times
          mulag=5
          siglag=3
          lag = apply(cbind(rep(0,N), rlnorm(N, meanlog = log(mulag)*rep(1,N), sdlog = siglag/mulag)), 1, max)
          
          # calculate times to clear
          TC = array(NA, nsen+nres)
          for (i in 1:nsen){
            TC[i]=(hsen[i]*(log(10^ip[i])-log(dlim))/(log(2)))+lag[i]
          }
          for (i in 1:nres){
            TC[nsen+i]=(hres[i]*(log(10^ip[nsen+i])-log(dlim))/(log(2)))+lag[i]
          }
          
          # DAY 1  
          nTC1 = array(NA, N)
          tob1 = rlnorm(1, meanlog = log(24), sdlog = 3.5/24)
          for (i in 1:N){
            if (TC[i]>=rnorm(1, mean = tob1, sd = 6)){ 
              nTC1[i]=1
            }else{	
              nTC1[i]=0	
            }
          }
          perc1day[rr,s]=sum(nTC1)
          
          # DAY 2
          nTC2 = array(NA, N)
          tob2=rlnorm(1, meanlog = log(48), sdlog = 3.5/48)
          for (i in 1:N){
            if (TC[i]>=rnorm(1, mean = tob2, sd = 6)){
              nTC2[i]=1
            }else{
              nTC2[i]=0
            }
          }
          perc2day[rr,s]=sum(nTC2) 
          
          # DAY 3
          nTC3 = array(NA, N)
          tob3=rlnorm(1, meanlog = log(72), sdlog = 3.5/72)
          for (i in 1:N){
            if (TC[i]>=rnorm(1, mean = tob3, sd = 6)){
              nTC3[i]=1
            }else{
              nTC3[i]=0
            }
          }
          perc3day[rr,s]=sum(nTC3)
        }
      }
      #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
      plot(100*apply(perc1day, 1, quantile, prob = cil/10)/N,percres,col = "black", lty = "dashed", lwd = 5, typ = "l", xlim = c(0,100),ylim = c(0,100),main = sprintf("N=%s, P0 = %s, Day 1", ipdat[nn,2], round(exp(ipdat[nn,3]))), xlab = "", ylab = "")
      lines(100*apply(perc1day, 1, quantile, prob = 0.5)/N,percres,col = "black", lty = "solid", lwd = 5, typ = "l")
      lines(100*apply(perc1day, 1, quantile, prob = ciu/100)/N,percres,col = "black", lty = "dashed", lwd = 5, typ = "l")
      lines(c(ipdat[nn,5],ipdat[nn,5]),c(0,100),col = "red", lty = "dashed", lwd = 5, typ = "l")
      #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
      plot(100*apply(perc2day, 1, quantile, prob = cil/100)/N,percres,col = "black", lty = "dashed", lwd = 5, typ = "l", xlim = c(0,100),ylim = c(0,100),main = sprintf("N=%s, P0 = %s, Day 2", ipdat[nn,2], round(exp(ipdat[nn,3]))),xlab = "", ylab = "")
      lines(100*apply(perc2day, 1, quantile, prob = 0.5)/N,percres,col = "black", lty = "solid", lwd = 5, typ = "l")
      lines(100*apply(perc2day, 1, quantile, prob = ciu/100)/N,percres,col = "black", lty = "dashed", lwd = 5, typ = "l")
      lines(c(ipdat[nn,6],ipdat[nn,6]),c(0,100),col = "red", lty = "dashed", lwd = 5, typ = "l")
      #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
      plot(100*apply(perc3day, 1, quantile, prob = cil/100)/N,percres,col = "black", lty = "dashed", lwd = 5, typ = "l", xlim = c(0,100),ylim = c(0,100),main = sprintf("N=%s, P0 = %s, Day 3", ipdat[nn,2], round(exp(ipdat[nn,3]))), xlab = "", ylab = "")
      lines(100*apply(perc3day, 1, quantile, prob = 0.5)/N,percres,col = "black", lty = "solid", lwd = 5, typ = "l")
      lines(100*apply(perc3day, 1, quantile, prob = ciu/100)/N,percres,col = "black", lty = "dashed", lwd = 5, typ = "l")
      lines(c(ipdat[nn,7],ipdat[nn,7]),c(0,100),col = "red", lty = "dashed", lwd = 5, typ = "l")
      #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
    }
  }
  pres3D[nn,1]<-approx(100*apply(perc3day, 1, quantile, prob = cil/100)/N,y=percres,xout=ipdat[nn,7])$y
  pres3D[nn,2]<-approx(100*apply(perc3day, 1, quantile, prob = 50/100)/N,y=percres,xout=ipdat[nn,7])$y
  pres3D[nn,3]<-approx(100*apply(perc3day, 1, quantile, prob = ciu/100)/N,y=percres,xout=ipdat[nn,7])$y
}

for(i in 1:length(pres3D[,1])){
  if(is.na(pres3D[i,1])==TRUE){
    pres3D[i,1]<-100
  }
  if(is.na(pres3D[i,2])==TRUE){
    if (pres3D[i,1]<100){
      pres3D[i,2]<-0      
    }
    else{
      pres3D[i,2]<-100
    }
        
  }
  if(is.na(pres3D[i,3])==TRUE){
    pres3D[i,3]<-0    
  }
}


# BOTTOM GRAPH OF FIGURE 2 (BLUE LINE) 
plot(seq(2001,2011),pres3D[,2]/100,type="l",ylim=c(0,1),xlab = "Time (years)",ylab="predicted proportion resistant",xlim = c(2001,2012))
par(new=T)
pres3D.low<-pres3D[,3]/100
pres3D.high<-pres3D[,1]/100
polygon(c(seq(2001,2011), rev(seq(2001,2011))), c(pres3D.high, rev(pres3D.low)),col=rgb(0,0,100,50,maxColorValue=255), border = NA,xlim = c(2001,2012))
par(new=T)
plot(seq(2001,2011),pres3D[,2]/100,type="l",ylim=c(0,1),xlab = "Time (years)",ylab="predicted proportion resistant",lwd=5,ps=20,col="blue4",xlim = c(2001,2012))
par(new=T)
plot(seq(2001,2012),pres.med,type="l",ylim=c(0,1),xlab = "Time (years)",ylab="predicted proportion resistant",xlim = c(2001,2012))
par(new=T)
polygon(c(seq(2001,2012), rev(seq(2001,2012))), c(pres.high, rev(pres.low)),col=rgb(100,100,100,70,maxColorValue=255), border = NA,xlim = c(2001,2012))
par(new=T)
plot(seq(2001,2012),pres.med,type="l",ylim=c(0,1),xlab = "Time (years)",ylab="predicted proportion resistant",lwd=5,ps=20,xlim = c(2001,2012))

plot(seq(2001,2011),pres3D[,2]/100,type="l",ylim=c(0,1),xlab = "Time (years)",ylab="predicted proportion resistant")
par(new=T)
pres3D.low<-pres3D[,3]/100
pres3D.high<-pres3D[,1]/100
polygon(c(seq(2001,2011), rev(seq(2001,2011))), c(pres3D.high, rev(pres3D.low)),col = "grey80", border = NA)
par(new=T)
plot(seq(2001,2011),pres3D[,2]/100,type="l",ylim=c(0,1),xlab = "Time (years)",ylab="predicted proportion resistant",lwd=5,ps=20)



}


# TEST APPROACH USING SIMULATED DATA
# ANALYSIS TO ACCOMPANY SUPPORTING INFORMATION FILE 1

mixdat <-read.csv('simulateddata50.csv')  # only samples with clones from SMRU

N<-ncol(mixdat)-1


output.mu <- matrix(NA,nrow=M,ncol=N)
output.sigma <- matrix(NA,nrow=M,ncol=N)
output.lambda <- matrix(NA,nrow=M,ncol=N)
output.loglik <- matrix(NA,nrow=M,ncol=N)
output.mu.se <- matrix(NA,nrow=M,ncol=N)
output.sigma.se <- matrix(NA,nrow=M,ncol=N)
output.lambda.se <- matrix(NA,nrow=M,ncol=N)
AIC<-matrix(0,nrow=M,ncol=N)
AICdelta<-matrix(0,nrow=M,ncol=N)

nb<-na.omit(mixdat[,N+1])

for (i in 1:N){
  # 1 COMPONENT LOG NORMAL
  nmixdat<-na.omit(mixdat[,i])
  lmixdat<- log(nmixdat)
  xll<-fitdistr(lmixdat,"normal")
  output.loglik[1,i]<- xll$loglik
  output.mu[1,i]<-xll$estimate[1]
  output.lambda[1,i]<-1
  output.sigma[1,i]<-xll$estimate[2]
  output.mu.se[1,i]<-xll$sd[1]
  output.sigma.se[1,i]<-xll$sd[2]
  output.lambda.se[1,i]<-0
  AIC[1,i]<-2*(3*1-1)-2*output.loglik[1,i]
  AICdelta[1,i]<-0
}
for (i in 1:N){
  nmixdat<-na.omit(mixdat[,i])
  lmixdat<- log(nmixdat)
  # >=2 COMPONENTS LOG NORMAL
  j<-1
  while((j<=M-1) && AICdelta[j,i]<=pval){
    j<-j+1
    res <- normalmixEM(lmixdat, lambda = matrix((1/j),nrow=1,ncol=j), mu = 2*(1:j)/j, sigma = 0.3*matrix(1,nrow=1,ncol=j))
    resboot <- boot.se(res, B = nboot)
    resboot[c("lambda.se", "mu.se", "sigma.se","loglik.se")]	
    output.loglik[j,i]<-res$loglik
    AIC[j,i]<-2*(3*j-1)-2*output.loglik[j,i]
    AICdelta[j,i]<-exp(-(AIC[j-1,i]-AIC[j,i])/2)
    if(AICdelta[j,i]<=pval){
      output.mu[1:j,i]<-res$mu
      output.sigma[1:j,i]<-res$sigma
      output.lambda[1:j,i]<-res$lambda
      output.mu.se[1:j,i]<-resboot$mu.se
      output.sigma.se[1:j,i]<-resboot$sigma.se
      output.lambda.se[1:j,i]<-resboot$lambda.se			
    }
  }
}


# GRAPH FOR SAMPLE SIZE 50
for (i in 1:6){
  plot((1:9)/10,1-output.lambda[1,((i-1)*9+1):(9*i)],type="l",xlim = c(0,1),ylim=c(0,1),xlab = "input proportion resistant",ylab="predicted proportion resistant",col=rgb(0,(7-i)*30,i*30,150,maxColorValue=255),lwd = 3)
  par(new=T)  
}
for (i in 7:7){
  plot((1:9)/10,1-output.lambda[1,((i-1)*9+1):(9*i)],type="l",xlim = c(0,1),ylim=c(0,1),xlab = "input proportion resistant",ylab="predicted proportion resistant",col=rgb(0,(7-i)*30,i*30,150,maxColorValue=255),lwd = 3,main="sample size = 50")
}
par(new=T)  
plot((1:9)/10,(1:9)/10,type="l",lty=3,lwd = 3,xlim = c(0,1),ylim=c(0,1),xlab = "input proportion resistant",ylab="predicted proportion resistant")


mixdat <-read.csv('simulateddata100.csv')  # only samples with clones from SMRU


N<-ncol(mixdat)-1


output.mu <- matrix(NA,nrow=M,ncol=N)
output.sigma <- matrix(NA,nrow=M,ncol=N)
output.lambda <- matrix(NA,nrow=M,ncol=N)
output.loglik <- matrix(NA,nrow=M,ncol=N)
output.mu.se <- matrix(NA,nrow=M,ncol=N)
output.sigma.se <- matrix(NA,nrow=M,ncol=N)
output.lambda.se <- matrix(NA,nrow=M,ncol=N)
AIC<-matrix(0,nrow=M,ncol=N)
AICdelta<-matrix(0,nrow=M,ncol=N)

nb<-na.omit(mixdat[,N+1])

for (i in 1:N){
  # 1 COMPONENT LOG NORMAL
  nmixdat<-na.omit(mixdat[,i])
  lmixdat<- log(nmixdat)
  xll<-fitdistr(lmixdat,"normal")
  output.loglik[1,i]<- xll$loglik
  output.mu[1,i]<-xll$estimate[1]
  output.lambda[1,i]<-1
  output.sigma[1,i]<-xll$estimate[2]
  output.mu.se[1,i]<-xll$sd[1]
  output.sigma.se[1,i]<-xll$sd[2]
  output.lambda.se[1,i]<-0
  AIC[1,i]<-2*(3*1-1)-2*output.loglik[1,i]
  AICdelta[1,i]<-0
}
for (i in 1:N){
  nmixdat<-na.omit(mixdat[,i])
  lmixdat<- log(nmixdat)
  # >=2 COMPONENTS LOG NORMAL
  j<-1
  while((j<=M-1) && AICdelta[j,i]<=pval){
    j<-j+1
    res <- normalmixEM(lmixdat, lambda = matrix((1/j),nrow=1,ncol=j), mu = 2*(1:j)/j, sigma = 0.3*matrix(1,nrow=1,ncol=j))
    resboot <- boot.se(res, B = nboot)
    resboot[c("lambda.se", "mu.se", "sigma.se","loglik.se")]  
    output.loglik[j,i]<-res$loglik
    AIC[j,i]<-2*(3*j-1)-2*output.loglik[j,i]
    AICdelta[j,i]<-exp(-(AIC[j-1,i]-AIC[j,i])/2)
    if(AICdelta[j,i]<=pval){
      output.mu[1:j,i]<-res$mu
      output.sigma[1:j,i]<-res$sigma
      output.lambda[1:j,i]<-res$lambda
      output.mu.se[1:j,i]<-resboot$mu.se
      output.sigma.se[1:j,i]<-resboot$sigma.se
      output.lambda.se[1:j,i]<-resboot$lambda.se			
    }
  }
}

# GRAPH FOR SAMPLE SIZE 100
for (i in 1:6){
  plot((1:9)/10,1-output.lambda[1,((i-1)*9+1):(9*i)],type="l",xlim = c(0,1),ylim=c(0,1),xlab = "input proportion resistant",ylab="predicted proportion resistant",col=rgb(0,(7-i)*30,i*30,150,maxColorValue=255),lwd = 3)
  par(new=T)  
}
for (i in 7:7){
  plot((1:9)/10,1-output.lambda[1,((i-1)*9+1):(9*i)],type="l",xlim = c(0,1),ylim=c(0,1),xlab = "input proportion resistant",ylab="predicted proportion resistant",col=rgb(0,(7-i)*30,i*30,150,maxColorValue=255),lwd = 3,main="sample size = 100")
}
par(new=T)  
plot((1:9)/10,(1:9)/10,type="l",lty=3,lwd = 3,xlim = c(0,1),ylim=c(0,1),xlab = "input proportion resistant",ylab="predicted proportion resistant")



mixdat <-read.csv('simulateddata200.csv')  # only samples with clones from SMRU


N<-ncol(mixdat)-1


output.mu <- matrix(NA,nrow=M,ncol=N)
output.sigma <- matrix(NA,nrow=M,ncol=N)
output.lambda <- matrix(NA,nrow=M,ncol=N)
output.loglik <- matrix(NA,nrow=M,ncol=N)
output.mu.se <- matrix(NA,nrow=M,ncol=N)
output.sigma.se <- matrix(NA,nrow=M,ncol=N)
output.lambda.se <- matrix(NA,nrow=M,ncol=N)
AIC<-matrix(0,nrow=M,ncol=N)
AICdelta<-matrix(0,nrow=M,ncol=N)

nb<-na.omit(mixdat[,N+1])

for (i in 1:N){
  # 1 COMPONENT LOG NORMAL
  nmixdat<-na.omit(mixdat[,i])
  lmixdat<- log(nmixdat)
  xll<-fitdistr(lmixdat,"normal")
  output.loglik[1,i]<- xll$loglik
  output.mu[1,i]<-xll$estimate[1]
  output.lambda[1,i]<-1
  output.sigma[1,i]<-xll$estimate[2]
  output.mu.se[1,i]<-xll$sd[1]
  output.sigma.se[1,i]<-xll$sd[2]
  output.lambda.se[1,i]<-0
  AIC[1,i]<-2*(3*1-1)-2*output.loglik[1,i]
  AICdelta[1,i]<-0
}
for (i in 1:N){
  nmixdat<-na.omit(mixdat[,i])
  lmixdat<- log(nmixdat)
  # >=2 COMPONENTS LOG NORMAL
  j<-1
  while((j<=M-1) && AICdelta[j,i]<=pval){
    j<-j+1
    res <- normalmixEM(lmixdat, lambda = matrix((1/j),nrow=1,ncol=j), mu = 2*(1:j)/j, sigma = 0.3*matrix(1,nrow=1,ncol=j))
    resboot <- boot.se(res, B = nboot)
    resboot[c("lambda.se", "mu.se", "sigma.se","loglik.se")]  
    output.loglik[j,i]<-res$loglik
    AIC[j,i]<-2*(3*j-1)-2*output.loglik[j,i]
    AICdelta[j,i]<-exp(-(AIC[j-1,i]-AIC[j,i])/2)
    if(AICdelta[j,i]<=pval){
      output.mu[1:j,i]<-res$mu
      output.sigma[1:j,i]<-res$sigma
      output.lambda[1:j,i]<-res$lambda
      output.mu.se[1:j,i]<-resboot$mu.se
      output.sigma.se[1:j,i]<-resboot$sigma.se
      output.lambda.se[1:j,i]<-resboot$lambda.se  		
    }
  }
}

# GRAPH FOR SAMPLE SIZE 200
for (i in 1:6){
  plot((1:9)/10,1-output.lambda[1,((i-1)*9+1):(9*i)],type="l",xlim = c(0,1),ylim=c(0,1),xlab = "input proportion resistant",ylab="predicted proportion resistant",col=rgb(0,(7-i)*30,i*30,150,maxColorValue=255),lwd = 3)
  par(new=T)  
}
for (i in 7:7){
  plot((1:9)/10,1-output.lambda[1,((i-1)*9+1):(9*i)],type="l",xlim = c(0,1),ylim=c(0,1),xlab = "input proportion resistant",ylab="predicted proportion resistant",col=rgb(0,(7-i)*30,i*30,150,maxColorValue=255),lwd = 3,main="sample size = 200")
}
par(new=T)  
plot((1:9)/10,(1:9)/10,type="l",lty=3,lwd = 3,xlim = c(0,1),ylim=c(0,1),xlab = "input proportion resistant",ylab="predicted proportion resistant")


mixdat <-read.csv('simulateddata1000.csv')  # only samples with clones from SMRU


N<-ncol(mixdat)-1


output.mu <- matrix(NA,nrow=M,ncol=N)
output.sigma <- matrix(NA,nrow=M,ncol=N)
output.lambda <- matrix(NA,nrow=M,ncol=N)
output.loglik <- matrix(NA,nrow=M,ncol=N)
output.mu.se <- matrix(NA,nrow=M,ncol=N)
output.sigma.se <- matrix(NA,nrow=M,ncol=N)
output.lambda.se <- matrix(NA,nrow=M,ncol=N)
AIC<-matrix(0,nrow=M,ncol=N)
AICdelta<-matrix(0,nrow=M,ncol=N)

nb<-na.omit(mixdat[,N+1])

for (i in 1:N){
  # 1 COMPONENT LOG NORMAL
  nmixdat<-na.omit(mixdat[,i])
  lmixdat<- log(nmixdat)
  xll<-fitdistr(lmixdat,"normal")
  output.loglik[1,i]<- xll$loglik
  output.mu[1,i]<-xll$estimate[1]
  output.lambda[1,i]<-1
  output.sigma[1,i]<-xll$estimate[2]
  output.mu.se[1,i]<-xll$sd[1]
  output.sigma.se[1,i]<-xll$sd[2]
  output.lambda.se[1,i]<-0
  AIC[1,i]<-2*(3*1-1)-2*output.loglik[1,i]
  AICdelta[1,i]<-0
}
for (i in 1:N){
  nmixdat<-na.omit(mixdat[,i])
  lmixdat<- log(nmixdat)
  # >=2 COMPONENTS LOG NORMAL
  j<-1
  while((j<=M-1) && AICdelta[j,i]<=pval){
    j<-j+1
    res <- normalmixEM(lmixdat, lambda = matrix((1/j),nrow=1,ncol=j), mu = 2*(1:j)/j, sigma = 0.3*matrix(1,nrow=1,ncol=j))
    resboot <- boot.se(res, B = nboot)
    resboot[c("lambda.se", "mu.se", "sigma.se","loglik.se")]  
    output.loglik[j,i]<-res$loglik
    AIC[j,i]<-2*(3*j-1)-2*output.loglik[j,i]
    AICdelta[j,i]<-exp(-(AIC[j-1,i]-AIC[j,i])/2)
    if(AICdelta[j,i]<=pval){
      output.mu[1:j,i]<-res$mu
      output.sigma[1:j,i]<-res$sigma
      output.lambda[1:j,i]<-res$lambda
      output.mu.se[1:j,i]<-resboot$mu.se
      output.sigma.se[1:j,i]<-resboot$sigma.se
      output.lambda.se[1:j,i]<-resboot$lambda.se    	
    }
  }
}

# GRAPH FOR SAMPLE SIZE 1000
for (i in 1:6){
  plot((1:9)/10,1-output.lambda[1,((i-1)*9+1):(9*i)],type="l",xlim = c(0,1),ylim=c(0,1),xlab = "input proportion resistant",ylab="predicted proportion resistant",col=rgb(0,(7-i)*30,i*30,150,maxColorValue=255),lwd = 3)
  par(new=T)  
}
for (i in 7:7){
  plot((1:9)/10,1-output.lambda[1,((i-1)*9+1):(9*i)],type="l",xlim = c(0,1),ylim=c(0,1),xlab = "input proportion resistant",ylab="predicted proportion resistant",col=rgb(0,(7-i)*30,i*30,150,maxColorValue=255),lwd = 3,main="sample size = 1000")
}
par(new=T)  
plot((1:9)/10,(1:9)/10,type="l",lty=3,lwd = 3,xlim = c(0,1),ylim=c(0,1),xlab = "input proportion resistant",ylab="predicted proportion resistant")



# legend for graphs
for (i in 1:6){
  plot((1:9)/10,((1:9)/(1:9))*(i-0.5)/10,type="l",xlim = c(0,1),ylim=c(0,1),col=rgb(0,(7-i)*30,i*30,150,maxColorValue=255),lwd = 3)
  par(new=T)  
}
for (i in 7:7){
  plot((1:9)/10,((1:9)/(1:9))*(i-0.5)/10,type="l",xlim = c(0,1),ylim=c(0,1),col=rgb(0,(7-i)*30,i*30,150,maxColorValue=255),lwd = 3)
}

# simulation experiment to test model selection algorithm forpredicting the number of components

nit<-5*1000 #number of mixture distributions for testing (must be divisible by 5)
pick<-matrix(NA,nrow=nit,ncol=20) # output matrix
maxcomp<-5 # maximum number of components
samplesize<-1000 # sample size for each distribution 
M<-maxcomp
pval<-0.1 #significance level for comparing models
nboot<-10 # number of iterations for bootstrap
maxmean<-10

# INPUTS
# col1 - mean 1
# col2 - mean 2
# col3 - mean 3
# col4 - mean 4
# col5 - mean 5
# col6 - propotion 1
# col7 - propotion 2
# col8 - propotion 3
# col9 - propotion 4
# col10 - propotion 5
# OUTPUTS
# col11 - mean 1
# col12 - mean 2
# col13 - mean 3
# col14 - mean 4
# col15 - mean 5
# col16 - propotion 1
# col17 - propotion 2
# col18 - propotion 3
# col19 - propotion 4
# col20 - propotion 5

# randomly select mean for each component distribution
for (i in 1:nit){
  k<-sample.int(5, size = 5, replace = FALSE, prob = NULL)
  for (j in 1:5){
    pick[i,j]<-runif(1,min=(k[j]-1)*maxmean/5,max=k[j]*maxmean/5)
  }
}

# randomly select an sd for each set of distributions
sdev<-matrix(NA,nrow=nit,ncol=1)
sdev<-runif(nit,min=0,max=0.2)

# number of components for each distribution
ncomp<-matrix(NA,nrow=1,ncol=nit)
for (j in 1:5){
  ncomp[((j-1)*(nit/5)+1):(j*(nit/5))]<-j  
}

# randomly select proportion for each component
pick[,6:10]<-0
for (i in 1:nit){
  if (ncomp[i]==1){
    pick[i,6]<-1
  }
  if (ncomp[i]>1){
    pick[i,6:(5+ncomp[i])]<-runif(ncomp[i],min=0,max=1)
    pick[i,6:(5+ncomp[i])]<-pick[i,6:(5+ncomp[i])]/sum(pick[i,6:(5+ncomp[i])])
  } 
}

mixdat<-matrix(NA,nrow=(samplesize+1),ncol=(nit+1)) #matrix of input distribution samples
mixdat[2:(nit+1),(nit+1)]<-samplesize

q<-matrix(NA, nrow=nit,ncol=samplesize)
pp<-matrix(NA,nrow=nit,ncol=5)
qmean<-matrix(NA, nrow=nit,ncol=samplesize)

for (i in 1:nit){
  for (k in 1:5){
    pp[i,k]<-sum(pick[i,6:(5+k)])
  }
  for (j in 1: samplesize){
    q[i,j]<-runif(1,min=0,max=1)
    qmean[i,j]<-(q[i,j]<=pp[i,1])*pick[i,1]+(q[i,j]>pp[i,1])*(q[i,j]<=pp[i,2])*pick[i,2]+(q[i,j]>pp[i,2])*(q[i,j]<=pp[i,3])*pick[i,3]+(q[i,j]>pp[i,3])*(q[i,j]<=pp[i,4])*pick[i,4]+(q[i,j]>pp[i,4])*pick[i,5]
    mixdat[j+1,i]<-exp(rnorm(1,mean=log(qmean[i,j]),sd=sdev[i]))
  }
}

N<-ncol(mixdat)-1

output.mu <- matrix(NA,nrow=M,ncol=N)
output.sigma <- matrix(NA,nrow=M,ncol=N)
output.lambda <- matrix(NA,nrow=M,ncol=N)
output.loglik <- matrix(NA,nrow=M,ncol=N)
output.mu.se <- matrix(NA,nrow=M,ncol=N)
output.sigma.se <- matrix(NA,nrow=M,ncol=N)
output.lambda.se <- matrix(NA,nrow=M,ncol=N)
AIC<-matrix(0,nrow=M,ncol=N)
AICdelta<-matrix(0,nrow=M,ncol=N)

nb<-na.omit(mixdat[,N+1])

for (i in 1:N){
  # 1 COMPONENT LOG NORMAL
  nmixdat<-na.omit(mixdat[,i])
  lmixdat<- log(nmixdat)
  xll<-fitdistr(lmixdat,"normal")
  output.loglik[1,i]<- xll$loglik
  output.mu[1,i]<-xll$estimate[1]
  output.lambda[1,i]<-1
  output.sigma[1,i]<-xll$estimate[2]
  output.mu.se[1,i]<-xll$sd[1]
  output.sigma.se[1,i]<-xll$sd[2]
  output.lambda.se[1,i]<-0
  AIC[1,i]<-2*(3*1-1)-2*output.loglik[1,i]
  AICdelta[1,i]<-0
}
for (i in 1:N){
  nmixdat<-na.omit(mixdat[,i])
  lmixdat<- log(nmixdat)
  # >=2 COMPONENTS LOG NORMAL
  j<-1
  while((j<=M-1) && AICdelta[j,i]<=pval){
    j<-j+1
    res <- normalmixEM(lmixdat, lambda = matrix((1/j),nrow=1,ncol=j), mu = 2*(1:j)/j, sigma = 0.3*matrix(1,nrow=1,ncol=j))
    resboot <- boot.se(res, B = nboot)
    resboot[c("lambda.se", "mu.se", "sigma.se","loglik.se")]  
    output.loglik[j,i]<-res$loglik
    AIC[j,i]<-2*(3*j-1)-2*output.loglik[j,i]
    AICdelta[j,i]<-exp(-(AIC[j-1,i]-AIC[j,i])/2)
    if(AICdelta[j,i]<=pval){
      output.mu[1:j,i]<-res$mu
      output.sigma[1:j,i]<-res$sigma
      output.lambda[1:j,i]<-res$lambda
      output.mu.se[1:j,i]<-resboot$mu.se
      output.sigma.se[1:j,i]<-resboot$sigma.se
      output.lambda.se[1:j,i]<-resboot$lambda.se  		
    }
  }
}

pick[,11:15]<-t(exp(output.mu))
pick[,16:20]<-t(output.lambda)

# assess model performance

predncomp<-matrix(NA,nrow=nit,ncol=1)
for (i in 1:nit){
  predncomp[i]<-5-sum(is.na(pick[i,11:16]))
}

predncompmed<-matrix(NA,nrow=1,ncol=5)
predncomplow<-matrix(NA,nrow=1,ncol=5)
predncomphigh<-matrix(NA,nrow=1,ncol=5)
predncompmean<-matrix(NA,nrow=1,ncol=5)
perccorrect<-matrix(NA,nrow=1,ncol=5)

accncomp<-t(ncomp)==predncomp

# plotting

for (j in 1:5){
  predncompmed[j]<-quantile(predncomp[((j-1)*(nit/5)+1):(j*(nit/5))],probs=0.5)
  predncomplow[j]<-quantile(predncomp[((j-1)*(nit/5)+1):(j*(nit/5))],probs=0.05)
  predncomphigh[j]<-quantile(predncomp[((j-1)*(nit/5)+1):(j*(nit/5))],probs=0.95) 
  predncompmean[j]<-mean(predncomp[((j-1)*(nit/5)+1):(j*(nit/5))])
  perccorrect[j]<-100*mean(accncomp[((j-1)*(nit/5)+1):(j*(nit/5))])
}


plot(seq(1:5),perccorrect, type="l",lwd=5,xlim=c(1,5),ylim=c(0,100),xlab = "True number of components",ylab="% correct prediction")

plot(0,xlim=c(1,5),xaxt = "n",ylim=c(1,5),xlab = "True number of components",ylab="predicted number of components")
axis(1,at=1:5,1:5)
par(new=T)
polygon(c(seq(1,5), rev(seq(1,5))), c(predncomphigh, rev(predncomplow)),col = "grey80", border = NA)
lines(seq(1:5),predncompmed,lwd=5,ps=20)
lines(seq(1:5),predncompmean,lwd=3,ps=20, col="grey10", lty=3)


for (j in 1:5){
  predncompmed[j]<-quantile(predncomp[((j-1)*(nit/5)+1):(j*(nit/5))],probs=0.5)
  predncomplow[j]<-quantile(predncomp[((j-1)*(nit/5)+1):(j*(nit/5))],probs=0.25)
  predncomphigh[j]<-quantile(predncomp[((j-1)*(nit/5)+1):(j*(nit/5))],probs=0.75) 
  predncompmean[j]<-mean(predncomp[((j-1)*(nit/5)+1):(j*(nit/5))])
  perccorrect[j]<-100*mean(accncomp[((j-1)*(nit/5)+1):(j*(nit/5))])
}


polygon(c(seq(1,5), rev(seq(1,5))), c(predncomphigh, rev(predncomplow)),col = "grey50", border = NA)
lines(seq(1:5),predncompmed,lwd=5,ps=20)
lines(seq(1:5),predncompmean,lwd=3,ps=20, col="grey10", lty=3)
