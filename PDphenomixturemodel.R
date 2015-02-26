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


#### SMRU BY YEAR ####
mixdat <-read.csv('simulated_cloneHLdata_SMRUbyyear.csv')  # only samples with clones from SMRU
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

# HISTOGRAM AND PDF
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

# PLOT MEAN OVER TIME
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



plot(seq(2001,2012),1-output.lambda[1,],type="l",ylim=c(0,1),xlab = "Time (years)",ylab="predicted proportion resistant")
par(new=T)
pres.low<-(1-output.lambda[1,])-output.lambda.se[1,]
pres.high<-(1-output.lambda[1,])+output.lambda.se[1,]
#polygon(c(seq(2001,2012), rev(seq(2001,2012))), c(pres.high, rev(pres.low)),col = "grey80", border = "black", lty="dashed")
polygon(c(seq(2001,2012), rev(seq(2001,2012))), c(pres.high, rev(pres.low)),col = "grey80", border = NA)
par(new=T)
plot(seq(2001,2012),1-output.lambda[1,],type="l",ylim=c(0,1),xlab = "Time (years)",ylab="predicted proportion resistant",lwd=5,ps=20)
pres.med<-1-output.lambda[1,]
