BSMOption<-function(PutCallFlag,S,T,K,r,sigma,Greek)
{
  if(PutCallFlag=="Call" && Greek=="Price" && T>0)
  {
    d1<-(log(S/K)+(r+sigma^2/2)*T) / (sigma*sqrt(T))
    d2<-d1 - sigma * sqrt(T)
    return(S * pnorm(d1, mean=0, sd=1) - K * exp(-r*T) * pnorm(d2, mean=0, sd=1))
  }
  if(PutCallFlag=="Call" && Greek=="Price" && T==0)
  {return(pmax(S-K,0))}
  if(PutCallFlag=="Call" && Greek=="Delta" && T>0)
  {
    d1<-(log(S/K)+(r+sigma^2/2)*T) / (sigma*sqrt(T))
    return(pnorm(d1, mean=0, sd=1))
  }
  if(PutCallFlag=="Call" && Greek=="Delta" && T==0)
  {return(0)}
  if(PutCallFlag=="Put" && T>0)
  {return(-1 * BSMOption("Call",S,T,K,r,-1 * sigma,Greek))}
  if(PutCallFlag=="Put" && Greek=="Price" && T==0)
  {return(pmax(K-S,0))}
}

### Parameter ###
PutCallFlag<-"Call" # Put or Call
S_0<-100
T<-1
K<-S_0 * 1.15
r<-0.05
mu<-0.1
sigma<-0.2
Npaths<-10000
Nhedgepoints<-252

dt<-T/Nhedgepoints
set.seed(123456)

# Results
S_t<-rep(0, times=Npaths)
Payoff<-rep(0, times=Npaths)
Hedge<-rep(0, times=Npaths)
Error<-rep(0, times=Npaths)

### Simulation Engine ###
# Loop over Npaths
for(j in 1:Npaths) 
{
  S<-S_0
  
  # DeltaHedge
  V<-BSMOption(PutCallFlag,S,T,K,r,sigma,"Price") # initial investment
  a<-BSMOption(PutCallFlag,S,T,K,r,sigma,"Delta") # stock position<-delta
  b<-V-a*S # rest in bank, self-fin. Cond.
  
  # Loop over Nhedgepoints
  for(i in 1:Nhedgepoints) 
  {
    eps<-rnorm(1) # eps=random.normal() simulate outcome of N(0,1)
    S<-S*exp((mu-0.5*sigma^2)*dt+sigma*sqrt(dt)*eps)
    
    # Delta-Hedge
    V<-a*S+b*exp(r*dt)
    if((PutCallFlag=="Call" && K<S) || ((PutCallFlag=="Put" && K>S))){
        a<-BSMOption(PutCallFlag,S,T-i*dt,K,r,sigma,"Delta")
    } else {
      a<-0
    }
    
    b<-V-a*S
  }
  
  S_t[j]<-S
  Hedge[j]<-V
  Payoff[j]<-BSMOption(PutCallFlag,S,0,K,r,sigma,"Price")
}

# Calc Error
Error<-Payoff-Hedge

### Plots ###
# Define Scatterplot PayOff vs Hedge
layout(rbind(c(1)))

# Plot PayOff vs Hedge
plot(S_t, Hedge, col="red", type="p", pch='.', cex=2, main="Results Hedging", ylab="Hedge / PayOff", xlab="")
points(S_t, Payoff, col='green', type="p", pch='.', cex=1)

# Add legend 
legend(ifelse(PutCallFlag=="Call","topleft","topright"), legend=c("Payoff","Hedge"), col=c("green", "red"), pch=c(19,19), cex=1)

# Define Histogramms Errors and Boxplot Errors
par(mfrow=c(2,1))

# Plot Hist Error & Boxplot
hist(Error, breaks=nclass.FD(Error), freq=F, col="gray", border="black", main='', xlab='Error Heding', ylab='') 
boxplot(Error, horizontal=TRUE, main='')
