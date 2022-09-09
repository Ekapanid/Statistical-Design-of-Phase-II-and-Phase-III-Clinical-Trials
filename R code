
install.packages("clinfun")
library(clinfun)

#Simon two stage design for 20% diff power 0.8
ph2simon(0.05,0.25,0.05,0.2)
ph2simon(0.10,0.30,0.05,0.2)
ph2simon(0.20,0.40,0.05,0.2)
ph2simon(0.3,0.5,0.05,0.2)
ph2simon(0.4,0.6,0.05,0.2)
ph2simon(0.6,0.8,0.05,0.2)
ph2simon(0.7,0.9,0.05,0.2)

#Simon two stage design for 20% diff power 0.9
ph2simon(0.05,0.25,0.05,0.1)
ph2simon(0.10,0.30,0.05,0.1)
ph2simon(0.20,0.40,0.05,0.1)
ph2simon(0.3,0.5,0.05,0.1)
ph2simon(0.4,0.6,0.05,0.1)
ph2simon(0.6,0.8,0.05,0.1)
ph2simon(0.7,0.9,0.05,0.1)


#Simon two stage design for 15% diff power 0.8
ph2simon(0.05,0.2,0.05,0.2)
ph2simon(0.10,0.25,0.05,0.2)
ph2simon(0.20,0.35,0.05,0.2)
ph2simon(0.30,0.45,0.05,0.2)
ph2simon(0.40,0.55,0.05,0.2)
ph2simon(0.50,0.65,0.05,0.2)
ph2simon(0.6,0.75,0.05,0.2)
ph2simon(0.7,0.85,0.05,0.2)
ph2simon(0.8,0.95,0.05,0.2)


#Simon two stage design for 15% diff power 0.9
ph2simon(0.05,0.2,0.05,0.1)
ph2simon(0.10,0.25,0.05,0.1)
ph2simon(0.20,0.35,0.05,0.1)
ph2simon(0.30,0.45,0.05,0.1)
ph2simon(0.40,0.55,0.05,0.1)
ph2simon(0.50,0.65,0.05,0.1)
ph2simon(0.6,0.75,0.05,0.1)
ph2simon(0.7,0.85,0.05,0.1)
ph2simon(0.8,0.95,0.05,0.1)



#######################################################################################


#Simon Pick the Winner Design Sample Size n=44 per arm K=3 arms 15% diff
pselect(44,c(0.2,0.2,0.35))
pselect(44,c(0.3,0.3,0.45))
pselect(44,c(0.4,0.2,0.35))
pselect(44,c(0.4,0.2,0.35))
pselect(44,c(0.4,0.2,0.35))
pselect(44,c(0.4,0.2,0.35))


#This is the probability of selecting the correct treatment without ties


#With the code below we create a new function that adds to the previous result the probability
#of selecting the correct treatment even when there are ties.
i<-0
j<-1:2
x<-0

sum(dbinom(i,44,0.35))

length(i)


#probability for ties can be written as sum(choose(2,j)*dbinom(i,44,0.20)**j*pbinom(i-1,44,0.20)**(2-j)/(j+1)) for specific i

for (i in 0:44 ){
  x<-x+dbinom(i,44,0.35)*sum(choose(2,j)*(dbinom(i,44,0.20)**j)*(pbinom(i-1,44,0.20)**(2-j))/(j+1))
}

x

simon_rand<-function(n,prob){
  result<-pselect(n,prob)
  prob_tot<-as.numeric(result$prob.selection[length(prob),length(prob)])
  j<-1:(length(prob)-1)
  for (i in 0:44 ){
    prob_tot<-prob_tot+dbinom(i,44,0.35)*sum(choose(2,j)*(dbinom(i,44,0.20)**j)*(pbinom(i-1,44,0.20)**(2-j))/(j+1))
  }
  return(prob_tot)
}



#Thus the total probability of correct selection is :
simon_rand(44,c(0.2,0.2,0.35))
simon_rand(44,c(0.3,0.3,0.45))
simon_rand(44,c(0.4,0.4,0.55))
simon_rand(44,c(0.5,0.5,0.65))
simon_rand(44,c(0.6,0.6,0.75))
simon_rand(44,c(0.7,0.7,0.85))
simon_rand(44,c(0.8,0.8,0.95))



#######################################################################################


#Distribution of MLE and UMVUE null p=0.3 alternative p=0.5 actual p=0.5

z<-rep(0,1000)
c<-rep(0,1000)
r<-rep(0,1000)
t<-rep(0,1000)



for (i in 1:1000){
  
  x<-rbinom(n=1,size=15,prob=0.50)
  y<-rbinom(n=1,size=31,prob=0.50)
  
  
  c[i]<-if(x<=5){
    c[i]<-x/15} else{
      c[i]<-(x+y)/46
    }
  
  
  z[i]<-twostage.inference(x+y,5,15,46,0.30,0.05)[1]
  
  x<-rbinom(n=1,size=19,prob=0.50)
  y<-rbinom(n=1,size=20,prob=0.50)
  
  t[i]<-if(x<=6){
    t[i]<-x/19} else{
      t[i]<-(x+y)/39
    }
  
  
  r[i]<-twostage.inference(x+y,6,19,39,0.30,0.05)[1]
  
}


par(mfrow=c(1,4))
hist(c,main='Distribution of MLE (Optimal)',xlab = 'Response Rate estimate')
hist(t,main='Distribution of MLE (Minimax)',xlab = 'Response Rate estimate')
hist(z,main='Distribution of UMVUE (Optimal)',xlab = 'Response Rate estimate')
hist(r,main='Distribution of UMVUE (Minimax)',xlab = 'Response Rate estimate')



#Distribution of MLE and UMVUE null p=0.3 alternative p=0.5 actual p=0.4


z<-rep(0,1000)
c<-rep(0,1000)
r<-rep(0,1000)
t<-rep(0,1000)




for (i in 1:1000){
  
  x<-rbinom(n=1,size=15,prob=0.40)
  y<-rbinom(n=1,size=31,prob=0.40)
  
  
  c[i]<-if(x<=5){
    c[i]<-x/15} else{
      c[i]<-(x+y)/46
    }
  
  
  z[i]<-twostage.inference(x+y,5,15,46,0.30,0.05)[1]
  
  x<-rbinom(n=1,size=19,prob=0.40)
  y<-rbinom(n=1,size=20,prob=0.40)
  
  t[i]<-if(x<=6){
    t[i]<-x/19} else{
      t[i]<-(x+y)/39
    }
  
  
  r[i]<-twostage.inference(x+y,6,19,39,0.30,0.05)[1]
  
}



par(mfrow=c(1,4))
hist(c,main='Distribution of MLE(Optimal)',xlab = 'Response Rate estimate')
hist(t,main='Distribution of MLE (Minimax)',xlab = 'Response Rate estimate')
hist(z,main='Distribution of UMVUE (Optimal)',xlab = 'Response Rate estimate')
hist(r,main='Distribution of UMVUE (Minimax)',xlab = 'Response Rate estimate')


#Distribution of MLE and UMVUE null p=0.3 alternative p=0.5 actual p=0.3


z<-rep(0,1000)
c<-rep(0,1000)
r<-rep(0,1000)
t<-rep(0,1000)




for (i in 1:1000){
  
  x<-rbinom(n=1,size=15,prob=0.30)
  y<-rbinom(n=1,size=31,prob=0.30)
  
  
  c[i]<-if(x<=5){
    c[i]<-x/15} else{
      c[i]<-(x+y)/46
    }
  
  
  z[i]<-twostage.inference(x+y,5,15,46,0.30,0.05)[1]
  
  x<-rbinom(n=1,size=19,prob=0.30)
  y<-rbinom(n=1,size=20,prob=0.30)
  
  t[i]<-if(x<=6){
    t[i]<-x/19} else{
      t[i]<-(x+y)/39
    }
  
  
  r[i]<-twostage.inference(x+y,6,19,39,0.30,0.05)[1]
  
}





par(mfrow=c(1,4))
hist(c,main='Distribution of MLE (Optimal)',xlab = 'Response Rate estimate')
hist(t,main='Distribution of MLE (Minimax)',xlab = 'Response Rate estimate')
hist(z,main='Distribution of UMVUE (Optimal)',xlab = 'Response Rate estimate')
hist(r,main='Distribution of UMVUE (Minimax)',xlab = 'Response Rate estimate')




#Distribution of MLE and UMVUE null p=0.3 alternative p=0.5 actual p=0.2

z<-rep(0,1000)
c<-rep(0,1000)
r<-rep(0,1000)
t<-rep(0,1000)


for (i in 1:1000){
  
  x<-rbinom(n=1,size=15,prob=0.20)
  y<-rbinom(n=1,size=31,prob=0.20)
  
  
  c[i]<-if(x<=5){
    c[i]<-x/15} else{
      c[i]<-(x+y)/46
    }
  
  
  z[i]<-twostage.inference(x+y,5,15,46,0.30,0.05)[1]
  
  x<-rbinom(n=1,size=19,prob=0.20)
  y<-rbinom(n=1,size=20,prob=0.20)
  
  t[i]<-if(x<=6){
    t[i]<-x/19} else{
      t[i]<-(x+y)/39
    }
  
  
  r[i]<-twostage.inference(x+y,6,19,39,0.30,0.05)[1]
  
}

par(mfrow=c(1,4))
hist(c,main='Distribution of MLE (Optimal)',xlab = 'Response Rate estimate')
hist(t,main='Distribution of MLE (Minimax)',xlab = 'Response Rate estimate')
hist(z,main='Distribution of UMVUE (Optimal)',xlab = 'Response Rate estimate')
hist(r,main='Distribution of UMVUE (Minimax)',xlab = 'Response Rate estimate')

#Distribution of MLE and UMVUE null p=0.3 alternative p=0.5 actual p=0.1

z<-rep(0,1000)
c<-rep(0,1000)
r<-rep(0,1000)
t<-rep(0,1000)





for (i in 1:1000){
  
  x<-rbinom(n=1,size=15,prob=0.10)
  y<-rbinom(n=1,size=31,prob=0.10)
  
  
  c[i]<-if(x<=5){
    c[i]<-x/15} else{
      c[i]<-(x+y)/46
    }
  
  
  z[i]<-twostage.inference(x+y,5,15,46,0.30,0.05)[1]
  
  x<-rbinom(n=1,size=19,prob=0.10)
  y<-rbinom(n=1,size=20,prob=0.10)
  
  t[i]<-if(x<=6){
    t[i]<-x/19} else{
      t[i]<-(x+y)/39
    }
  
  
  r[i]<-twostage.inference(x+y,6,19,39,0.30,0.05)[1]
  
}


par(mfrow=c(1,4))
hist(c,breaks=5,main='Distribution of MLE (Optimal)',xlab = 'Response Rate estimate')
hist(t,breaks=5,main='Distribution of MLE (Minimax)',xlab = 'Response Rate estimate')
hist(z,breaks=5,main='Distribution of UMVUE (Optimal)',xlab = 'Response Rate estimate')
hist(r,breaks=5,main='Distribution of UMVUE (Minimax)',xlab = 'Response Rate estimate')









######################################################################################


#Distribution of MLE and UMVUE null p=0.2 alternative p=0.4 actual p=0.5

for (i in 1:1000){
  x<-rbinom(n=1,size=13,prob=0.50)
  y<-rbinom(n=1,size=30,prob=0.50)
  z[i]<-twostage.inference(x+y,3,13,43,0.20,0.05)[1]
}


for (i in 1:1000){
  x<-rbinom(n=1,size=13,prob=0.50)
  y<-rbinom(n=1,size=30,prob=0.50)
  r[i]<-twostage.inference(x+y,3,13,43,0.20,0.05)[1]
}



for (i in 1:1000){
  x<-rbinom(n=1,size=13,prob=0.50)
  y<-rbinom(n=1,size=30,prob=0.50)
  c[i]<-if(x<4){
    c[i]<-x/18} else{
      c[i]<-(x+y)/33
    }
}

hist(c)
hist(z)

par(mfrow=c(1,3))
hist(c,main='Distribution of MLE',xlab = 'Response Rate estimate')
hist(z,main='Distribution of UMVUE (Optimal)',xlab = 'Response Rate estimate')
hist(r,main='Distribution of UMVUE (Minimax)',xlab = 'Response Rate estimate')


#Distribution of MLE and UMVUE null p=0.2 alternative p=0.4 actual p=0.4

for (i in 1:1000){
  x<-rbinom(n=1,size=13,prob=0.40)
  y<-rbinom(n=1,size=30,prob=0.40)
  z[i]<-twostage.inference(x+y,3,13,43,0.20,0.05)[1]
}


for (i in 1:1000){
  x<-rbinom(n=1,size=13,prob=0.40)
  y<-rbinom(n=1,size=30,prob=0.40)
  r[i]<-twostage.inference(x+y,3,13,43,0.20,0.05)[1]
}



for (i in 1:1000){
  x<-rbinom(n=1,size=18,prob=0.40)
  y<-rbinom(n=1,size=15,prob=0.40)
  c[i]<-if(x<4){
    c[i]<-x/18} else{
      c[i]<-(x+y)/33
    }
}

hist(c)
hist(z)

par(mfrow=c(1,3))
hist(c,main='Distribution of MLE',xlab = 'Response Rate estimate')
hist(z,main='Distribution of UMVUE (Optimal)',xlab = 'Response Rate estimate')
hist(r,main='Distribution of UMVUE (Minimax)',xlab = 'Response Rate estimate')



#Distribution of MLE and UMVUE null p=0.2 alternative p=0.4 actual p=0.3

for (i in 1:1000){
  x<-rbinom(n=1,size=13,prob=0.30)
  y<-rbinom(n=1,size=30,prob=0.30)
  z[i]<-twostage.inference(x+y,3,13,43,0.20,0.05)[1]
}


for (i in 1:1000){
  x<-rbinom(n=1,size=13,prob=0.30)
  y<-rbinom(n=1,size=30,prob=0.30)
  r[i]<-twostage.inference(x+y,3,13,43,0.20,0.05)[1]
}



for (i in 1:1000){
  x<-rbinom(n=1,size=13,prob=0.30)
  y<-rbinom(n=1,size=30,prob=0.30)
  c[i]<-if(x<4){
    c[i]<-x/18} else{
      c[i]<-(x+y)/33
    }
}



par(mfrow=c(1,3))
hist(c,main='Distribution of MLE',xlab = 'Response Rate estimate')
hist(z,main='Distribution of UMVUE (Optimal)',xlab = 'Response Rate estimate')
hist(r,main='Distribution of UMVUE (Minimax)',xlab = 'Response Rate estimate')


#Distribution of MLE and UMVUE null p=0.2 alternative p=0.4 actual p=0.2

for (i in 1:1000){
  x<-rbinom(n=1,size=13,prob=0.20)
  y<-rbinom(n=1,size=30,prob=0.20)
  z[i]<-twostage.inference(x+y,3,13,43,0.20,0.05)[1]
}


for (i in 1:1000){
  x<-rbinom(n=1,size=13,prob=0.20)
  y<-rbinom(n=1,size=30,prob=0.20)
  r[i]<-twostage.inference(x+y,3,13,43,0.20,0.05)[1]
}



for (i in 1:1000){
  x<-rbinom(n=1,size=13,prob=0.20)
  y<-rbinom(n=1,size=30,prob=0.20)
  c[i]<-if(x<4){
    c[i]<-x/18} else{
      c[i]<-(x+y)/33
    }
}

hist(c)
hist(z)

par(mfrow=c(1,3))
hist(c,main='Distribution of MLE',xlab = 'Response Rate estimate')
hist(z,main='Distribution of UMVUE (Optimal)',xlab = 'Response Rate estimate')
hist(r,main='Distribution of UMVUE (Minimax)',xlab = 'Response Rate estimate')


#Distribution of MLE and UMVUE null p=0.2 alternative p=0.4 actual p=0.1

for (i in 1:1000){
  x<-rbinom(n=1,size=13,prob=0.10)
  y<-rbinom(n=1,size=30,prob=0.10)
  z[i]<-twostage.inference(x+y,3,13,43,0.20,0.05)[1]
}


for (i in 1:1000){
  x<-rbinom(n=1,size=13,prob=0.10)
  y<-rbinom(n=1,size=30,prob=0.10)
  r[i]<-twostage.inference(x+y,3,13,43,0.20,0.05)[1]
}



for (i in 1:1000){
  x<-rbinom(n=1,size=13,prob=0.10)
  y<-rbinom(n=1,size=30,prob=0.10)
  c[i]<-if(x<4){
    c[i]<-x/18} else{
      c[i]<-(x+y)/33
    }
}

hist(c)
hist(z)

par(mfrow=c(1,3))
hist(c,main='Distribution of MLE',xlab = 'Response Rate estimate')
hist(z,main='Distribution of UMVUE (Optimal)',xlab = 'Response Rate estimate')
hist(r,main='Distribution of UMVUE (Minimax)',xlab = 'Response Rate estimate')

#Power in Non-inferiority,Equivalence and Superiority trials


#Mean difference is 20

pwr_sup<-rep(0,2000)
pwr_ninf<-rep(0,2000)
  pwr_equiv<-rep(0,2000)
i<-0
q<-0
r<-0

for (n in 0:2000){
  pwr_sup[n]<-1-pnorm(1.64-(20*sqrt(n))/sqrt(7200))
  i<-i+1
}

index<-1:2000


for (n in 0:2000){
  pwr_ninf[n]<-pnorm(-1.64+(10*sqrt(n))/sqrt(7200))
  q<-q+1
}

for (n in 0:2000){
  pwr_equiv[n]<-pnorm(-1.64 +(10*sqrt(n))/sqrt(7200))-pnorm(1.64 -(10*sqrt(n))/sqrt(7200))
}

pnorm(3)

pwr_equiv

par(mfrow=c(1,1))
plot(index,pwr_sup,type="l",col="red",xlab = 'sample size',ylab = 'power')
lines(index,pwr_ninf,type="l",lty=2)
lines(index,pwr_equiv,col="blue",type="l",lty=3)
legend(1500,0.85, legend=c("Superiority", "Non-Inferiority","Equivalency"),
       col=c("red", "black","blue"),lty=1:3,  cex=0.8)

#Mean difference is 10
pwr_sup<-rep(0,2000)
pwr_ninf<-rep(0,2000)
pwr_equiv<-rep(0,2000)
i<-0
q<-0
r<-0

for (n in 0:2000){
  pwr_sup[n]<-pnorm(-1.64+(10*sqrt(n))/sqrt(7200))
}

index<-1:2000


for (n in 0:2000){
  pwr_ninf[n]<-pnorm(-1.64+(5*sqrt(n))/sqrt(7200))
}

for (n in 0:2000){
  pwr_equiv[n]<-pnorm(-1.64 +(5*sqrt(n))/sqrt(7200))-pnorm(1.64-(5*sqrt(n))/sqrt(7200))
}


plot(index,pwr_sup,type="l",col="red",xlab = 'sample size',ylab = "power")
lines(index,pwr_ninf,lty=2)
lines(index,pwr_equiv,col="blue",lty=3)
legend(1500,0.30, legend=c("Superiority", "Non-Inferiority","Equivalency"),
       col=c("red", "black","blue"),lty=1:3,  cex=0.8)



#Survival analysis
pwr_surv<-rep(0,2000)
r<-0


for (d in 0:2000){
pwr_surv[d]<-pnorm(-1.96-0.5*log(0.75)*sqrt(d))
}

pwr_surv


index<-1:2000

par(mfrow=c(1,1))
plot(index,pwr_surv,type="l",col="red",xlab = 'Number of Events',ylab = "power")



#Group Sequential designs

install.packages("gsDesign")
library(gsDesign)
library(ggpubr)

a<-gsDesign(k=5,test.type=2,alpha=0.025,beta=0.2,
            timing = c(0.2,0.4,0.6,0.8,1),sfu="Pocock",n.fix=285)
b<-gsDesign(k=5,test.type=2,alpha=0.025,beta=0.2,timing = c(0.2,0.4,0.6,0.8,1),sfu="OF",n.fix=285)

a$upper

boundary_pocock<-rep(2.413176,5)
boundary_fleming<-c(4.561743, 3.225639, 2.633723, 2.280871, 2.040073)


index<-1:5

plot(index,boundary_pocock,type="b",lty=1,xlab = 'stage',ylab = "threshold",ylim=c(-5,5))
lines(index,boundary_fleming,type="b",lty=2)
lines(index,-boundary_pocock,type="b",lty=1)
lines(index,-boundary_fleming,type="b",lty=2)
legend(3.2,5, legend=c("Pocock Boundary","O' Brien and Fleming Boundary"),
       lty=1:2,  cex=0.8)




c<-gsDesign(k=2,test.type=2,alpha=0.025,beta=0.2,timing = c(0.5,1),sfu="Pocock",n.fix=285)
d<-gsDesign(k=2,test.type=2,alpha=0.025,beta=0.2,timing = c(0.5,1),sfu="OF",n.fix=285)

c$upper

e<-gsDesign(k=8,test.type=2,alpha=0.025,beta=0.2,timing = c(0.125,0.25,0.375,0.500,0.625,0.75,0.875,1),sfu="Pocock",n.fix=285)
f<-gsDesign(k=8,test.type=2,alpha=0.025,beta=0.2,timing = c(0.125,0.25,0.375,0.500,0.625,0.75,0.875,1),sfu="OF",n.fix=285)


e$upper

g<-gsDesign(k=10,test.type=2,alpha=0.025,beta=0.2,timing =c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),sfu="Pocock",n.fix=285)
h<-gsDesign(k=10,test.type=2,alpha=0.025,beta=0.2,timing = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),sfu="OF",n.fix=285)


a
b
c
d
e
f
g
h

par

inflation_factor_ob_stage_0.025<-c(1.008,1.028,1.037,1.040)

inflation_factor_po_stage_0.025<-c(1.110,1.229,1.279,1.301)


index<-1:4

par(mfrow=c(1,1))
plot(index,inflation_factor_ob_stage,type="b",lty=1,xlab = 'Number of Stages',ylab = "",ylim=c(0,2))
lines(index,inflation_factor_po_stage,type="b",lty=2)
legend(3.2,5, legend=c("O' Brien Fleming Inflation Factor","Pocock Inflation Factor"),
       lty=1:2,  cex=0.8)
       
########################################################### DESIGNING A SURVIVAL CLINICAL TRIAL #################################################

install.packages("asaur")
library(asaur) 

install.packages("survival")
library(survival)

install.packages("survminer")
library(survminer) 

install.packages("KMsurv")
library(KMsurv)

library(dplyr)


data(larynx)

str(larynx)

ggsurvplot(survfit(Surv(time,delta) ~1,larynx),conf.int = F,palette="black",legend.title="")

larynx$stage<-as.factor(larynx$stage)



cox_mod_larynx_2<-coxph(Surv(time,delta)~stage+age,data=larynx)

cox_mod_larynx_2

logLik(cox_mod_larynx_2)


par(mfrow=c(1,1))
#SCHOENFELD
test.ph = cox.zph(cox_mod_larynx_2,terms=FALSE, singledf=FALSE, global=TRUE)
test.ph



ggcoxzph(test.ph)


cox_mod_larynx<-coxph(Surv(time,delta)~as.factor(stage),data=larynx)

logLik(cox_mod_larynx)

-2*(-logLik(cox_mod_larynx_2)+logLik(cox_mod_larynx))

pchisq(1.826852, df=1, lower.tail=F)




#SCHOENFELD
test.ph =  cox.zph(cox_mod_larynx,terms=FALSE, singledf=FALSE, global=TRUE)
test.ph

ggcoxzph(test.ph)




#COX-SNELL RES
larynx$res3<- residuals(cox_mod_larynx, type="martingale",data=larynx) #see martingale section in [50]
larynx$csres<-larynx$delta-larynx$res3
surv3<-survfit(Surv(csres,delta)~1,type="fleming-harrington",data=larynx)
summary(surv3)

plot(surv3$time,-log(surv3$surv),ylab="Cumulative Hazard of CS-Residuals",xlab="CS-Residuals")
abline(a=0,b=1)

ggsurvplot(survfit(Surv(time,delta) ~as.factor(stage),larynx))

survival<-survfit(Surv(time,delta) ~as.factor(stage),larynx)

summary(survival)
       
