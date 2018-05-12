#===========================================
## @knitr code_part0
#  TITLE:  computational economics: pset 6
#  AUTHOR: John Eric Humphries
#  Replicate part of Bresnahan and Reiss 1991
#
#  Date: May-3-2014
#================================


#========================
## @knitr code_part0
# Section 0: setup
#========================

rm(list=ls())           # Clear the workspace
set.seed(907) 
library(ggplot2)
library(numDeriv)
library(nloptr)
library(Rsolnp)

# Readiing the data
br <- read.csv("BresnahanAndReiss1991_DATA.csv")

#View(BresnahanAndReiss1991_DATA)


#==================
## @knitr code_part1
#  Section 1: setup
#===================

V <- function (br,n, alpha=rep(1,5) ,bet=rep(1,4) ) {
    profit <- alpha[1] + bet[1]*br$ELD + bet[2]*br$PINC + bet[3]*br$LNHDD + bet[4]*br$FFRAC    
    if (n>=2){
        for (i in 2:n)
            profit <- profit - alpha[i]
    }
    return(profit)
}

F <- function (br,n, gam=rep(1,6) ) {
    F<- gam[1] + gam[6] *br$LANDV    
    if (n>=2){
        for (i in 2:n)
            F <- F  + gam[i]
    }
    return(F)
} 

# theta values from their paper
theta.true <- c(-.53, 2.25,.34,.23,-.49,-.03,.004,-.02,.86, .03, .15, 0, .08, .53, .76, .46, .60, .12, -.74)

# Liklihood function
lLik <- function(theta=rep(1,19),br=br) {
    lam   <- theta[1:4]
    bet   <- theta[5:8]
    alpha <- theta[9:13]
    gam   <- theta[14:19]
    
    S  <-  br$TPOP + lam[1]*br$OPOP + lam[2]*br$NGRW + lam[3]*br$PGRW + lam[4]*br$OCTY
    P <- list()
    P[[1]] <- log(1 - pnorm(S*V(br,1,alpha=alpha,bet=bet) - F(br,1,gam=gam)))
    P[[2]] <- log(  pnorm(S*V(br,1,alpha=alpha,bet=bet) - F(br,1,gam=gam)) -  pnorm(S*V(br,2,alpha=alpha,bet=bet) - F(br,2,gam=gam))  )
    P[[3]] <- log(  pnorm(S*V(br,2,alpha=alpha,bet=bet) - F(br,2,gam=gam)) -  pnorm(S*V(br,3,alpha=alpha,bet=bet) - F(br,3,gam=gam))  )
    P[[4]] <- log(  pnorm(S*V(br,3,alpha=alpha,bet=bet) - F(br,3,gam=gam)) -  pnorm(S*V(br,4,alpha=alpha,bet=bet) - F(br,4,gam=gam))  )
    P[[5]] <- log(  pnorm(S*V(br,4,alpha=alpha,bet=bet) - F(br,4,gam=gam)) -  pnorm(S*V(br,5,alpha=alpha,bet=bet) - F(br,5,gam=gam))  )
    P[[6]] <- log(  pnorm(S*V(br,5,alpha=alpha,bet=bet) - F(br,5,gam=gam)) )
    for (i in 1:6)
        P[[i]][P[[i]]==-Inf] <- -100000000
    
    lLik <-  sum(P[[1]][br$TIRE==0]) + sum(P[[2]][br$TIRE==1]) + sum(P[[3]][br$TIRE==2]) +
                 sum(P[[4]][br$TIRE==3]) + sum(P[[5]][br$TIRE==4]) + sum(P[[6]][br$TIRE>=5])
    return(-1* lLik)
}


lLik(theta=theta.true, br=br)

#==================
## @knitr code_part2
#  Section 2: numerical optimization
#===================

# What I want them to do
out1 <- optim(par=rep(.1,19), lLik, br=br, lower = c(rep(-Inf,8),rep(0,10),-Inf), upper = rep(Inf,19),method="L-BFGS-B", control = list(maxit=500))

#==================
## @knitr code_part3
#  Section 3: using numerical gradient
#===================

# Defining a function that returns the numeric gradient.
lLik.grad <-  function(theta=rep(1,19),br=br) {
    grad(x=theta,lLik, br= br)
}

out2 <- optim(par=rep(.1,19), lLik, gr =lLik.grad, br=br, lower = c(rep(-Inf,8),rep(0,10),-Inf), upper = rep(Inf,19),method="L-BFGS-B")

# Calculating Standard errors using hessian (not  a great method ofcourse since we are inverting a numerical hessian etc etc.)
out.grad= grad(lLik, x=out2$par, br= br)
out.hessian = hessian(lLik, x=out2$par, br=br)
out.se <- sqrt(diag(solve(out.hessian)))

cbind(out2$par,out1$par, theta.true, out2$par - theta.true, out1$par - theta.true )

results <- cbind(out2$par,out.se)
colnames(results) <- c('estimate','se')
rownames(results) <- c('lambda1','lambda2','lambda3','lambda4','beta1','beta2','beta3','beta4','alpha1','alpha2','alpha3','alpha4','alpha5','gamma1','gamma2','gamma3','gamma4','gamma5','gamma6')
stargazer(results,out='n6_paperReplication/results.tex')

#==================
## @knitr code_part4
#  Section 4: Calculate thresholds
#===================

calc.S <- function(br,theta){
    lam   <- theta[1:4]
    bet   <- theta[5:8]
    alpha <- theta[9:13]
    gam   <- theta[14:19]
    br.m  <- data.frame(t(colMeans(br)))
    
    S <- rep(0,5)
    for(i in 1:5)
        S[i] <- F(br.m,i,gam)/V(br.m,i,alpha,bet)
    
    return(S)
}

S <- calc.S(br,out$par)
SN.S5ratio <- (S[5]*(1:5))/(S*5)
qplot(1:5,SN.S5ratio)
ggsave(file="SNS5ratio.pdf")

#==================
## @knitr code_part5
#  Section 4: Additional optimization methods
#===================

# out3 <- nloptr(x0 = rep(.1,19), lLik, lLik.grad, br=br, lb = c(rep(-Inf,8),rep(0,10),-Inf), ub = rep(Inf,19),
#                opts= list("algorithm"="NLOPT_LD_LBFGS", "maxtime"=500,"maxeval"=500, "print_level" = 1))
# out4 <- nloptr(x0 = rep(.1,19), lLik, lLik.grad, br=br, lb = c(rep(-Inf,8),rep(0,10),-Inf), ub = rep(Inf,19),
#                opts= list("algorithm"="NLOPT_LD_VAR1", "maxtime"=500,"maxeval"=500, "print_level" = 1))
# 
# out5 <- nloptr(x0 = rep(.1,19), lLik, lLik.grad, br=br, lb = c(rep(-Inf,8),rep(0,10),-Inf), ub = rep(Inf,19),
#                opts= list("algorithm"="NLOPT_LD_VAR2", "maxtime"=500,"maxeval"=500, "print_level" = 1))
# 
# local_opts <- list( "algorithm" = "NLOPT_LD_LBFGS",
#                     "xtol_rel"  = 1.0e-10 )
# out6 <- nloptr(x0 = rep(.1,19), lLik, lLik.grad, br=br, lb = c(rep(-Inf,8),rep(0,10),-Inf), ub = rep(Inf,19),
#                opts= list("algorithm"="NLOPT_LD_AUGLAG", "maxtime"=500,"maxeval"=500, "print_level" = 1,"local_opts"=local_opts, "xtol_rel"  = 1.0e-10))
# 
# 
# 




