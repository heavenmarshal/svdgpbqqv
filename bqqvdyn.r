#dyn.load("bqqvdyn.so")
findlevel <- function(x)
{
    return(length(unique(x)))
}
boundAddDyn <- function(d,q,resp,levels)
{
    lentht <- sum(levels*(levels-1)/2)
    lb <- c(rep(1e-10,d*q+q+1),rep(0,lentht))
    ub <- c(rep(10,d*q),.5,rep(max(apply(resp,1,var)),q),rep(pi,lentht))
    return(list(ub=ub,lb=lb,lentht=lentht))
}

loglikAddDyn <- function(param,levels,X,Z,tv,resp,n,L,d,q)
{
    phis <- head(param,d*q)
    phit <- param[(d*q+1)]
    vsigma2 <- param[(d*q+1):(d*q+1+q)]
    vtheta <- param[-(1:(d*q + q + 1))]
    out <- .C("loglikAddDyn_R", as.double(phis), as.double(phit), as.double(vtheta),
              as.double(vsigma2), as.integer(levels),as.double(X), as.integer(Z),
              as.double(tv), as.double(resp), as.integer(n), as.integer(L), as.integer(d),
              as.integer(q), ans = double(1))
    return(out$ans)
}

estloglikAddDyn <- function(X,Z,resp,ninit=200)
{
    n <- nrow(X)
    L <- nrow(resp)
    d <- ncol(X)
    q <- ncol(Z)
    bk1 <- q*d
    bk2 <- 1+bk1
    bk3 <- q+bk2

    tv <- seq(0,1,length=L)
    levels <- apply(Z,2,findlevel)
    bound <- boundAddDyn(d,q,resp,levels)
    nvar <- bk3+bound$lentht
    ## initmat <- t(lhs::maximinLHS(ninit,d*q+q+1+bound$lentht))*(bound$ub-bound$lb)+bound$lb
    ## initlik <- apply(initmat,2,loglikAddDyn,levels,X,Z,tv,resp,n,L,d,q)
    ## start <- initmat[,which.max(initlik)]
    ## opt <- optim(start,loglikAddDyn,NULL,levels,X,Z,tv,resp,n,L,d,q,
    ##              method="L-BFGS-B",lower=bound$lb, upper=bound$ub,
    ##              control=list(fnscale=-1))
    likwrap <- function(param)loglikAddDyn(param,levels,X,Z,tv,resp,n,L,d,q)
    opt <- rgenoud::genoud(likwrap,nvar,max=TRUE,
                           Domains=cbind(bound$lb,bound$ub),boundary.enforcement=2,
                           print.level=0)
    param <- opt$par
    phis <- matrix(param[1:bk1],nrow=d)
    phit <- param[(bk1+1):bk2]
    vsigma2 <- param[(bk2+1):bk3]
    vtheta <- param[-(1:bk3)]
    ret <- list(#stat=opt$convergence,
        phis=phis,phit=phit,vsigma2=vsigma2,
        vtheta=vtheta, levels=levels, X=X, Z=Z, tv = tv,resp=resp,
        n=n,L=L,d=d,q=q)
    class(ret) <- "bqqvaddDyn"
    return(ret)
}
predict.bqqvaddDyn <- function(md,newX,newZ)
{
    n0 <- nrow(newX)
    out <- .C("addDynPred_R",as.double(newX), as.integer(newZ),
              as.double(md$X), as.integer(md$Z), as.double(md$tv), as.double(md$resp),
              as.double(md$phis), as.double(md$phit), as.double(md$vtheta), as.double(md$vsigma2),
              as.integer(md$levels), as.integer(n0), as.integer(md$n),
              as.integer(md$L),as.integer(md$d), as.integer(md$q),
              pmean=double(n0*md$L))
    pmean <- matrix(out$pmean,nrow=md$L)
    return(list(pmean=pmean))
}    
