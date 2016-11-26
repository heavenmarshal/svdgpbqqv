#dyn.load("/home/hpc3104/work/Rscripts/Rinclude/bqqv.so")
library("lhs")
findlevel <- function(x)
{
    return(length(unique(x)))
}
boundProd <- function(d,levels)
{
    lentht <- sum(levels*(levels-1)/2)
    lb <- c(rep(1e-10,d),rep(0,lentht))
    ub <- c(rep(10,d),rep(pi,lentht))
    return(list(lb=lb,ub=ub,lentht=lentht))
}
boundAdd <- function(d,q,resp,levels)
{
    lentht <- sum(levels*(levels-1)/2)
    lb <- c(rep(1e-10,d*q+q),rep(0,lentht))
    ub <- c(rep(10,d*q),rep(var(resp),q),rep(pi,lentht))
    return(list(ub=ub,lb=lb,lentht=lentht))
}
genInit <- function(bound,control,covtype,levels,X,Z,resp,n,d,q)
{
    psearch <- control[1]
    ppercent <- control[2]
    pclust <- control[3]
    nparam <- length(bound$ub)
    likfun <- get(paste("loglik",covtype,sep=""))
    initmat1 <- t(lhs::maximinLHS(psearch,nparam))*(bound$ub-bound$lb)+bound$lb
    initlik1 <- apply(initmat1,2,likfun,levels,X,Z,resp,n,d,q)
    initidx1 <- order(initlik1,decreasing=TRUE)[1:ppercent]
    
    initmat2 <- t(initmat1[,initidx1])
    initlik2 <- initlik1[initidx1]
    kk <- replicate(5,kmeans(initmat2,pclust),simplify=FALSE)
    wsssum <- sapply(sapply(kk,"[","withinss"),sum)
    kk <- kk[[which.max(wsssum)]]       #should it be the smallest?
    initidx2 <- initlik2%in%tapply(initlik2,kk$cluster,max)
    initmat <- initmat2[initidx2,]
}
loglikProd <- function(param,levels,X,Z,resp,n,d,q)
{
    phi <-  head(param,d)
    vtheta <-  param[-(1:d)]
    out <- .C("loglikProd_R",as.double(phi), as.double(vtheta),
              as.integer(levels), as.double(X), as.integer(Z),
              as.double(resp), as.integer(n),as.integer(d),
              as.integer(q), ans = double(1))
    return(out$ans)
}

loglikAdd <- function(param,levels,X,Z,resp,n,d,q)
{
    phi <- head(param,d*q)
    vsigma2 <- param[(d*q+1):(d*q+q)]
    vtheta <- param[-(1:(d*q+q))]
    out <- .C("loglikAdd_R", as.double(phi), as.double(vtheta),
              as.double(vsigma2), as.integer(levels),as.double(X), as.integer(Z),
              as.double(resp), as.integer(n), as.integer(d),
              as.integer(q), ans = double(1))
    return(out$ans)
}
estloglikProd <- function(resp,X,Z,control=c(200*pvar,80*pvar,2*pvar))
{
    n <- nrow(X)
    d <- ncol(X)
    q <- ncol(Z)
    pvar <- d+q
    levels <- apply(Z,2,findlevel)
    bound <- boundProd(d,levels)
    matinit <- genInit(bound,control,"Prod",levels,X,Z,resp,n,d,q)
    ninit <- nrow(matinit)
    maxlik <- -Inf
    for(i in 1:ninit)
    {
        topt <- optim(matinit[i,],loglikProd,NULL,levels,X,Z,resp,n,d,q,
                     method="L-BFGS-B",lower=bound$lb, upper=bound$ub,
                     control=list(fnscale=-1))
        if(topt$value>maxlik)
        {
            maxlik <- topt$value
            opt <- topt
        }
    }
    param <- opt$par
    phi <- param[1:d]
    vtheta <- param[-(1:d)]
    ret <- list(converge=opt$convergence,phi=phi,vtheta=vtheta,
                levels=levels,X=X,Z=Z,resp=resp,n=n,d=d,q=q)
    class(ret) <- "bqqvprod"
    return(ret)
}
predict.bqqvprod <- function(md,newX,newZ)
{
    n0 <- nrow(newX)
    out <- .C("prodPredict_R", as.double(newX), as.integer(newZ),
              as.double(md$X), as.integer(md$Z), as.double(md$resp),
              as.double(md$phi), as.double(md$vtheta), as.integer(md$levels),
              as.integer(n0), as.integer(md$n), as.integer(md$d),
              as.integer(md$q), pmean=double(n0), psig2=double(n0))
    ret <- list(pmean=out$pmean,psig2=out$psig2)
    return(ret)
}
estloglikAdd <- function(resp,X,Z,control=c(200*pvar,80*pvar,2*pvar))
{
    n <- nrow(X)
    d <- ncol(X)
    q <- ncol(Z)
    pvar <- d+q
    bk1 <- q*d
    bk2 <- q+bk1
    levels <- apply(Z,2,findlevel)
    bound <- boundAdd(d,q,resp,levels)
    matinit <- genInit(bound,control,"Add",levels,X,Z,resp,n,d,q)
    ninit <- nrow(matinit)
    maxlik <- -Inf
    for(i in 1:ninit)
    {
        topt <- optim(matinit[i,],loglikAdd,NULL,levels,X,Z,resp,n,d,q,
                      method="L-BFGS-B",lower=bound$lb, upper=bound$ub,
                      control=list(fnscale=-1))
        if(topt$value>maxlik)
        {
            maxlik <- topt$value
            opt <- topt
        }
    }
    param <- opt$par
    phi <- matrix(param[1:bk1],nrow=d)
    vsigma2 <- param[(bk1+1):bk2]
    vtheta <- param[-(1:bk2)]
    ret <- list(converge=opt$convergence,phi=phi,vsigma2=vsigma2,vtheta=vtheta,
                levels=levels,X=X,Z=Z,resp=resp,n=n,d=d,q=q)
    class(ret) <- "bqqvadd"
    return(ret)
}
predict.bqqvadd <- function(md,newX,newZ)
{
    n0 <- nrow(newX)
    out <- .C("addPredict_R", as.double(newX), as.integer(newZ),
              as.double(md$X), as.integer(md$Z), as.double(md$resp),
              as.double(md$phi), as.double(md$vtheta), as.double(md$vsigma2),
              as.integer(md$levels), as.integer(n0), as.integer(md$n),
              as.integer(md$d), as.integer(md$q),
              pmean=double(n0), psig2=double(n0))
    ret <- list(pmean=out$pmean,psig2=out$psig2)
    return(ret)
}
