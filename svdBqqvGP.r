source("bqqv.r")
buildBasis <- function(response,percent=.95,numbas=NULL)
{
    svdm <- svd(response)
    if(is.null(numbas))
    {
        cumd <- cumsum(svdm$d)/sum(svdm$d)
        numbas <- min(which(cumd>=percent))
    }
    basis <- svdm$u[,1:numbas,drop=FALSE]          #basis matrix TT by numbas
    redd <- svdm$d[1:numbas]            #reduced d vector length numbas
    redv <- svdm$v[,1:numbas,drop=FALSE]           #reduced v matrix nn by mumbas
    basis <- t(t(basis)*redd)            #coefficient matrix nn by numbas
    coeff <- redv
    if(!is.matrix(coeff))
        coeff <- matrix(coeff,ncol=numbas)
    ret <- list(basis=basis,redd=redd,coeff=coeff,
                numbas=numbas)
}
svdBqqvGP_worker <- function(xx,zz,response,lbasis,nn,TT,covtype,nthread=2)
{
    basis <- lbasis$basis
    numbas <- ncol(basis)
    coeff <- lbasis$coeff
    resid <- response-basis%*%t(coeff)
    varres <- var(as.vector(resid))
    mdlist <- vector(mode="list",length=numbas)
    nthread <- min(nthread,numbas)
    estim <- get(paste("estloglik",covtype,sep=""))
    cl <- parallel::makeCluster(nthread)
    parallel::clusterEvalQ(cl,{source("bqqv.r");
        dyn.load("bqqv.so")})
    mdlist <- tryCatch(parallel::parApply(cl,coeff,2,estim,xx,zz),
                       finally=parallel::stopCluster(cl))
    ret <- list(mdlist=mdlist,basis=basis,numbas=numbas,TT=TT,varres=varres)
    class(ret) <- "svdBqqvGP"
    return(ret)
}
## svdBqqvGPst_worker <- function(xx,zz,response,lbasis,nn,TT,covtype)
## {
##     basis <- lbasis$basis
##     numbas <- ncol(basis)
##     coeff <- lbasis$coeff
##     resid <- response-basis%*%t(coeff)
##     varres <- var(as.vector(resid))
##     estim <- get(paste("estloglik",covtype,sep=""))
##     mdlist <- apply(coeff,2,estim,xx,zz)
##     ret <- list(mdlist=mdlist,basis=basis,numbas=numbas,TT=TT,varres=varres)
##     class(ret) <- "svdBqqvGP"
##     return(ret)
## }

svdBqqvGP <- function(xx,zz,response,covtype=c("Prod","Add","AddHom"),percent=.95,numbas=NULL,nthread=2)
{
    covtype=match.arg(covtype)
    xx <- as.matrix(xx)
    response <- as.matrix(response)
    nn <- nrow(xx)
    if(ncol(response)!=nn)
        stop("inconsistent number of design points!")
    TT <- nrow(response)
    lbasis <- buildBasis(response,percent,numbas)
    ret <- svdBqqvGP_worker(xx,zz,response,lbasis,nn,TT,covtype,nthread)
    return(ret)
}
predict.svdBqqvGP <- function(model,newxx,newzz,...)
{
    newxx <- as.matrix(newxx)
    numbas <- model$numbas
    nt <- nrow(newxx)
    coeff <- matrix(nrow=nt,ncol=numbas)
    varmat <- coeff
    for(i in 1:numbas)
    {
        md <- model$mdlist[[i]]
        pred <- predict(md,newxx,newzz)
        coeff[,i] <- pred$pmean
        varmat[,i] <- pred$psig2
    }
    ave <- model$basis%*%t(coeff)
    sd <- sqrt(model$basis^2%*%t(varmat)+model$varres)
    varsvd <- sum(model$basis^2%*%t(varmat))
    res <- list(mean=ave,sd=sd,coefvar=varmat,coefmean=coeff,varsvd=varsvd)
    return(res)
}
