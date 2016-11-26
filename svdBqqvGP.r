source("bqqv.r")
buildBasis <- function(response,percent,fracns)
{
    svdm <- svd(response)
    cumd <- cumsum(svdm$d)/sum(svdm$d)
    numbas <- min(which(cumd>=percent))
    numns <- min(which(cumd>=fracns))
    
    basis <- svdm$u[,1:numbas,drop=FALSE]          #basis matrix TT by numbas
    redd <- svdm$d[1:numbas]            #reduced d vector length numbas
    redv <- svdm$v[,1:numbas,drop=FALSE]           #reduced v matrix nn by mumbas
    basis <- t(t(basis)*redd)            #coefficient matrix nn by numbas
    coeff <- redv
    if(!is.matrix(coeff))
        coeff <- matrix(coeff,ncol=numbas)
    ret <- list(basis=basis,redd=redd,coeff=coeff,
                numbas=numbas,numns=numns)
}
svdBqqvGP_worker <- function(xx,zz,coeff,numbas,covtype,nthread)
{
    mdlist <- vector(mode="list",length=numbas)
    nthread <- min(nthread,numbas)
    estim <- get(paste("estloglik",covtype,sep=""))
    cl <- parallel::makeCluster(nthread)
    parallel::clusterEvalQ(cl,{source("bqqv.r");
        dyn.load("bqqv.so")})
    mdlist <- tryCatch(parallel::parApply(cl,coeff,2,estim,xx,zz),
                       finally=parallel::stopCluster(cl))
    return(mdlist)
}
svdgp_worker <- function(xx,coeff,numbas,nthread)
{
    nthread <- min(nthread,numbas)
    cl <- parallel::makeCluster(nthread)
    parallel::clusterEvalQ(cl,{library("mlegp")})
    parallel::clusterExport(cl,"xx",envir=environment())
    mdlist <- tryCatch(parallel::parApply(cl,coeff,2,
                                          function(y)mlegp::mlegp(xx,y,verbose=0)),
                       finally=parallel::stopCluster(cl)
                       )
    return(mdlist)
}
svdBqqvGP <- function(xx,zz,response,covtype=c("Prod","Add"),percent=.95,
                      fracns=.5,nthread=2)
{
    covtype=match.arg(covtype)
    xx <- as.matrix(xx)
    response <- as.matrix(response)
    nn <- nrow(xx)
    if(ncol(response)!=nn)
        stop("inconsistent number of design points!")
    TT <- nrow(response)
    lbasis <- buildBasis(response,percent,fracns)
    resid <- response-lbasis$basis%*%t(lbasis$coeff)
    varres <- var(as.vector(resid))
    numbas <- lbasis$numbas
    numns <- lbasis$numns
    coeffns <- lbasis$coeff[,1:numns]
    mdlistns <- svdBqqvGP_worker(xx,zz,coeffns,numns,covtype,nthread)
    mdlistst <- NULL
    numst <- numbas-numns
    if(numst>0)
    {

        coeffst <- lbasis$coeff[,-(1:numns)]
        mdlistst <- svdgp_worker(xx,coeffst,numst,nthread)
    }
    ret <- list(mdlistns=mdlistns,mdlistst=mdlistst,basis=lbasis$basis,
                numbas=numbas,numns=numns,numst=numst,TT=TT,varres=varres)
    class(ret) <- "svdBqqvGP"
    return(ret)
}
predict.svdBqqvGP <- function(model,newxx,newzz,...)
{
    newxx <- as.matrix(newxx)
    numbas <- model$numbas
    numns <- model$numns
    numst <- model$numst
    nt <- nrow(newxx)
    coeffns <- matrix(nrow=nt,ncol=numns)
    varmatns <- coeffns
    for(i in 1:numns)
    {
        md <- model$mdlistns[[i]]
        pred <- predict(md,newxx,newzz)
        coeffns[,i] <- pred$pmean
        varmatns[,i] <- pred$psig2
    }
    coeffst <- varmatst <- NULL
    if(numst>0)
    {
        coeffst <- matrix(nrow=nt,ncol=numst)
        varmatst <- coeffst
        for(i in 1:numst)
        {
            md <- model$mdlistst[[i]]
            pred <- mlegp::predict.gp(md,newxx,se.fit=TRUE)
            coeffst[,i] <- pred$fit
            varmatst[,i] <- (pred$se.fit)^2
        }
    }
    coeff <- cbind(coeffns,coeffst)
    varmat <- cbind(varmatns,varmatst)
    ave <- model$basis%*%t(coeff)
    sd <- sqrt(model$basis^2%*%t(varmat)+model$varres)
    varsvd <- sum(model$basis^2%*%t(varmat))
    res <- list(mean=ave,sd=sd,coefvar=varmat,coefmean=coeff,varsvd=varsvd)
    return(res)
}
