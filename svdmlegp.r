buildBasis <- function(response,percent=.95,numbas=NULL)
{
    svdm <- svd(response)
    if(is.null(numbas))
    {
        cumd <- cumsum(svdm$d)/sum(svdm$d)
        numbas <- min(which(cumd>=percent))
    }
    basis <- svdm$u[,1:numbas]          #basis matrix TT by numbas
    redd <- svdm$d[1:numbas]            #reduced d vector length numbas
    redv <- svdm$v[,1:numbas]           #reduced v matrix nn by mumbas
    basis <- t(t(basis)*redd)            #coefficient matrix nn by numbas
    coeff <- redv
    if(!is.matrix(coeff))
        coeff <- matrix(coeff,ncol=numbas)
    ret <- list(basis=basis,redd=redd,coeff=coeff,numbas=numbas)
}
mlegpsvdgp_worker<- function(design,response,lbasis,nn,TT)
{
    basis <- lbasis$basis
    numbas <- ncol(basis)
    coeff <- lbasis$coeff
    resid <- response-basis%*%t(coeff)
    varres <- var(as.vector(resid))
    mdlist <- vector(mode="list",length=numbas)
    for(i in 1:numbas)
    {
        yy <- coeff[,i]
        sink("NUL")
        md <- mlegp::mlegp(design,yy)
        sink()
        mdlist[[i]] <- md
    }
    ret <- list(mdlist=mdlist,basis=basis,numbas=numbas,TT=TT,varres=varres)
    class(ret) <- "svdgp"
    return(ret)
}
kmsvdgp_worker <- function(design,response,lbasis,nn,TT)
{
    basis <- lbasis$basis
    numbas <- ncol(basis)
    coeff <- lbasis$coeff
    resid <- response-basis%*%t(coeff)
    varres <- var(as.vector(resid))
    mdlist <- vector(mode="list",length=numbas)
    for(i in 1:numbas)
    {
        yy <- coeff[,i]
        md <- DiceKriging::km(~1,design,yy,"gauss",nugget.estim=TRUE,control=list(trace=FALSE))
        mdlist[[i]] <- md
    }
    ret <- list(mdlist=mdlist,basis=basis,numbas=numbas,TT=TT,varres=varres)
    class(ret) <- "kmsvdgp"
    return(ret)
}
gpfsvdgp_worker <- function(design,response,lbasis,nn,TT)
{
    basis <- lbasis$basis
    numbas <- ncol(basis)
    coeff <- lbasis$coeff
    resid <- response-basis%*%t(coeff)
    varres <- var(as.vector(resid))
    mdlist <- vector(mode="list",length=numbas)
    for(i in 1:numbas)
    {
        yy <- coeff[,i]
        md <- GPfit::GP_fit(design,yy,corr=list(type="exponential",power=2))
        mdlist[[i]] <- md
    }
    ret <- list(mdlist=mdlist,basis=basis,numbas=numbas,TT=TT,varres=varres)
    class(ret) <- "gpfsvdgp"
}
svdgp <- function(design,response,percent=.95,numbas=NULL,method=c("mlegp","km","gpf"))
{
    method <- match.arg(method)
    design <- as.matrix(design)
    response <- as.matrix(response)
    nn <- nrow(design)
    if(ncol(response)!=nn)
        stop("inconsistent number of design points!")
    TT <- nrow(response)
    fname <- paste(method,"svdgp_worker",sep="")
    worker <- get(fname)
    lbasis <- buildBasis(response,percent,numbas)
    ret <- worker(design,response,lbasis,nn,TT)
    return(ret)
}

predict.svdgp <- function(model,newdata)
{
    newdata <- as.matrix(newdata)
    numbas <- model$numbas
    nt <- nrow(newdata)
    coeff <- matrix(nrow=nt,ncol=numbas)
    varmat <- coeff
    for(i in 1:numbas)
    {
        md <- model$mdlist[[i]]
        pred <- mlegp::predict.gp(md,newdata,se.fit=TRUE)
        coeff[,i] <- pred$fit
        varmat[,i] <- pred$se.fit^2
    }
    ave <- model$basis%*%t(coeff)
    sd <- sqrt(model$basis^2%*%t(varmat)+model$varres)
    ## debug code
    varsvd <- sum(model$basis^2%*%t(varmat))
    ## finish debug
    res <- list(mean=ave,sd=sd,coefvar=varmat,coefmean=coeff,varsvd=varsvd)
    return(res)
}
predict.kmsvdgp <- function(model,newdata)
{
    newdata <- as.matrix(newdata)
    numbas <- model$numbas
    nt <- nrow(newdata)
    coeff <- matrix(nrow=nt,ncol=numbas)
    varmat <- coeff
    for(i in 1:numbas)
    {
        md <- model$mdlist[[i]]
        pred <- DiceKriging::predict.km(md,newdata,type="UK",checkNames=FALSE)
        coeff[,i] <- pred$mean
        varmat[,i] <- pred$sd^2
    }
    ave <- model$basis%*%t(coeff)
    sd <- sqrt(model$basis^2%*%t(varmat)+model$varres)
    ## debug code
    varsvd <- sum(model$basis^2%*%t(varmat))
    ## finish debug
    res <- list(mean=ave,sd=sd,coefvar=varmat,coefmean=coeff,varsvd=varsvd)
    return(res)
}
predict.gpfsvdgp <- function(model,newdata)
{
    newdata <- as.matrix(newdata)
    numbas <- model$numbas
    nt <- nrow(newdata)
    coeff <- matrix(nrow=nt,ncol=numbas)
    varmat <- coeff
    for(i in 1:numbas)
    {
        md <- model$mdlist[[i]]
        pred <- GPfit::predict.GP(md,newdata)
        coeff[,i] <- pred$Y_hat
        varmat[,i] <- pred$MSE
    }
    ave <- model$basis%*%t(coeff)
    sd <- sqrt(model$basis^2%*%t(varmat)+model$varres)
    ## debug code
    varsvd <- sum(model$basis^2%*%t(varmat))
    ## finish debug
    res <- list(mean=ave,sd=sd,coefvar=varmat,coefmean=coeff,varsvd=varsvd)
    return(res)
}
