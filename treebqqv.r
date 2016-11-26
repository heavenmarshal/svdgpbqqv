splitlevels <- function(resp,zz,uniq,minobs=0)
{
    mn <- apply(resp,1,mean)
    minsse <- sum((resp-mn)^2)
    q <- ncol(zz)
    divar <- leftidx <- rightidx <- NULL
    leftlevel <- rightlevel <- NULL
    for(i in 1:q)
    {
        lvec <- uniq[[i]]
        len <- length(lvec)
        if(len<=1) next
        ndic <- floor(len/2)
        if(ndic==len/2)
        {
            lv <- combn(lvec,ndic)
            nlv <- ncol(lv)
            half <- nlv/2
            leftlv <- lv[,1:half,drop=FALSE]
            rightlv <- lv[,nlv:(half+1),drop=FALSE]
            for(k in 1:half)
            {
                lidx <- which(zz[,i]%in%leftlv[,k])
                ridx <- which(zz[,i]%in%rightlv[,k])
                if(length(lidx)<minobs || length(ridx)<minobs)
                    next
                leftresp <- resp[,lidx,drop=FALSE]
                rightresp <- resp[,ridx,drop=FALSE]
                leftmn <- apply(leftresp,1,mean)
                sse <- sum((leftresp-leftmn)^2)
                rightmn <- apply(rightresp,1,mean)
                sse <- sse+sum((rightresp-rightmn)^2)
                if(sse<minsse)
                {
                    leftlevel <- leftlv[,k]
                    rightlevel <- rightlv[,k]
                    minsse <- sse
                    divar <- i
                    leftidx <- lidx
                    rightidx <- ridx
                }
            }
            ndic=ndic-1
        }
        while(ndic>0)
        {
            leftlv <- combn(lvec,ndic)
            half <- ncol(leftlv)
            rightlv <- combn(lvec,len-ndic)
            rightlv <- rightlv[,half:1]
            for(k in 1:half)
            {
                lidx <- which(zz[,i]%in%leftlv[,k])
                ridx <- which(zz[,i]%in%rightlv[,k])
                if(length(lidx)<minobs || length(ridx)<minobs)
                    next
                leftresp <- resp[,lidx,drop=FALSE]
                rightresp <- resp[,ridx,drop=FALSE]
                leftmn <- apply(leftresp,1,mean)
                sse <- sum((leftresp-leftmn)^2)
                rightmn <- apply(rightresp,1,mean)
                sse <- sse+sum((rightresp-rightmn)^2)
                if(sse<minsse)
                {
                    leftlevel <- leftlv[,k]
                    rightlevel <- rightlv[,k]
                    minsse <- sse
                    divar <- i
                    leftidx <- lidx
                    rightidx <- ridx
                }
            }
            ndic=ndic-1
        }
    }
    ret <- if(!is.null(divar)) list(minsse=minsse,divar=divar,
                                    leftidx=leftidx,rightidx=rightidx,
                                    leftlevel=leftlevel,rightlevel=rightlevel)
    return(ret)
}
getlevels <- function(x)
{
    ret <- sort(unique(x))
}
getperiod <- function(x)
{
    ret <- spec.pgram(x,plot=FALSE)$spec
}
buildTree <- function(zz,resp,numNodeBas=1,minNodeObs=10)
{
    tree <- list()
    nn <- nrow(zz)
    pos <- 1
    splitthres <- minNodeObs*2
    clevels <- apply(zz,2,getlevels)
    if(is.matrix(clevels))              #change to list for consistency
        clevels <- lapply(apply(clevels,2,as.list),unlist)
    nlevels <- sapply(clevels,length)
    splitable <- any(nlevels>1)
    modelable <- nlevels>1
    lbasis <- buildBasis(resp,numbas=numNodeBas)
    resid <- resp-lbasis$basis%*%t(lbasis$coeff)
    residperiod <- apply(resid,2,getperiod)
    node <- list(pos=pos,resp=resp,splitable=splitable,
                 modelable=modelable,
                 lbasis=lbasis,resid=resid,residperiod=residperiod,clevels=clevels,
                 nlevels=nlevels,index=1:nn)
    tree[[1]] <- node
    while(pos <= length(tree))
    {
        node <- tree[[pos]]
        if(is.null(node) || !node$splitable)
        {
            pos <- pos+1
            next
        }
        sp <- splitlevels(node$resid,zz[node$index,,drop=FALSE],node$clevels,minNodeObs)
        if(!is.null(sp))
        {
            divar <- sp$divar
            tree[[pos]]$divar <- divar
            lindex <- node$index[sp$leftidx]
            ldlevels <- sp$leftlevel
            lnumobs <- length(lindex)
            lpos <- pos*2
            lclevels <- apply(zz[lindex,,drop=FALSE],2,getlevels)
            if(is.matrix(lclevels))
                lclevels <- lapply(apply(lclevels,2,as.list),unlist)
            lnlevels <- sapply(lclevels,length)
            lsplitable <- any(lnlevels>1) && lnumobs>=splitthres
            lmodelable <- node$modelable & (lnlevels==node$nlevels) #not degenerated because of split
            lmodelable[divar] <- node$modelable[divar] && lnlevels[divar]>1 #variable for split is different
            lresp <- node$resid[,sp$leftidx]
            llbasis <- buildBasis(lresp,numbas=numNodeBas)
            lresid <- lresp-llbasis$basis%*%t(llbasis$coeff)
            lresidperiod <- apply(lresid,2,getperiod)
            lnode <- list(pos=lpos,resp=lresp,splitable=lsplitable,
                          modelable=lmodelable,lbasis=llbasis,
                          resid=lresid,residperiod=lresidperiod,clevels=lclevels,dlevels=ldlevels,
                          nlevels=lnlevels,index=lindex)
            tree[[lpos]] <- lnode

            rindex <- node$index[sp$rightidx]
            rdlevels <- sp$rightlevel
            rpos <- lpos+1
            rnumobs <- length(rindex)
            rclevels <- apply(zz[rindex,,drop=FALSE],2,getlevels)
            if(is.matrix(rclevels))
                rclevels <- lapply(apply(rclevels,2,as.list),unlist)
            rnlevels <- sapply(rclevels,length)
            rsplitable <- any(rnlevels>1) && rnumobs>=splitthres
            rmodelable <- node$modelable & (rnlevels==node$nlevels) #not degenerated because of split
            rmodelable[divar] <- node$modelable[divar] && rnlevels[divar]>1
            rresp <- node$resid[,sp$rightidx]
            rlbasis <- buildBasis(rresp,numbas=numNodeBas)
            rresid <- rresp-rlbasis$basis%*%t(rlbasis$coeff)
            rresidperiod <- apply(rresid,2,getperiod)
            rnode <- list(pos=rpos,resp=rresp,splitable=rsplitable,
                          modelable=rmodelable,lbasis=rlbasis,
                          resid=rresid,residperiod=rresidperiod,clevels=rclevels,dlevels=ldlevels,
                          nlevels=rnlevels,index=rindex)
            tree[[rpos]] <- rnode
        }
        pos <- pos+1
    }
    return(tree)
}
adjlevel <- function(x)
{
    ret <- as.integer(factor(x))-1
}
builder <- function(node,xx,zz,covtype)
{
    xxin <- xx[node$index,]
    nn <- nrow(xxin)
    TT <- nrow(node$resp)
    if(!any(node$modelable))
        return(svdmlegp_worker(xxin,node$resp,node$lbasis,nn,TT))

    zzin <- zz[node$index,node$modelab,drop=FALSE]
    zzin <- apply(zzin,2,adjlevel)
    return(svdBqqvGPst_worker(xxin,zzin,node$resp,node$lbasis,nn,TT,covtype))
}
buildModel <- function(tree,xx,zz,covtype=c("Prod","Add"),nthread=4)
{
    covtype <- match.arg(covtype)
    cl <- parallel::makeCluster(nthread)
    parallel::clusterEvalQ(cl,{dyn.load("bqqv.so");
        source("treebqqv.r");
        source("svdmlegp.r");
        source("svdBqqvGP.r")})
    model <- tryCatch(parallel::parLapply(cl,tree,builder,xx,zz,covtype),
                      finally=parallel::stopCluster(cl))
}
associate <- function(tree,newzz)
{
    treelen <- length(tree)
    n0 <- nrow(newzz)
    asslist <- vector("list",treelen)
    asslist[[1]] <- 1:n0
    for(pos in 2:treelen)
    {
        if(is.null(tree[[pos]])) next
        ppos <- floor(pos/2)
        divar <- tree[[ppos]]$divar
        dlevels <- tree[[pos]]$dlevels
        idx <- which(newzz[,divar]%in%dlevels)
        asslist[[pos]] <- idx
    }
    return(asslist)
}

treePredict <- function(tree,model,newxx,newzz)
{
    asslist <- associate(tree,newzz)
    treelen <- length(tree)
    TT <- nrow(tree[[1]]$resp)
    n0 <- nrow(newxx)
    pmean <- matrix(0,nrow=TT,ncol=n0)
    for(pos in 1:treelen)
    {
        cnxx <- newxx[asslist[[pos]],]
        modelable <- tree[[pos]]$modelable
        cnzz <- NULL
        if(any(modelable))
        {
            cnzz <- newzz[asslist[[pos]],tree[[pos]]$modelable,drop=FALSE]
            cnzz <- apply(cnzz,2,adjlevel)
        }
        py <- predict(model[[pos]],cnxx,cnzz)
        pmean[,asslist[[pos]]] <- pmean[,asslist[[pos]]]+py$mean
    }
    return(pmean)
}
