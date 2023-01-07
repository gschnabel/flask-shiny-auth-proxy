library(MCMCpack)
library(mvtnorm)

setClass("probdist",
         contains="VIRTUAL")

#' Draw a sample
#'
#' @param dist a distribution object 
#' @param n number of samples to draw from distribution
#'
#' @return a \code{p x n} matrix, each column contains an observation vector
#' @export
#' 
setGeneric("getSample", 
           function(dist,n)
             standardGeneric("getSample"))

#' Get density
#'
#' @param dist a distribution object
#' @param x contains \code{p x n} matrix where columns are observation vectors
#' @param log if TRUE return logarithm of density
#'
#' @return a vector with the (log) densities for the observation vectors
#' @export
#'
setGeneric("getDensity", 
           function(dist, x, log=TRUE, vars=NULL)
             standardGeneric("getDensity"))

#' Maximum Likelihood Estimate
#'
#' @param dist a distribution object
#' @param x contains \code{p x n} matrix where columns are observation vectors
#' @param weights importance of each observation
#' @param ... distribution specific tuning parameters
#'
#' @return a distribution object with optimized distribution parameters
#' @export
#'
setGeneric("getMLEstimate",
           function(dist, x, weights=rep(1,ncol(x)), ...)
             standardGeneric("getMLEstimate")
)


#' Get Posterior Sample
#'
#' @param dist a distribution object
#' @param x contains \code{p x n} matrix where columns are observation vectors
#' @param n number of samples
#' @param ... distribution specific tuning parameters
#'
#' @return a sample from the posterior distribution (via Gibbs sampling)
#' @export
#'
setGeneric("getPosteriorSample",
           function(dist, x, n, ...)
             standardGeneric("getPosteriorSample")
)


#' Evaluate Bayesian information criterion (BIC)
#'
#' @param dist a distribution object
#' @param x a \code{p x n} matrix where each column is an observation vector
#' @param weights importance of each observation
#'
#' @return the value of the Bayesian information criterion
#' @note If not all weights are equal one, the effective sample size (ESS)
#'       is used as the number of data points
#' @export
#'
setGeneric("getBIC",
           function(dist, x, weights=rep(1,ncol(x)))
             standardGeneric("getBIC")
)


#' Get membership probability
#'
#' @param dist a distribution object
#' @param x contains \code{p x n} matrix where columns are observation vectors
#' @param log if TRUE returns logarithm of membership probabilities
#' 
#' @return a \code{m x n} matrix where each column is associated with an observation vector
#'         and contains for each of the \code{m} components the probability of membership 
#'      
#' @note Only distributions which are mixtures of distributions implement this function.
#' @export
#'
setGeneric("getMembership",
           function(dist, x, log=TRUE) 
             standardGeneric("getMembership")
)

#' Create mixture distribution
#'
#' @param prop vector of probabilities for the components
#' @param compList list of distributions
#'
#' @return a distribution object containing a mixture distribution
#' @export
#'
createDist_Mix <- function(prop=numeric(0),compList=list(),
                           prior = list()) {
  
  if (length(compList)==0) {
    compDim <- 0L
    prop <- 0
  }
  else {
    compDim <- nrow(getSample(compList[[1]],1L))
    prop <- prop / sum(prop)
  }
  new("mixdist",prop=prop, comp=compList, dim=compDim,
      prior = prior)
}


setClass("mixdist",
         slots=list(
           prop="numeric",
           comp="list",
           dim="numeric",
           prior = 'list'
         ),
         prototype=list(
           prop=numeric(0),
           comp=list(),
           dim=0L,
           prior = list()
         ),
         validity=function(object) {
           isTRUE(length(object@prop)==length(object@comp) &&
                    all(diff(sapply(object@comp,function(x) nrow(getSample(x,1L))))==0))
         },
         contains = "probdist")


setMethod("getSample",
          signature=list(
            dist = "mixdist",
            n = "numeric"
          ),
          definition=function(dist, n) {
            
            numComp <- length(dist@prop)
            compIdx <- sample(numComp,n,prob=dist@prop,replace=TRUE)
            res <- matrix(0,nrow=dist@dim, ncol=n)
            for (curComp in seq(numComp)) {
              curIdx <- which(compIdx==curComp)
              if (length(curIdx)>0)
                res[,curIdx] <- getSample(dist@comp[[curComp]],length(curIdx))
            }
            res
          })

setMethod("getDensity",
          signature=list(
            dist = "mixdist",
            x = "matrix"
          ),
          definition=function(dist, x, log=TRUE, vars=NULL) {
            numComp <- length(dist@prop)
            logProp <- log(dist@prop)
            res <- matrix(0,nrow=numComp,ncol=ncol(x))
            for (curComp in seq(numComp)) {
              res[curComp,] <- logProp[curComp] + getDensity(dist@comp[[curComp]], x, log=TRUE, vars)
            }
            res <- apply(res,2,function(x) {
              xmax <- max(x)
              log(sum(exp(x - xmax))) + xmax
            })
            if (!isTRUE(log)) res <- exp(res)
            res
          })

setMethod("getMLEstimate",
          signature=list(
            dist = "mixdist",
            x = "matrix"
          ),
          definition=function(dist, x, weights, maxIter=50) {
            
            numComp <- length(dist@prop)
            for (curIter in seq(maxIter)) {
              # expectation
              memShip <- getMembership(dist, x, log=FALSE)
              memShip <- t(weights*t(memShip))
              dist@prop <- rowSums(memShip) / sum(weights)
              # maximization
              for (curComp in seq(numComp))
                dist@comp[[curComp]] <- getMLEstimate(dist@comp[[curComp]], x, memShip[curComp,])
            }
            dist
          })


setMethod("getMembership",
          signature=list(
            dist = "mixdist",
            x = "matrix"
          ),
          definition=function(dist, x, log=TRUE) {
            numComp <- length(dist@prop)
            logProp <- log(dist@prop)
            res <- matrix(0,nrow=numComp,ncol=ncol(x))
            for (curComp in seq(numComp))
              res[curComp,] <- logProp[curComp] + getDensity(dist@comp[[curComp]], x, log=TRUE)
            res <- apply(res,2,function(x) x-max(x))
            
            if (isTRUE(log)) {
              res <- apply(res,2,function(x) x - log(sum(exp(x))))
            }
            else
              res <- apply(exp(res),2,function(x) x/sum(x))
            res
          })


setMethod("getPosteriorSample",
          signature=list(
            dist = "mixdist",
            x = "matrix",
            n = "numeric"
          ),
          definition=function(dist, x, n) {
            
            alpha <- dist@prior$alpha
            numComp <- length(dist@prop)
            distList <- replicate(n, NULL, simplify = FALSE)
            for (i in seq(n)) {
              # sample membership
              memShip <- getMembership(dist, x, log=FALSE)
              numComp <- nrow(memShip)
              z <- rep(0, ncol(x))
              for (j in seq_along(z)) {
                z[j] <- sample.int(numComp, 1, prob = memShip[,j])
              }
              # sample mvn components
              for (curComp in seq(numComp)) {
                selectedPoints <- x[,z==curComp,drop=FALSE]
                dist@comp[[curComp]] <- getPosteriorSample(dist@comp[[curComp]], 
                                                           selectedPoints, 1)[[1]]
              }
              # sample from dirichlet posterior
              counts <- colSums(outer(z, 1:numComp, `==`))
              dist@prop <- as.vector(rdirichlet(1, alpha + counts))
              # save
              distList[[i]] <- dist
            }
            distList
          })


#' Create Multivariate Normal Distribution
#'
#' @param mean the center vector
#' @param sigma the covariance matrix
#'
#' @return a distribution object of a multivariate normal distribution
#' @import MCMCpack mvtnorm
#' @export
#'
createDist_MVN <- function(mean, sigma, prior = list()) {
  L <- chol(sigma)
  logDetSigma <- 2*sum(log(diag(L)))
  invL <- backsolve(chol(sigma),diag(ncol(sigma)))
  new("mvndist",mean=mean, sigma=sigma, L=L, invL=invL,
      logDetSigma=logDetSigma, prior=prior)
}


setClass("mvndist",
         slots=list(
           mean="numeric",
           sigma="matrix",
           L="matrix",
           invL="matrix",
           logDetSigma="numeric",
           prior = 'list'
         ),
         prototype=list(
           mean=0,
           sigma=matrix(1),
           L=matrix(1),
           invL=matrix(1),
           logDetSigma=1,
           prior = list()
         ),
         validity=function(object) {
           all(dim(object@sigma)==length(object@mean))
         },
         contains="probdist")



setMethod("getSample",
          signature=list(
            dist = "mvndist",
            n = "numeric"
          ),
          definition=function(dist, n) {
            numDim <- length(dist@mean)
            randVecs <- matrix(rnorm(numDim*n),ncol=n)
            t(dist@L) %*% randVecs + dist@mean
          })


setMethod("getDensity",
          signature=list(
            dist = "mvndist",
            x = "matrix"
          ),
          definition=function(dist, x, log=TRUE, vars=NULL) {
            if (is.null(vars)) {
              numDim <- length(dist@mean)
              stopifnot(nrow(x)==numDim)
              logNormConst <- (-numDim/2)*log(2*pi) - 0.5*dist@logDetSigma
              res <- logNormConst - 0.5*colSums((t(dist@invL) %*% (x-dist@mean))^2)
              if (!log) res <- exp(res)
              res
            }
            else
            {
              subX <- x
              subMean <- dist@mean[vars]
              subSigma <- dist@sigma[vars,vars,drop=FALSE]
              marDist <- createDist_MVN(subMean,subSigma)
              getDensity(marDist, subX, log, vars=NULL)
            }
          })



setMethod("getMLEstimate",
          signature=list(
            dist = "mvndist",
            x = "matrix"
          ),
          definition=function(dist, x, weights) {
            pars <- cov.wt(t(x), wt=weights)
            createDist_MVN(pars$center, pars$cov)
          })

setMethod("getPosteriorSample",
          signature=list(
            dist = "mvndist",
            x = "matrix",
            n = "numeric"
          ),
          definition=function(dist, x, n) {
            
            # prior
            prior <- dist@prior
            mu0 <- prior$mu0
            kappa0 <- prior$kappa0
            nu0 <- prior$nu0
            phi <- prior$phi
            # data
            numObs <- ncol(x)
            D <- cov.wt(t(x), method = "ML")
            # posterior
            mu1 <- (kappa0*mu0 + numObs*D$center) / (kappa0 + numObs)
            kappa1 <- kappa0 + numObs
            nu1 <- nu0 + numObs
            phi1 <- phi + numObs*D$cov + (kappa0*numObs) / (kappa0+numObs) * 
              tcrossprod(D$center - mu0)
            # create samples
            distList <- replicate(n, NULL, simplify = FALSE)
            for (i in seq(n)) {
              sampleCov <- riwish(nu1, phi1)
              sampleMean <- as.vector(rmvnorm(1, mu1, 1/kappa1 * sampleCov))
              distList[[i]] <- createDist_MVN(sampleMean, sampleCov, prior = dist@prior)
            }
            distList
          })


#' Create Constant Distribution
#'
#' @return a distribution object of a constant distribution
#' @export
#'
createDist_Const <- function(box) {
  new("constdist",box=box)
}


setClass("constdist",
         slots=list(
           box = 'numeric'
         ),
         prototype=list(
           box = numeric(0)
         ),
         validity = function(object) {
           length(object@box) %% 2 == 0
         },
         contains="probdist")



setMethod("getSample",
          signature=list(
            dist = "constdist",
            n = "numeric"
          ),
          definition=function(dist, n) {
            numDim <- length(dist@box) / 2
            matrix(runif(n*numDim,
                      dist@box[1:numDim], 
                      dist@box[(numDim+1):(2*numDim)]),
                   nrow = numDim)
          })


setMethod("getDensity",
          signature=list(
            dist = "constdist",
            x = "matrix"
          ),
          definition=function(dist, x, log=TRUE, vars=NULL) {
            numDim <- length(dist@box) / 2
            mins <- dist@box[1:numDim]
            maxs <- dist@box[(numDim+1):(2*numDim)]
            lens <- abs(maxs - mins)
            if (is.null(vars)) vars <- seq_len(numDim)
            logfun <- get("log", pos = "package:base")
            if (isTRUE(log))
              -sum(logfun(lens[vars]))
            else
              1/prod(lens[vars])
          })



setMethod("getMLEstimate",
          signature=list(
            dist = "constdist",
            x = "matrix"
          ),
          definition=function(dist, x, weights) {
            createDist_Const(dist@box)
          })


setMethod("getPosteriorSample",
          signature=list(
            dist = "constdist",
            x = "matrix",
            n = "numeric"
          ),
          definition=function(dist, x, n) {
            
            # create samples
            distList <- replicate(n, NULL, simplify = FALSE)
            for (i in seq(n)) {
              distList[[i]] <- dist
            }
            distList
          })
