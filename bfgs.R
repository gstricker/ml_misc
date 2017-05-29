## Data
#######

library(fastGenoGAM)
## specify folder and experiment design path
folder <- "/s/project/coreProm/Michi/thornton/align_STAR"
config <- "/data/ouga/home/ag_gagneur/strickeg/workspace/analysis/diffBinding/config.txt"

## specify chunk and overhang size, bpk being the knot spacing
bpknots <- 20
chunkSize <- 20000
overhangSize <- 50*bpknots

## parallel backend
##BiocParallel::register(BiocParallel::SnowParam(workers = 4, progressbar = TRUE))

## build the GenoGAMDataSet
ggd <- GenoGAMDataSet(
  config, directory = folder,
  chunkSize = chunkSize, overhangSize = overhangSize,
  design = ~ s(x) + s(x, by = genotype)
)

## compute size factors
ggd <- computeSizeFactors(ggd)

## restrict to the actual small data region (usually not needed, just for vignette purposes)
subggd <- subset(ggd, seqnames == "chrXIV" & pos >= 300001 & pos <= 360000)

lambda <- 4601
theta <- 4.51
family <- "nb"
H <- 0
order <- 2
m <- 2

## compute model without parameter estimation to save time in vignette
coords <- fastGenoGAM:::.getCoordinates(subggd)
chunks <- fastGenoGAM:::.getChunkCoords(coords)
subggs <- fastGenoGAM:::setupGenoGAM(subggd, lambda = lambda, theta = theta, family = family, 
                    H = H, bpknots = bpknots, order = order,
                    penorder = m)

## Functions
############

.ll_nb <- function(beta,X,y,offset,theta,lambda,S){
  n <- dim(X)[1]
  eta <- offset + X%*%beta
  mu <- exp(eta)
  aux1 <- theta + y
  aux2 <- theta + mu
  ## pull the log inside to use gamma and factorial in log space due to 
  ## possibly very high numbers
  l <- sum(lgamma(aux1) - (lfactorial(y) + lgamma(theta))) + t(y) %*% eta + n*theta*log(theta) - t(aux1) %*% log(aux2)
  pen <- t(beta) %*% S %*% beta
  return(l[1]-lambda*pen[1,1])
}  

.grad_ll_nb <- function(beta,X,y,offset,theta,lambda,S){
  eta <- offset + X%*%beta
  mu <- exp(eta)
  z <- (y-mu)/(1+mu/theta)
  res <- t(X)%*%z
  pen <- S %*% beta
  return (res[,1]-2*lambda*pen[,1])
}

setup <- fastGenoGAM:::.initiate(subggd, subggs, coords, 1)

betas <- slot(setup, "beta")
X <- slot(setup, "designMatrix")
y <- slot(setup, "response")
offset <- slot(setup, "offset")
params <- slot(setup, "params")
S <- slot(setup, "penaltyMatrix")

grad <- .grad_ll_nb(betas, X, y, offset, params$theta, params$lambda, S)
ll <- .ll_nb(betas, X, y, offset, params$theta, params$lambda, S)

## BFGS implementation
######################

library(Matrix)

Hupdate <- function(H, s, y) {
    ro <- as.numeric(1/(t(y) %*% s))
    I <- diag(nrow(H))
    l <- (I - ro * s %*% t(y))
    r <- (I - ro * y %*% t(s))
    res <- l %*% H %*% r + ro * s %*% t(s)
    return(res)
}

Bupdate <- function(B, s, y) {
    m <- (B %*% s %*% t(s) %*% B)/as.vector((t(s) %*% B %*% s))
    r <- (y %*% t(y))/as.vector((t(y) %*% s))
    res <- B - m + r
    return(res)
}

linesearch <- function(x, alpha, p, ...) {

    fun <- function(alpha, pvec, beta, ...){
          n <- dim(X)[1]
          eta <- offset + X%*%(beta + alpha * pvec)
          mu <- exp(eta)
          aux1 <- theta + y
          aux2 <- theta + mu
          ## pull the log inside to use gamma and factorial in log space due to 
          ## possibly very high numbers
          l <- sum(lgamma(aux1) - (lfactorial(y) + lgamma(theta))) + t(y) %*% eta + n*theta*log(theta) - t(aux1) %*% log(aux2)
          pen <- t(beta + alpha * pvec) %*% S %*% (beta + alpha * pvec)
          return(l[1]-lambda*pen[1,1])
    }
    
    ## res <- optimize(f = fun, interval = c(1e-12, 2), pvec = p, beta = x, maximum = TRUE, tol = .Machine$double.eps, ...)
    ## return(res$maximum)
    res <- optim(alpha, f = fun, pvec = p, beta = x, method = "Nelder-Mead", control = list(fnscale = -1), ...)
    return(res$par)
}

bfgs <- function(x0, H0 = diag(length(x0)), epsilon, ffull, gradient, maxiter, alpha = 1, ...) {
    k <- 0
    x <- x0
    H <- H0
    lllast <- 0
    llk <- epsilon + 1
    traceDF <- data.frame()
    
    grad <- gradient(beta = x, ...)/(-1)
    lenGrad <- as.numeric(sqrt(t(grad) %*% grad))
    lldiff <- abs(llk - lllast)
    
    while(lldiff > epsilon) {
        cat("Iteration: ", k, "\n")
        cat("||fprime|| = ", lenGrad, "\n")
        grad <- matrix(gradient(beta = x, ...)/(-1), ncol = 1)
        p <- (-1) * H %*% grad
        alphak <- linesearch(x, alpha, p, ...)
        cat("Alpha = ", alphak, "\n")
        xnext <- x + alphak*p
        llk <- ffull(xnext, ...)
        lldiff <- abs(llk - lllast)
        cat("Likelihood = ", llk, "\n")
        
        tracek <- data.frame(iteration = k, fprime = lenGrad, alpha = alphak, ll = llk)
        traceDF <- rbind(traceDF, tracek)
        
        s = xnext - x
        gradNext <- matrix(gradient(beta = xnext, ...)/(-1), ncol = 1)
        yk <- gradNext - grad
        newH <- Hupdate(H, s, yk)
        k <- k + 1
        H <- newH
        x <- xnext
        alpha <- alphak
        lllast <- llk
        lenGrad <- as.numeric(sqrt(t(gradNext) %*% gradNext))
    }
    
    return(list(par = x, H = H, trace = traceDF))
}

initializeBeta <- function(id, data, init, coords) {
    tile <- coords[id]
    y <- assay(ggd)[IRanges::start(tile):IRanges::end(tile),]
    yhat <- endoapply(y, function(x) {
        runmean(x, k = 21, endrule = "constant")   
    })
    yhat <- unname(unlist(as.data.frame(yhat)))
    X <- slot(init, "designMatrix")
    XX <- t(X)%*%X
    b <- t(log(yhat + 1) %*% X)
    beta <- Matrix::solve(XX, b, sparse = TRUE)
    return(beta)
}

## apply bfgs

H <- t(X) %*% X
eigenDecompose <- eigen(H)
lipschitz <- 1/max(eigenDecompose$values)
Hinv <- Matrix::solve(H)

initBetas <- initializeBeta(1, subggd, subggs, coords)
slot(setup, "beta") <- as.matrix(initBetas, ncol = 1)
slot(setup, "fits") <- fastGenoGAM:::.getFits(setup)
H <- fastGenoGAM:::.compute_hessian_negbin(setup)
Hinv <- Matrix::solve(H)

res <- bfgs(x0 = initBetas, H0 = as.matrix(Hinv), epsilon = 1e-8, ffull = .ll_nb, gradient = .grad_ll_nb, maxiter = 10, alpha = lipschitz/2200,
            X = X, y = y, offset = offset, theta = params$theta, lambda = params$lambda, S = S)

## save(res_identity, file = "bfgs_Hidentity.rda")
## res_identity = initialize H with identity matrix
## res_inv = initialize H with H inverse matrix
## res_init = initialize betas with partial regression and H with the inverse of H
