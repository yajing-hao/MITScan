#' Perform score-based genome-wide association analysis using zero-inflated negative binomial models
#'
#' MITScan_ZINB() performs score-based genome-wide association analysis
#' across all OTU abundances and host genes using zero-inflated negative binomial models,
#' and returns a matrix where the (i,j) element of the matrix represents
#' the score test statistic for the i-th host gene and the j-th OTU abundance.
#'
#' @param Ymat A matrix of OTU abundances from 16S rRNA sequencing.
#' Rows and columns represent OTU names and samples, respectively.
#' @param X A matrix of covariates related to the count component.
#' @param Z A matrix of covariates related to the zero component.
#' @param A A matrix of host gene expression from RNA-Seq.
#' Rows and columns represent host gene names and samples, respectively.
#' @param perm An integer indicating the number of permutations for genomic control.
#' @param no_cores An integer indicating the number of cores for parallel calculations.
#' @param seed An integer specifying the seed for random number generation.
#'
#' @return A matrix of score test statistics across all OTU abundances and host genes.
#' Rows and columns represent host gene names and OTU names, respectively.
#' The last row contains the median of permuted score test statistics for each OTU.
#'
#' @examples
#' OTU = IBD_OTU
#' hostgene = IBD_hostgene
#' cvrt = IBD_cvrt
#' zinb1 = MITScan_ZINB(Ymat=OTU, X=cvrt, Z=NULL, A=hostgene, perm=2, no_cores=8, seed=1)
#'
#' @import pscl
#' @import parallel
#'
#' @export
MITScan_ZINB = function(Ymat, X=NULL, Z=NULL, A, perm, no_cores=64, seed=1) {
  library(parallel)
  Ymat = as.matrix(round(Ymat))
  if (is.null(X)) {
    X = matrix(1, nrow=ncol(Ymat), ncol=1)
  } else {
    X = model.matrix(~., data.frame(X))
  }
  if (is.null(Z)) {
    Z = matrix(1, nrow=ncol(Ymat), ncol=1)
  } else {
    Z = model.matrix(~., data.frame(Z))
  }
  A = as.matrix(A)
  A = scale(t(A), center=TRUE, scale=FALSE)
  set.seed(seed)
  permA = function(A) {A[sample(1:nrow(A), replace=FALSE),]}
  Ap = do.call(cbind, replicate(perm, permA(A), simplify=FALSE))

  clust = makeCluster(no_cores)
  clusterEvalQ(clust, library(pscl))
  res = parApply(clust, Ymat, 1, scorep.zinb, X, Z, A, Ap)
  stopCluster(clust)
  rownames(res) = c(unlist(colnames(A)), "GC")
  return(res)
}



#' Compute score test statistics for a given OTU across host genes using zero-inflated negative binomial models
#'
#' scorep.zinb() computes score test statistics for a given OTU
#' across all host genes using zero-inflated negative binomial models,
#' and return a vector where the i-th element of the vector represents
#' the score test statistic for the i-th host gene and the given OTU.
#'
#' @param Y A vector of abundances for a given OTU.
#' @param X A matrix of covariates related to the count component.
#' @param Z A matrix of covariates related to the zero component.
#' @param A A matrix of host gene expression from RNA-Seq.
#' Rows and columns represent host gene names and samples, respectively.
#' @param Ap A matrix of permuted host gene expression.
#'
#' @return A vector of score test statistics for a given OTU across all host genes
#' and the median of permuted score test statistics.
#'
#' @import pscl
#'
#' @export
scorep.zinb = function(Y, X, Z, A, Ap) {
  error0 =  tryCatch({
    zinb1 = zeroinfl(Y ~ 0 + X | 0 + Z, dist="negbin", link="logit")
  }, error = function(e) {e})

  if (inherits(error0, "error")) {
    return(rep(NA, ncol(A)+1))
  } else {
    xp = ncol(X)
    zp = ncol(Z)
    mui = as.numeric(exp(X %*% coef(zinb1)[1:xp]))
    lambdai = as.numeric(exp(Z %*% coef(zinb1)[(xp+1):(xp+zp)]))
    theta1 = zinb1[["theta"]]
    t1 = theta1 / (theta1 + mui)
    t2 = t1^theta1

    Sbeta = (Y==0) * (-mui*t1*t2/(lambdai+t2)) +  (Y>0) * (Y-mui)*t1
    S = colSums(Sbeta*A)

    Dbeta = ifelse(Y==0, -((mui-1)*lambdai*t1^2*t2-t1^2*t2^2)*mui/(lambdai+t2)^2,
                   theta1*(theta1+Y)*mui/(theta1+mui)^2)
    Dbg = ifelse(Y==0, -t1*t2*mui*lambdai/(lambdai+t2)^2, 0)
    Dbt = ifelse(Y==0, -(lambdai*t1*t2*log(1/t1) - mui/theta1*t1*t2*(lambdai*(theta1+1)/(theta1+mui) + t1*t2/theta1))*mui/(lambdai+t2)^2,
                 -(Y-mui)*mui/(theta1+mui)^2)
    Dbeta = ifelse(is.na(Dbeta), 0, Dbeta)
    Dbg = ifelse(is.na(Dbg), 0, Dbg)
    Dbt = ifelse(is.na(Dbt), 0, Dbt)

    i1 = -solve(zinb1[["optim"]][["hessian"]])
    i2 = colSums(Dbeta*A^2)
    i3 = t(X) %*% diag(Dbeta) %*% A
    i4 = t(Z) %*% diag(Dbg) %*% A
    i5 = colSums(theta1*Dbt*A)
    i6 = rbind(i3,i4,i5)
    i7 = apply(i6, 2, function(x) {t(x) %*% i1 %*% x})
    I = i2 - i7
    I = ifelse(I>0, I, NA)

    Sp = colSums(Sbeta*Ap)
    i2p = colSums(Dbeta*Ap^2)
    i3p = t(X) %*% diag(Dbeta) %*% Ap
    i4p = t(Z) %*% diag(Dbg) %*% Ap
    i5p = colSums(theta1*Dbt*Ap)
    i6p = rbind(i3p,i4p,i5p)
    i7p = apply(i6p, 2, function(x) {t(x) %*% i1 %*% x})
    Ip = i2p - i7p
    Ip = ifelse(Ip>0, Ip, NA)

    m1 = S^2/I
    m2 = Sp^2/Ip
    m2 = m2[!is.na(m2)]
    return(c(round(m1, 2), median(m2)))
  }
}
