#' Risk Stratification from Polygenic Risk Score (PRS) based on Variance Explained
#'
#' This function is used to build absolute risk models and apply them to estimate absolute risks.
#'
#' @param apply.age.start single integer or vector of integer ages for the start of the interval over which to compute absolute risk.
#' @param apply.age.interval.length single integer or vector of integer years over which absolute risk should be computed.
#' @param model.disease.incidence.rates two column matrix [ integer ages, incidence rates] or three column matrix [start age, end age, rate] with incidence rate of disease. Must fully cover age interval for estimation.
#' @param use.c.code binary indicator of whether to run the c program for fast computation.
#' @param model.competing.incidence.rates two column matrix [ integer ages, incidence rates] or three column matrix [start age, end age, rate] with incidence rate of competing events. Must fully cover age interval for estimation.
#' @param variance.explained proportion of variance explained by PRS
#' @param N number of simulated PRS values
#' @param predict.stand.vals standardized values for which to make a risk prediction
#' @return This function returns risks based on simulated PRS based on SNPs with given \code{variance.explained}.
#' @export
#' @examples
#' print("Hello world")
PRS_research <- function(apply.age.start, apply.age.interval.length, model.disease.incidence.rates, use.c.code = 1,
                         model.competing.incidence.rates = NULL, variance.explained, N = 10000,
                         predict.stand.vals = NULL){

  if(length(apply.age.start)!=1 || length(apply.age.interval.length)!=1){
    return_and_print("ERROR:  You may only choose one value for apply.age.start and one value for apply.age.interval.length.")
    stop()
  }
  apply.age.start      <- rep(apply.age.start, N)
  apply.age.interval.length <- rep(apply.age.interval.length, N)
  lambda       <- model.disease.incidence.rates
  if(is.null(model.competing.incidence.rates)){  model.competing.incidence.rates= cbind(lambda[,1], rep(0, length(lambda[,1])))   }

  if(ncol(lambda)!=2){ return_and_print("ERROR: model.disease.incidence.rates should have 2 columns: [Ages,Rates]"); stop() }

  if(ncol(model.competing.incidence.rates)!=2){ return_and_print("ERROR: model.competing.incidence.rates should have 2 columns: [Ages,Rates]"); stop() }

  if(sum(lambda[,1]%%1)!=0){ return_and_print("ERROR: The first column of model.disease.incidence.rates should be integer ages."); stop()  }

  if(sum(model.competing.incidence.rates[,1]%%1)!=0){ return_and_print("ERROR: The first column of model.competing.incidence.rates should be integer ages."); stop()  }
  if( prod(is.element(seq(range(apply.age.start)[1], range(apply.age.start+apply.age.interval.length)[2]), lambda[,1])) == 0){
    return_and_print("ERROR: The 'model.disease.incidence.rates' input must have age-specific rates for each integer age covered by the prediction intervals defined by 'apply.age.start' and 'apply.age.interval.length.'  You must make these inputs consistent with one another to proceed.")
    stop()
  }
  if( prod(is.element(seq(range(apply.age.start)[1], range(apply.age.start+apply.age.interval.length)[2]), model.competing.incidence.rates[,1])) == 0){
    return_and_print("ERROR: The 'model.competing.incidence.rates' input must have age-specific rates for each integer age covered by the prediction intervals defined by 'apply.age.start' and 'apply.age.interval.length.'  You must make these inputs consistent with one another to proceed.")
    stop()
  }

  PRS.         <- rnorm(N)
  log.OR       <- variance.explained^0.5
  pop.dist.mat <- cbind(PRS.)
  beta_est     <- as.matrix(nrow = length(c( log.OR)), ncol = 1, c( log.OR) )
  Z_new        <- t(pop.dist.mat)
  pop.weights  <- rep(1/nrow(pop.dist.mat), length(pop.dist.mat))

  approx_expectation_rr = weighted.mean(  exp( pop.dist.mat%*% beta_est), w = pop.weights ) ## all equal weights
  calc = precise_lambda0(lambda, approx_expectation_rr, beta_est, pop.dist.mat, pop.weights)
  lambda_0 = calc[[1]]
  precise_expectation_rr = calc[[2]]
  pop.dist.mat = t(pop.dist.mat )

  ###### Compute A_j  ##
  final_risks <- comp_Aj(Z_new, apply.age.start, apply.age.interval.length, lambda_0, beta_est, model.competing.incidence.rates)

  if(is.null(predict.stand.vals)==0){

    if(is.vector(predict.stand.vals)==0){
      return_and_print("ERROR: predict.stand.vals must be a vector")
      stop()
    }else{
      temp = t(cbind(predict.stand.vals))
      specific_risks <- comp_Aj(temp, apply.age.start[1:length(predict.stand.vals)], 
          apply.age.interval.length[1:length (predict.stand.vals)], lambda_0, beta_est, model.competing.incidence.rates)
    }
  }
  results      <- list()
  results$risk <- final_risks
  if(is.null(predict.stand.vals)==0){
    results$predictions.for.stand.vals <- specific_risks
  }
  results
}

package_results <- function(final_risks, Z_new, covs_in_model, handle.snps, apply.age.start, apply.age.interval.length, apply.cov.profile ,
                            model.log.RR, beta_est, apply.snp.profile, snps.names, return.lp, lps ){
  result <-list()
  result$risk           <- cbind(as.vector(final_risks))
  colnames(result$risk) <- c("Risk_Estimate")

  if(covs_in_model){
    if(handle.snps==0){
      info = cbind(as.vector(apply.age.start), as.vector(apply.age.start + apply.age.interval.length), as.vector(final_risks), apply.cov.profile )
      colnames(info) = c("Int_Start", "Int_End", "Risk_Estimate", colnames(apply.cov.profile ))
      beta.names = rownames(model.log.RR)
    }else{
      info = cbind(as.vector(apply.age.start), as.vector(apply.age.start + apply.age.interval.length), as.vector(final_risks), cbind( apply.snp.profile, apply.cov.profile ))
      colnames(info) = c("Int_Start", "Int_End", "Risk_Estimate", snps.names, colnames(apply.cov.profile ))
      beta.names = c( snps.names, rownames(model.log.RR) )
    }
  }else{
    info = cbind(as.vector(apply.age.start), as.vector(apply.age.start + apply.age.interval.length), as.vector(final_risks), apply.snp.profile)
    colnames(info) = c("Int_Start", "Int_End", "Risk_Estimate", snps.names)
    beta.names =  snps.names
  }
  result$details = info
  beta.used = beta_est
  rownames(beta.used) <- beta.names
  colnames(beta.used) <- "log OR used"
  result$beta.used    <- beta.used

  if(return.lp==TRUE){
    result$lps <- lps
  }
  result
}

  
