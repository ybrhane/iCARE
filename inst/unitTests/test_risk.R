test_risk <- function() {

  risks <- c(0.0814961, 0.0771146, 0.1295784)

  data(bc_data, package="iCARE")
  set.seed(1)
  results <- computeAbsoluteRisk(model.formula     = caco ~ famhist + as.factor(parity), 
                                         model.cov.info    = bc_model_cov_info,
                                         model.snp.info    = bc_15_snps,
                                         model.log.RR      = bc_model_log_or,
                                         model.ref.dataset = ref_cov_dat,
                                         model.disease.incidence.rates   = bc_inc,
                                         model.competing.incidence.rates = mort_inc, 
                                         model.bin.fh.name = "famhist",
                                         apply.age.start    = 50, 
                                         apply.age.interval.length = 30,
                                         apply.cov.profile  = new_cov_prof,
                                         apply.snp.profile  = new_snp_prof, 
                                         return.refs.risk   = TRUE)
  ret <- results$risk
  
  
  checkEqualsNumeric(ret, risks, tolerance=1.0e-4)
  
}
