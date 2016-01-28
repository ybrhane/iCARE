% --- Source file: bc_data.Rd ---
\name{bc_data}
\alias{bc_data}
\alias{bc_15_snps}
\alias{bc_inc}
\alias{bc_model_cov_info}
\alias{bc_model_log_or}
\alias{bc_model_log_or_post_50}
\alias{mort_inc}
\alias{new_cov_prof}
\alias{ref_cov_dat}
\alias{new_snp_prof}
\alias{ref_cov_dat_post_50}


%\docType{package}
\title{
  Data for examples
}
\description{
 Data for \code{\link{computeAbsoluteRisk}} and \code{\link{computeAbsoluteRiskSplitInterval}} 
}
\details{ 
"mort_inc" contains age-specific incidence rates of all-cause mortality from reference (1) below \cr
"bc_15_snps" contains published SNP information from reference (2)  \cr
 "bc_inc" contains age-specific incidence rates of breast cancer from reference (3) \cr 
"ref_cov_dat" contains a subsample of data imputed using reference (4) and (5) 
}
\examples{

 temp <- data(bc_data, package="iCARE")

 # Display the object names
 temp

}
\keyword{data}

\references{
(1)  Centers for Disease Control and Prevention (CDC), National Center for Health Statistics (NCHS). Underlying Cause of Death 1999-2011 on CDC WONDER Online Database, released 2014. Data are from the Multiple Cause of Death Files, 1999-2011, as compiled from data provided by the 57 vital statistics jurisdictions through the Vital Statistics Cooperative Program. \cr
 Accessed at http://wonder.cdc.gov/ucd-icd10.html on Aug 26, 2014.

(2)  Michailidou K, Beesley J, Lindstrom S, et al. Genome-wide association analysis of more than 120,000 individuals identifies 15 new susceptibility loci for breast cancer. Nature genetics 2015;47:373-80.

(3)  Surveillance, Epidemiology, and End Results (SEER) Program SEER*Stat Database: Incidence - SEER 18 Regs Research Data, Nov 2011 Sub, Vintage 2009 Pops (2000-2009) <Katrina/Rita Population Adjustment> - Linked To County Attributes - Total U.S., 1969-2010 Counties. In: National Cancer Institute D, Surveillance Research Program, Surveillance Systems Branch, ed. SEER18 ed.

(4)  2010 National Health Interview Survey (NHIS) Public Use Data Release, NHIS Survey Description. 2011. \cr
Accessed at ftp://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/NHIS/2010/srvydesc.pdf.)

(5)  Centers for Disease Control and Prevention (CDC). National Center for Health Statistics (NCHS). National Health and Nutrition Examination Survey Questionnaire. Hyattsville, MD: U.S. Department of Health and Human Services, Centers for Disease Control and Prevention; 2010.

}


