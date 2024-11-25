library(simulateGP)
library(susieR)
library(here)
library(dplyr)

#' Convert log Bayes Factor to summary stats
#'
#' @param lbf p-vector of log Bayes Factors for each SNP
#' @param n Overall sample size
#' @param af p-vector of allele frequencies for each SNP
#' @param prior_v Variance of prior distribution. SuSiE uses 50
#'
#' @return tibble with lbf, af, beta, se, z
lbf_to_z_cont <- function(lbf, n, af, prior_v=50)
{
  se = sqrt(1 / (2 * n * af * (1-af)))
  r = prior_v / (prior_v + se^2)
  z = sqrt((2 * lbf - log(sqrt(1-r)))/r)
  beta <- z * se
  return(tibble(lbf, af, z, beta, se))
}

# Read in example LD matrix from simulateGP repository

map <- readRDS(url("https://github.com/explodecomputer/simulateGP/raw/refs/heads/master/inst/extdata/ldobj_5_141345062_141478055.rds", "rb"))
glimpse(map)

# Generate summary statistics for a single causal variant 

set.seed(1234)
ss <- map$map %>%
  generate_gwas_params(h2=0.003, Pi=1/nrow(.)) %>%
  generate_gwas_ss(50000, ld=map$ld)
table(ss$beta == 0)

# Run SuSiE

sout <- susie_rss(ss$bhat / ss$se, R = map$ld, n = 50000, bhat = ss$bhat, var_y=1)

summary(sout)

# Get z scores from lbf

a <- lbf_to_z_cont(sout$lbf_variable[1,], 50000, ss$af, prior_v = 50)
a

plot(z ~ lbf, a)

# new z vs original z

plot(a$z^2 ~ ss$fval)

# check concordance

summary(lm(a$z^2 ~ ss$fval))

# compare betas

plot(a$beta ~ ss$bhat, xlab="Original betas", ylab="Inferred betas from log Bayes factors")


