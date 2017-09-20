##RcppExports
##https://github.com/NMML/kfdnm/blob/8d1be2235a8d40a41d09abc9c423f77c23369a6f/R/RcppExports.R

N_trans <- function(fromN, omega, gamma, R, N_max) {
  .Call('kfdnm_N_trans', PACKAGE = 'kfdnm', fromN, omega, gamma, R, N_max)
}

N_trans_mat <- function(omega, gamma, R, N_max) {
  .Call('kfdnm_N_trans_mat', PACKAGE = 'kfdnm', omega, gamma, R, N_max)
}

dnm_hmm <- function(n, R, id, omega_dnm, gamma, p, N_max, back_sample) {
  .Call('kfdnm_dnm_hmm', PACKAGE = 'kfdnm', n, R, id, omega_dnm, gamma, p, N_max, back_sample)
}

kf_hmm <- function(Y, omega, id) {
  .Call('kfdnm_kf_hmm', PACKAGE = 'kfdnm', Y, omega, id)
}