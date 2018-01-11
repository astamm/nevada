// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// dist_hamming_impl
double dist_hamming_impl(const arma::mat& x, const arma::mat& y);
RcppExport SEXP _nevada_dist_hamming_impl(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(dist_hamming_impl(x, y));
    return rcpp_result_gen;
END_RCPP
}
// dist_frobenius_impl
double dist_frobenius_impl(const arma::mat& x, const arma::mat& y);
RcppExport SEXP _nevada_dist_frobenius_impl(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(dist_frobenius_impl(x, y));
    return rcpp_result_gen;
END_RCPP
}
// dist_spectral_impl
double dist_spectral_impl(const arma::mat& x, const arma::mat& y);
RcppExport SEXP _nevada_dist_spectral_impl(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(dist_spectral_impl(x, y));
    return rcpp_result_gen;
END_RCPP
}
// dist_root_euclidean_impl
double dist_root_euclidean_impl(const arma::mat& x, const arma::mat& y);
RcppExport SEXP _nevada_dist_root_euclidean_impl(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(dist_root_euclidean_impl(x, y));
    return rcpp_result_gen;
END_RCPP
}
// dist_nvd_impl
arma::mat dist_nvd_impl(const Rcpp::List& z, const std::string distance);
RcppExport SEXP _nevada_dist_nvd_impl(SEXP zSEXP, SEXP distanceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const std::string >::type distance(distanceSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_nvd_impl(z, distance));
    return rcpp_result_gen;
END_RCPP
}
// repr_adjacency_impl
arma::mat repr_adjacency_impl(const unsigned int numberOfVertices, const arma::mat& edgeList, const arma::vec& weights);
RcppExport SEXP _nevada_repr_adjacency_impl(SEXP numberOfVerticesSEXP, SEXP edgeListSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type numberOfVertices(numberOfVerticesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type edgeList(edgeListSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(repr_adjacency_impl(numberOfVertices, edgeList, weights));
    return rcpp_result_gen;
END_RCPP
}
// stat_lot_impl
double stat_lot_impl(const arma::mat& distanceMatrix, const arma::vec& firstGroupIndices, const arma::vec& secondGroupIndices);
RcppExport SEXP _nevada_stat_lot_impl(SEXP distanceMatrixSEXP, SEXP firstGroupIndicesSEXP, SEXP secondGroupIndicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type distanceMatrix(distanceMatrixSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type firstGroupIndices(firstGroupIndicesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type secondGroupIndices(secondGroupIndicesSEXP);
    rcpp_result_gen = Rcpp::wrap(stat_lot_impl(distanceMatrix, firstGroupIndices, secondGroupIndices));
    return rcpp_result_gen;
END_RCPP
}
// stat_sot_impl
double stat_sot_impl(const arma::mat& distanceMatrix, const arma::vec& firstGroupIndices, const arma::vec& secondGroupIndices);
RcppExport SEXP _nevada_stat_sot_impl(SEXP distanceMatrixSEXP, SEXP firstGroupIndicesSEXP, SEXP secondGroupIndicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type distanceMatrix(distanceMatrixSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type firstGroupIndices(firstGroupIndicesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type secondGroupIndices(secondGroupIndicesSEXP);
    rcpp_result_gen = Rcpp::wrap(stat_sot_impl(distanceMatrix, firstGroupIndices, secondGroupIndices));
    return rcpp_result_gen;
END_RCPP
}
// stat_biswas_impl
double stat_biswas_impl(const arma::mat& distanceMatrix, const arma::vec& firstGroupIndices, const arma::vec& secondGroupIndices);
RcppExport SEXP _nevada_stat_biswas_impl(SEXP distanceMatrixSEXP, SEXP firstGroupIndicesSEXP, SEXP secondGroupIndicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type distanceMatrix(distanceMatrixSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type firstGroupIndices(firstGroupIndicesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type secondGroupIndices(secondGroupIndicesSEXP);
    rcpp_result_gen = Rcpp::wrap(stat_biswas_impl(distanceMatrix, firstGroupIndices, secondGroupIndices));
    return rcpp_result_gen;
END_RCPP
}
// stat_energy_impl
double stat_energy_impl(const arma::mat& distanceMatrix, const arma::vec& firstGroupIndices, const arma::vec& secondGroupIndices, const unsigned int alpha);
RcppExport SEXP _nevada_stat_energy_impl(SEXP distanceMatrixSEXP, SEXP firstGroupIndicesSEXP, SEXP secondGroupIndicesSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type distanceMatrix(distanceMatrixSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type firstGroupIndices(firstGroupIndicesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type secondGroupIndices(secondGroupIndicesSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(stat_energy_impl(distanceMatrix, firstGroupIndices, secondGroupIndices, alpha));
    return rcpp_result_gen;
END_RCPP
}
// stat_t_euclidean_impl
double stat_t_euclidean_impl(const Rcpp::List& x, const Rcpp::List& y, const bool pooled);
RcppExport SEXP _nevada_stat_t_euclidean_impl(SEXP xSEXP, SEXP ySEXP, SEXP pooledSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const bool >::type pooled(pooledSEXP);
    rcpp_result_gen = Rcpp::wrap(stat_t_euclidean_impl(x, y, pooled));
    return rcpp_result_gen;
END_RCPP
}
// stat_edge_count_impl
arma::vec stat_edge_count_impl(const arma::mat& E, const arma::vec& indices);
RcppExport SEXP _nevada_stat_edge_count_impl(SEXP ESEXP, SEXP indicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type E(ESEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type indices(indicesSEXP);
    rcpp_result_gen = Rcpp::wrap(stat_edge_count_impl(E, indices));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_nevada_dist_hamming_impl", (DL_FUNC) &_nevada_dist_hamming_impl, 2},
    {"_nevada_dist_frobenius_impl", (DL_FUNC) &_nevada_dist_frobenius_impl, 2},
    {"_nevada_dist_spectral_impl", (DL_FUNC) &_nevada_dist_spectral_impl, 2},
    {"_nevada_dist_root_euclidean_impl", (DL_FUNC) &_nevada_dist_root_euclidean_impl, 2},
    {"_nevada_dist_nvd_impl", (DL_FUNC) &_nevada_dist_nvd_impl, 2},
    {"_nevada_repr_adjacency_impl", (DL_FUNC) &_nevada_repr_adjacency_impl, 3},
    {"_nevada_stat_lot_impl", (DL_FUNC) &_nevada_stat_lot_impl, 3},
    {"_nevada_stat_sot_impl", (DL_FUNC) &_nevada_stat_sot_impl, 3},
    {"_nevada_stat_biswas_impl", (DL_FUNC) &_nevada_stat_biswas_impl, 3},
    {"_nevada_stat_energy_impl", (DL_FUNC) &_nevada_stat_energy_impl, 4},
    {"_nevada_stat_t_euclidean_impl", (DL_FUNC) &_nevada_stat_t_euclidean_impl, 3},
    {"_nevada_stat_edge_count_impl", (DL_FUNC) &_nevada_stat_edge_count_impl, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_nevada(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
