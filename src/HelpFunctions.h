#ifndef HelpFunctions_H
#define HelpFunctions_H


#ifndef ARMA_USE_CXX11
#define ARMA_USE_CXX11
#endif

#ifndef ARMA_DONT_PRINT_ERRORS
#define ARMA_DONT_PRINT_ERRORS
#endif


#include "RcppArmadillo.h"

//////////////////////////
// Small help functions //
//////////////////////////


struct orderRank
{
  arma::uvec rankVector;
  arma::uvec orderVector;
};

orderRank rank(arma::vec& v);


orderRank rankwTiebreak(arma::vec& v, arma::vec& tieBreaker);

arma::umat getInterPerm(const arma::vec& X, const arma::vec& Y,
                        double theta_lo, double theta_hi, int open);

#endif
