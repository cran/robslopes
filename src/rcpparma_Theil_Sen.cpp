// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
// #include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// include functions to construct permutations
#include "HelpFunctions.h"

//////////////////////
// Count inversions //
//////////////////////

void merge_TS(arma::uvec& y, arma::uword left, arma::uword middle, arma::uword right, arma::uword& invCount)
{
  // c is number of inversions
  arma::uword i, j, k;
  const arma::uword n1 = middle - left + 1;
  const arma::uword n2 = right - middle;
  
  arma::uvec Left = y.subvec(left, left + n1 - 1);
  arma::uvec Right = y.subvec(middle + 1, middle + n2);
  
  
  i = 0, j = 0, k = left;
  while (i < n1 && j < n2)
  {
    if (Left(i) <= Right(j))
    {
      y(k++) = Left(i++);
    }
    else
    {//copy from the right an increase inversion count
      invCount+= n1 - i;
      y(k++) = Right(j++);
    }
  }
  
  while (i < n1)
  {
    y(k++) = Left(i++);
  }
  
  while (j < n2)
  {
    y(k++) = Right(j++);
  }
}


void mergeSort_TS(arma::uvec& y, arma::uword left, arma::uword right, arma::uword& invCount)
{
  // assumes y is an vector of distinct integers 0, 1, ... n-1 without gaps
  // 
  if (left < right)
  {
    arma::uword middle = left + (right - left) / 2;
    mergeSort_TS(y, left, middle, invCount);
    mergeSort_TS(y, middle + 1, right, invCount);
    merge_TS(y, left, middle, right, invCount);
  }
}


arma::uword countInversions(arma::uvec y)
{
  // calculate number of inversions in vector y
  // y is assumed to contain distinct integers from 1 to n without gaps
  // algorithm uses a modification of merge sort
  
  int n = y.size();
  arma::uword result = 0;
  mergeSort_TS(y, 0, n - 1, result);
  return(result);
}




//////////////////////////////////
// Sample intersection abscissa //
//////////////////////////////////

void merge2_TS(arma::uvec& y,arma::uword left, arma::uword middle, arma::uword right,
               arma::vec& sampleInds,
               const arma::vec& X, const arma::vec& Y, arma::uword& invCount,
               arma::uword& sampleInds_idx,
               const arma::uvec& valToline,
               double theta_lo, double theta_hi)
{
  // c is number of inversions
  arma::uword i, j, k;
  const arma::uword n1 = middle - left + 1;
  const arma::uword n2 = right - middle;
  
  arma::uvec Left = y.subvec(left, left + n1 - 1);
  arma::uvec Right = y.subvec(middle + 1, middle + n2);
  arma::uword newinvCount = 0;
  
  i = 0, j = 0, k = left;
  while (i < n1 && j < n2)
  {
    if (Left(i) <= Right(j))
    {
      y(k++) = Left(i++);
    }
    else
    {//copy from the right and increase inversion count
      newinvCount = invCount + n1 - i;
      
      while(sampleInds(sampleInds_idx) < newinvCount) {
        
        int itemp = i + sampleInds(sampleInds_idx) - invCount;
        
        double intersection;
        if (arma::is_finite(theta_lo) && arma::is_finite(theta_hi))
        {
          if (X(valToline(Right(j))) == X(valToline(Left(itemp)))) 
          {
            intersection = arma::datum::inf;
          }
          else
          {
            // use more robust formula to compute the sampled intersection ordinate:
            double delta_lo = (theta_lo * X(valToline(Right(j))) -
                               Y(valToline(Right(j)))) -
                               (theta_lo * X(valToline(Left(itemp))) -
                               Y(valToline(Left(itemp))));
            double delta_hi = (theta_hi * X(valToline(Right(j))) -
                               Y(valToline(Right(j)))) -
                               (theta_hi * X(valToline(Left(itemp))) -
                               Y(valToline(Left(itemp))));
            intersection = theta_lo + delta_lo / (delta_lo - delta_hi) * (theta_hi - theta_lo);
          }
        }
        else
        { // with infinite bounds, use non-robust formula
          if (X(valToline(Right(j))) == X(valToline(Left(itemp)))) 
          {
            intersection = arma::datum::inf;
          }
          else
          {
            intersection =(Y(valToline(Right(j))) - Y(valToline(Left(itemp)))) /
              (X(valToline(Right(j))) - X(valToline(Left(itemp))));
          } 
        }
        
        sampleInds(sampleInds_idx) = intersection;
        sampleInds_idx++;
      }
      
      invCount = newinvCount;
      y(k++) = Right(j++);
    }
  }
  
  while (i < n1)
  {
    y(k++) = Left(i++);
  }
  
  while (j < n2)
  {
    y(k++) = Right(j++);
  }
}



void mergeSort2_TS(arma::uvec& y, arma::uword left, arma::uword right, arma::vec& sampleInds,
                   const arma::vec& X, const arma::vec& Y, arma::uword& invCount,
                   arma::uword& sampleInds_idx, const arma::uvec& valToline,
                   double theta_lo, double theta_hi)
{
  // assumes y is an vector of distinct integers 0, 1, ... n-1 without gaps
  // 
  if (left < right)
  {
    arma::uword middle = left + (right - left) / 2;
    mergeSort2_TS(y, left, middle, sampleInds, X, Y, invCount,
                  sampleInds_idx, valToline, theta_lo, theta_hi);
    mergeSort2_TS(y, middle + 1, right, sampleInds, X, Y, invCount,
                  sampleInds_idx, valToline, theta_lo, theta_hi);
    merge2_TS(y, left, middle, right, sampleInds, X, Y, invCount,
              sampleInds_idx, valToline, theta_lo, theta_hi);
  }
}


void sampleIA(arma::uvec y, arma::vec& sampleInds, 
              const arma::vec& X, const arma::vec& Y,
              const arma::uvec valToline,
              double theta_lo, double theta_hi)
{
  // sampleInds contains the slopes after 
  int n = y.size();
  arma::uword invCount = 0;
  arma::uword sampleInds_idx = 0;
  mergeSort2_TS(y, 0, n - 1, sampleInds, X, Y, invCount, 
                sampleInds_idx, valToline, theta_lo, theta_hi);
}


///////////////
// Theil Sen //
///////////////

// [[Rcpp::export]]
arma::vec rcpp_TheilSen(const arma::vec X, const arma::vec Y, const bool verbose,
                        const arma::uword medind, const arma::uword medind2) {
  // Computes repeated median slope and intercept
  
  // global constants
  const arma::uword n = X.size();
  const arma::uword c = 20; // stop contraction iterations when C<= c * n
  const arma::uword r = n; // number of random ordinates to sample each iteration
  const arma::uword t = 3; // tuning parameter used in k_lo and k_hi
  
  // Initialization
  int termination = 0;
  double theta_lo = - arma::datum::inf;
  double theta_hi =  arma::datum::inf;
  double thetaP_lo, thetaP_hi = 0;
  arma::uword C = (arma::uword)(n * (n - 1) / 2) ; // number of ordinates in (theta_lo,theta_hi]
  
  arma::uvec IOcounter(2); // used to check if C goes still goes down
  IOcounter(1) = C;
  IOcounter(0) = 2 * IOcounter(1) ;
  arma::uword L = 0; // elements to the left of our interval
  arma::vec inds;
  arma::umat permutation(n, 3);
  
  int iterationCounter = 0;
  
  if (verbose)
  {
    Rcpp::Rcout << "Initialization finished, starting interval contraction." << std::endl;
  }
  
  while (C > (arma::uword)(c * n)) //switch to brute force otherwise
  {
    // sample ordinates
    
    // generate a sorted list of n random numbers, with replacement
    // we avoid arma::randi, since it has maximum_integer as upper bound
    inds = arma::sort(arma::randu<arma::vec>(r + 1)) * double(C);
    inds = arma::floor(inds);
    inds(r) = C + 1; //sentinel value
    permutation = getInterPerm(X, Y, theta_lo, theta_hi, 0);
    sampleIA(permutation.col(0), inds, X, Y, permutation.col(1), theta_lo, theta_hi);
    inds = arma::sort(inds);
    arma::uword k = medind - L;
    arma::uword k_lo = std::max(1.0, std::floor((double)k * (double)n / (double)C - t*std::sqrt(n) / 2.0)) - 1;
    arma::uword k_hi = std::min(n + 0.0, std::ceil((double)k * (double)n / (double)C + t*std::sqrt(n) / 2.0)) - 1;
    
    thetaP_lo = inds(k_lo);
    thetaP_hi = inds(k_hi);
    
    if (thetaP_lo == thetaP_hi)
    {
      termination = 1;
      break;
    }
    
    //count the number of IA :
    // count1 = in (theta_lo, thetaP_lo] = I1
    // count2 = in (thetaP_lo, thetaP_hi) = I2
    // count3 = in [thetaP_hi, thetaP_hi] = I3
    // count4 = in (thetaP_hi, theta_hi] = I4
    permutation = getInterPerm(X, Y, theta_lo, thetaP_lo, 0);
    arma::uword count1 = countInversions(permutation.col(0));
    permutation = getInterPerm(X, Y, thetaP_lo, thetaP_hi, 1); // open interval
    arma::uword count2 = countInversions(permutation.col(0));
    arma::uword count3 = arma::sum(permutation.col(3)) / 2; 
    arma::uword count4 = C - (count1 + count2 + count3);
    
    if (L + count1 >= medind)
    {
      theta_hi = thetaP_lo;
      C = count1;
    }
    else if (L + count1 + count2  >= medind) {
      theta_lo = thetaP_lo;
      theta_hi = thetaP_hi;
      L = L + count1;
      C = count2 + count3;
    }
    else if (L + count1 + count2 + count3 >= medind)
    {
      termination = 2;
      break;
    }
    else
    {
      theta_lo = thetaP_hi;
      L = L + count1 + count2 + count3;
      C = count4;
    }
    
    // check if number of intersection ordinates in the interval 
    // went down over the last 2 iterations
    if ((double) IOcounter(0) / (double)  C < 1.1)
    {
      Rcpp::Rcout << "Warning: algorithm switched early" <<
        " to brute-force enumeration" << std::endl;
      break;
    }
    IOcounter(0) = IOcounter(1);
    IOcounter(1) = C;   
    iterationCounter++;
  }
  
  if (verbose) 
  {
    Rcpp::Rcout << "Interval contraction ended after " << iterationCounter
                << " iterations."<< std::endl;
  }
  
  arma::vec result(2, arma::fill::zeros);
  
  if (termination > 0)
  { // thetaP_lo == thetaP_hi or thetaP_hi degeneracy
    if (verbose) 
    {
      Rcpp::Rcout << "Interval contraction ended in degenarcy."<<
        " No brute force computation needed." << std::endl;
    }
    result(1) = thetaP_hi;
  }
  else
  { // no special termination: brute force computation
    if (verbose) 
    {
      Rcpp::Rcout << "Now starting brute force computation." << std::endl;
    }
    
    // generate a sorted list of all remaining IAs.
    // last element C is a "sentinel value"
    inds = arma::regspace<arma::vec>(0, C);
    inds(C) = C + 1; // sentinel value
    permutation = getInterPerm(X, Y, theta_lo, theta_hi, 0);;
    sampleIA(permutation.col(0), inds, X, Y, permutation.col(1), theta_lo, theta_hi);
    inds = arma::sort(inds);
    result(1) = inds(medind - L - 1);
  }
  
  // determine TS intercept
  arma::vec residuals = Y - result(1) * X;
  std::nth_element(residuals.begin(),
                   residuals.begin() + medind2 - 1,
                   residuals.end());
  
  result(0) = residuals(medind2 - 1); // TS intercept
  
  if (verbose) 
  {
    Rcpp::Rcout << "Algorithm finished" << std::endl;
  }
  return(result); 
}
