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

void merge_PB(arma::uvec& y, arma::uword left, arma::uword middle, arma::uword right, arma::uword& invCount)
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


void mergeSort_PB(arma::uvec& y, arma::uword left, arma::uword right, arma::uword& invCount)
{
  // assumes y is an vector of distinct integers 0, 1, ... n-1 without gaps
  // 
  if (left < right)
  {
    arma::uword middle = left + (right - left) / 2;
    mergeSort_PB(y, left, middle, invCount);
    mergeSort_PB(y, middle + 1, right, invCount);
    merge_PB(y, left, middle, right, invCount);
  }
}


arma::uword countInversions_PB(arma::uvec y)
{
  // calculate number of inversions in vector y
  // y is assumed to contain distinct integers from 1 to n without gaps
  // algorithm uses a modification of merge sort
  
  int n = y.size();
  arma::uword result = 0;
  mergeSort_PB(y, 0, n - 1, result);
  return(result);
}




//////////////////////////////////
// Sample intersection abscissa //
//////////////////////////////////

void merge2_PB(arma::uvec& y,arma::uword left, arma::uword middle, arma::uword right,
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



void mergeSort2_PB(arma::uvec& y, arma::uword left, arma::uword right, arma::vec& sampleInds,
                   const arma::vec& X, const arma::vec& Y, arma::uword& invCount,
                   arma::uword& sampleInds_idx, const arma::uvec& valToline,
                   double theta_lo, double theta_hi)
{
  // assumes y is an vector of distinct integers 0, 1, ... n-1 without gaps
  // 
  if (left < right)
  {
    arma::uword middle = left + (right - left) / 2;
    mergeSort2_PB(y, left, middle, sampleInds, X, Y, invCount,
                  sampleInds_idx, valToline, theta_lo, theta_hi);
    mergeSort2_PB(y, middle + 1, right, sampleInds, X, Y, invCount,
                  sampleInds_idx, valToline, theta_lo, theta_hi);
    merge2_PB(y, left, middle, right, sampleInds, X, Y, invCount,
              sampleInds_idx, valToline, theta_lo, theta_hi);
  }
}


void sampleIA_PB(arma::uvec y, arma::vec& sampleInds, 
              const arma::vec& X, const arma::vec& Y,
              const arma::uvec valToline,
              double theta_lo, double theta_hi)
{
  // sampleInds contains the slopes after 
  int n = y.size();
  arma::uword invCount = 0;
  arma::uword sampleInds_idx = 0;
  mergeSort2_PB(y, 0, n - 1, sampleInds, X, Y, invCount, 
                sampleInds_idx, valToline, theta_lo, theta_hi);
}


////////////////////
// Passing Bablok //
////////////////////


//[[Rcpp::export]]
arma::vec rcpp_PassingBablok(const arma::vec X,
                             const arma::vec Y,
                             const bool verbose,
                             const arma::uword medind0,
                             const arma::uword medind1) {
  // Computes median absolute slope and intercept
  
  // global constants
  const arma::uword n = X.size();
  const arma::uword c = 20; // stop contraction iterations when C<= c * n
  const arma::uword t = 3; // tuning parameter used in k_lo and k_hi
  const arma::uword r = n;
  // Initialization
  int termination = 0;
  double theta_lo = 0;
  double theta_hi =  arma::datum::inf;
  double thetaP_lo, thetaP_hi = 0;
  arma::uword C = (arma::uword)(n * (n - 1) / 2);
  // C = number of absolute ordinates in (theta_lo,theta_hi]
  
  arma::uvec IOcounter(2); // used to check if C goes still goes down
  IOcounter(1) = C;
  IOcounter(0) = 2 * IOcounter(1);
  arma::umat permutation_l(n, 3);
  arma::umat permutation_h(n, 3);
  int iterationCounter = 0;
  
  // initialize counts
  arma::uword C_m = 0; // number of absolute slopes in [0,\theta_lo]
  // should be initiated as sum_j choose(nbtimes_j, 2) where j ranges over Y duplicates
  arma::uword C_l = 0;// number of slopes in [-\theta_hi,-\theta_lo)
  arma::uword C_h = 0;// number of slopes in (\theta_lo,\theta_hi]
  
  permutation_l = getInterPerm(X, Y, -arma::datum::inf,
                               0, 1); // open interval
  C_l = countInversions_PB(permutation_l.col(0));
  C_m = arma::sum(permutation_l.col(3)) / 2; 
  C_h = C - C_l - C_m;
  C   = C_h + C_l;
  // C will get updated each iteration as C <- C_h + C_l

  if (C_m > medind1) {// degenerate case in which estimated slope is 0
    termination  = 3;
  }
  
  arma::uword r_l = std::floor((double)r * ((double)C_l) / ((double)C_l + (double)C_h)) ; // number of lower random ordinates to sample each iteration
  arma::uword r_h = std::ceil((double)r * ((double)C_h) / ((double)C_l + (double)C_h)) ; // number of lower random ordinates to sample each iteration
  
  arma::vec inds(r);
  arma::vec inds_l;
  arma::vec inds_h;
  
  
  if (verbose)
  {
    Rcpp::Rcout << "Initialization finished, starting interval contraction." << std::endl;
  }
  
  while ((C > (arma::uword)(c * n)) && termination == 0) //switch to brute force otherwise
  {
    
    
    // sample ordinates
    
    // now we need to generate 2 sorted lists of n random numbers in total
    // with replacement
    // we avoid arma::randi, since it has maximum_integer as upper bound
    inds_l = arma::floor(arma::sort(arma::randu<arma::vec>(r_l + 1))* double(C_l));
    inds_h = arma::floor(arma::sort(arma::randu<arma::vec>(r_h + 1))* double(C_h));
    inds_l(r_l) = (arma::uword)(n * (n - 1) / 2) + 1; // sentinel value
    inds_h(r_h) = (arma::uword)(n * (n - 1) / 2) + 1; // sentinel value
    
    permutation_l = getInterPerm(X, -Y, theta_lo, theta_hi, 0);
    
    permutation_h = getInterPerm(X, Y, theta_lo, theta_hi, 0);
    
    sampleIA_PB(permutation_l.col(0), inds_l, X, -Y,
             permutation_l.col(1), theta_lo, theta_hi);
    
    sampleIA_PB(permutation_h.col(0), inds_h, X, Y,
             permutation_h.col(1),theta_lo, theta_hi);
    
    inds.head(r_l) = inds_l.head(r_l); // remove sentinel values
    inds.tail(r_h) = inds_h.head(r_h); // remove sentinel values
    inds = arma::sort(arma::abs(inds)); // we take the absolute value here
    
    arma::uword k = medind1 - C_m;
    
    arma::uword k_lo = std::max(1.0, std::floor((double)k * (double)n /
      (double)C - t*std::sqrt(n) / 2.0)) - 1;
    arma::uword k_hi = std::min(n + 0.0, std::ceil((double)k * (double)n /
      (double)C + t*std::sqrt(n) / 2.0)) - 1;
    
    // candidate bounds:
    thetaP_lo = inds(k_lo);
    thetaP_hi = inds(k_hi);
    
    if (thetaP_lo == thetaP_hi)
    {// Candidate contraction bounds form a singleton.
      if (theta_lo != thetaP_lo) {
        thetaP_lo = theta_lo;
      } else if (theta_hi != thetaP_hi) {
        thetaP_hi = theta_hi;
      } else {
        termination  = 2;
        break;
      }
    }
    
    //count the number of IA :
    // count1 = in [-\lb', -\lb) \cup (\lb, \lb'] 
    // count2 = in (-\ub', -\lb') \cup (\lb', \ub')
    // count3 = in {-\ub', \ub}
    // count4 = in (-\ub, -\ub') \cup (\ub', \ub)
    
    
    permutation_l = getInterPerm(X, -Y, theta_lo, thetaP_lo, 0);
    permutation_h = getInterPerm(X, Y, theta_lo, thetaP_lo, 0);
    arma::uword count1_l = countInversions_PB(permutation_l.col(0));
    arma::uword count1_h = countInversions_PB(permutation_h.col(0));
    arma::uword count1 = count1_l + count1_h;
    
    permutation_l = getInterPerm(X, -Y, thetaP_lo, thetaP_hi, 1);
    permutation_h = getInterPerm(X, Y, thetaP_lo, thetaP_hi, 1);
    arma::uword count2_l = countInversions_PB(permutation_l.col(0));
    arma::uword count2_h = countInversions_PB(permutation_h.col(0));
    arma::uword count2 = count2_l + count2_h;
    
    arma::uword count3_l = 0;
    if (thetaP_hi != arma::datum::inf) { 
      count3_l = arma::sum(permutation_l.col(3)) / 2;
    }
    arma::uword count3_h = arma::sum(permutation_h.col(3)) / 2;
    arma::uword count3 = count3_l + count3_h;
    
    
    arma::uword count4_l = C_l - (count1_l + count2_l + count3_l);
    arma::uword count4_h = C_h - (count1_h + count2_h + count3_h);
    // arma::uword count4 = count4_l + count4_h;
    
    
    if (C_m + count1 >= medind1)
    {
      theta_hi = thetaP_lo;
      C_l = count1_l;
      C_h = count1_h;
    }
    else if (C_m + count1 + count2  >= medind1) {
      theta_lo = thetaP_lo;
      theta_hi = thetaP_hi;
      C_m = C_m + count1;
      C_l = count2_l + count3_l;
      C_h = count2_h + count3_h;
    }
    else if (C_m + count1 + count2 + count3 >= medind1)
    {
      termination = 2;
      break;
    }
    else
    {
      theta_lo = thetaP_hi;
      C_m = C_m + count1 + count2 + count3;
      C_l = count4_l;
      C_h = count4_h;
    }
    // updates after updating C_h and C_l
    C = C_h + C_l;
    r_l = std::floor((double)r * ((double)C_l) /
      ((double)C_l + (double)C_h)) ; // number of lower random ordinates to sample each iteration
    r_h = std::ceil((double)r * ((double)C_h) /
      ((double)C_l + (double)C_h)) ; // number of upper random ordinates to sample each iteration
    
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
  
  if (termination == 2)
  { // thetaP_lo == thetaP_hi or thetaP_hi degeneracy
    if (verbose) 
    {
      Rcpp::Rcout << "Interval contraction ended in degeneracy."<<
        " No brute-force computation needed." << std::endl;
    }
    result(1) = thetaP_hi;
  }
  else if (termination == 3) {
    if (verbose) 
    {
      Rcpp::Rcout << "Careful, many zero-slopes in the sample." << std::endl;
    }
    result(1) = 0;
  }
  else
  { // no special termination: brute force computation
    if (verbose) 
    {
      Rcpp::Rcout << "Now starting brute-force computation." << std::endl;
    }
    
    // generate a sorted list of all remaining IAs.
    // last element C is a "sentinel value"
    inds_l = arma::regspace<arma::vec>(0, C_l);
    inds_h = arma::regspace<arma::vec>(0, C_h);
    inds_l(C_l) = (arma::uword)(n * (n - 1) / 2) + 1;
    inds_h(C_h) = (arma::uword)(n * (n - 1) / 2) + 1;
    
    if (theta_hi == arma::datum::inf) { // do not count duplicates in this case
      permutation_l = getInterPerm(X, -Y, theta_lo, theta_hi, 1);
    } else {
      permutation_l = getInterPerm(X, -Y, theta_lo, theta_hi, 0);
    }
    
    permutation_h = getInterPerm(X, Y, theta_lo, theta_hi, 0);
    
    sampleIA_PB(permutation_l.col(0), inds_l, X, -Y,
             permutation_l.col(1), theta_lo, theta_hi);
    sampleIA_PB(permutation_h.col(0), inds_h, X, Y,
             permutation_h.col(1), theta_lo, theta_hi);
    
  
    inds_l = inds_l.head(C_l); // remove sentinel values
    inds_h = inds_h.head(C_h); // remove sentinel values
    inds = arma::join_cols(inds_l, inds_h);
    inds = arma::sort(arma::abs(inds)); // we take the absolute value here
    
    result(1) = inds(medind1 - C_m - 1);
  }
  
  // determine PB intercept
  arma::vec residuals = Y - result(1) * X;
  std::nth_element(residuals.begin(),
                   residuals.begin() + medind0 - 1,
                   residuals.end());
  
  result(0) = residuals(medind0 - 1); // PB intercept
  
  if (verbose) 
  {
    Rcpp::Rcout << "Algorithm finished" << std::endl;
  }
  return(result); 
}



