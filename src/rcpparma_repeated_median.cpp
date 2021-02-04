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

void merge_RM(arma::uvec& y, int left, int middle, int right, arma::uvec& invCount)
{
  int i, j, k;
  const int n1 = middle - left + 1;
  const int n2 = right - middle;
  
  arma::uvec Left = y.subvec(left, left + n1 - 1);
  arma::uvec Right = y.subvec(middle + 1, middle + n2);
  
  i = 0, j = 0, k = left;
  while (i < n1 && j < n2)
  {
    if (Left(i) <= Right(j))
    {
      invCount(Left(i)) += j;
      y(k++) = Left(i++);
    }
    else
    {
      invCount(Right(j)) += n1 - i; 
      y(k++) = Right(j++);
    }
  }
  
  while (i < n1)
  {
    invCount(Left(i)) += j;
    y(k++) = Left(i++);
  }
  
  while (j < n2)
  {
    y(k++) = Right(j++);
  }
}



void mergeSort_RM(arma::uvec& y, int left, int right, arma::uvec& invCount)
{
  // assumes y is an vector of distinct integers 0, 1, ... n-1 without gaps
  // 
  if (left < right)
  {
    int middle = left + (right - left) / 2;
    mergeSort_RM(y, left, middle, invCount);
    mergeSort_RM(y, middle + 1, right, invCount);
    merge_RM(y, left, middle, right, invCount);
  }
}

arma::uvec countInversions_RM(arma::uvec y,  const arma::uvec lineToval)
{
  // calculate number of inversions in vector y
  // y is assumed to contain distinct integers from 1 to n without gaps
  // algorithm uses a modification of merge sort
  // note that for every element, it counts the number of 
  // inversions in which that element is involved. So they are all counted double!
  
  int n = y.n_elem;
  arma::uvec result(n, arma::fill::zeros);
  mergeSort_RM(y, 0, n - 1, result);
  result = result(lineToval);
  return(result);
}


//////////////////////////////////
// Sample intersection abscissa //
//////////////////////////////////


void merge2_RM(arma::uvec& y, int left, int middle, int right, arma::uvec& counts,
               arma::mat& sampleInds, arma::uvec& indicator, const arma::vec& X,
               const arma::vec& Y, const arma::uvec& valToline,
               const std::unordered_map<int, std::pair<int, int>>& lineKey,
               arma::uvec& IAcounter, double theta_lo, double theta_hi)
{
  int i, j, k;
  const int n1 = middle - left + 1;
  const int n2 = right - middle;
  
  // temporay subvectors
  arma::vec Left = arma::conv_to<arma::vec>::from(y.subvec(left, left + n1 - 1));
  arma::vec Right = arma::conv_to<arma::vec>::from(y.subvec(middle + 1, middle + n2));
  
  i = 0, j = 0, k = left;
  while (i < n1 && j < n2)
  {
    if (Left(i) <= Right(j))
    {
      int countTemp = counts(Left(i)) + j;
      
      if(indicator(Left(i)) == 1)
      { // now we know this line was among the sampled ones
        // now we have to iterate through the duplicates of this line, usually only one
        arma::uword firstInd = lineKey.at(valToline(Left(i))).first;// first
        arma::uword nbDups = lineKey.at(valToline(Left(i))).second;
        for (arma::uword dupNumber = 0; dupNumber < nbDups; dupNumber++)
        {
          while ((IAcounter(firstInd + dupNumber) < sampleInds.n_rows) &&
                 (sampleInds(IAcounter(firstInd + dupNumber),
                             firstInd + dupNumber) < countTemp))
          {
            int jIndex  = sampleInds(IAcounter(firstInd + dupNumber),
                                     firstInd + dupNumber) - counts(Left(i)); // has to be within 0, j
            double intersection;
            if (arma::is_finite(theta_lo) && arma::is_finite(theta_hi))
            {
              if (X(valToline(Right(jIndex))) == X(valToline(Left(i)))) 
              {
                intersection = arma::datum::inf;
              }
              else
              {
                // use more robust formula to compute the sampled intersection ordinate:
                double delta_lo = (theta_lo * X(valToline(Right(jIndex))) -
                                   Y(valToline(Right(jIndex)))) -
                                   (theta_lo * X(valToline(Left(i))) -
                                   Y(valToline(Left(i))));
                double delta_hi = (theta_hi * X(valToline(Right(jIndex))) -
                                   Y(valToline(Right(jIndex)))) -
                                   (theta_hi * X(valToline(Left(i))) -
                                   Y(valToline(Left(i))));
                intersection = theta_lo + delta_lo / (delta_lo - delta_hi) * (theta_hi - theta_lo);
              }
            }
            else
            { // with infinite bounds, use non-robust formula
              
              if (X(valToline(Right(jIndex))) == X(valToline(Left(i)))) 
              {
                intersection = arma::datum::inf;
              }
              else
              {
                intersection =(Y(valToline(Right(jIndex))) - Y(valToline(Left(i)))) /
                  (X(valToline(Right(jIndex))) - X(valToline(Left(i))));
              } 
            }
            sampleInds(IAcounter(firstInd + dupNumber), firstInd + dupNumber) = intersection;
            IAcounter(firstInd + dupNumber) += 1;
          }
        }
      }
      
      counts(Left(i)) = countTemp;
      y(k++) = Left(i++);
    }
    else
    {
      int countTemp = counts(Right(j)) + (n1 - i);
      if(indicator(Right(j)) == 1)
      {// now we know this line was among the sampled ones
        // now we have to iterate through the duplicates of this line, usually only one
        arma::uword firstInd = lineKey.at(valToline(Right(j))).first;// first
        arma::uword nbDups = lineKey.at(valToline(Right(j))).second;
        
        for (arma::uword dupNumber = 0; dupNumber < nbDups; dupNumber++)
        {
          
          while ((IAcounter(firstInd + dupNumber) < sampleInds.n_rows) &&
                 (sampleInds(IAcounter(firstInd + dupNumber),
                             firstInd + dupNumber) < countTemp))
          {
            int iIndex = i + sampleInds(IAcounter(firstInd + dupNumber),
                                        firstInd + dupNumber) - counts(Right(j));
            
            double intersection;
            if (arma::is_finite(theta_lo) && arma::is_finite(theta_hi))
            {
              // use more robust formula to compute the sampled intersection ordinate:
              double delta_lo = (theta_lo * X(valToline(Right(j))) -
                                 Y(valToline(Right(j)))) -
                                 (theta_lo *X(valToline(Left(iIndex))) -
                                 Y(valToline(Left(iIndex))));
              double delta_hi = (theta_hi * X(valToline(Right(j))) -
                                 Y(valToline(Right(j)))) -
                                 (theta_hi * X(valToline(Left(iIndex))) -
                                 Y(valToline(Left(iIndex))));
              intersection = theta_lo + delta_lo / (delta_lo - delta_hi)* (theta_hi - theta_lo);
            }
            else
            { 
              intersection = (Y(valToline(Right(j))) - Y(valToline(Left(iIndex)))) /
                (X(valToline(Right(j))) - X(valToline(Left(iIndex))));
            }
            sampleInds(IAcounter(firstInd + dupNumber), firstInd + dupNumber) = intersection;
            IAcounter(firstInd + dupNumber) += 1;
          }
        }
      }
      counts(Right(j)) = countTemp;
      y(k++) = Right(j++);
    }
  }
  while (i < n1)
  {
    int countTemp = counts(Left(i)) + j;
    if(indicator(Left(i)) == 1)
    { // now we know this line was among the sampled ones
      // now we have to iterate through the duplicates of this line, usually only one
      arma::uword firstInd = lineKey.at(valToline(Left(i))).first;// first
      arma::uword nbDups = lineKey.at(valToline(Left(i))).second;
      for (arma::uword dupNumber = 0; dupNumber < nbDups; dupNumber++)
      {
        while ((IAcounter(firstInd + dupNumber) < sampleInds.n_rows) &&
               (sampleInds(IAcounter(firstInd + dupNumber),
                           firstInd + dupNumber) < countTemp))
        {
          int jIndex  = sampleInds(IAcounter(firstInd + dupNumber),
                                   firstInd + dupNumber) - counts(Left(i)); // has to be within 0, j
          // use more robust formula to compute the sampled intersection ordinate:
          double intersection;
          if (arma::is_finite(theta_lo) && arma::is_finite(theta_hi))
          {
            double delta_lo = (theta_lo * X(valToline(Right(jIndex))) -
                               Y(valToline(Right(jIndex)))) -
                               (theta_lo * X(valToline(Left(i))) -
                               Y(valToline(Left(i))));
            double delta_hi = (theta_hi * X(valToline(Right(jIndex))) -
                               Y(valToline(Right(jIndex)))) -
                               (theta_hi * X(valToline(Left(i))) -
                               Y(valToline(Left(i))));
            intersection = theta_lo + delta_lo / (delta_lo - delta_hi) * (theta_hi - theta_lo);
          }
          else
          {
            intersection =(Y(valToline(Right(jIndex))) - Y(valToline(Left(i)))) /
              (X(valToline(Right(jIndex))) - X(valToline(Left(i))));
          }
          sampleInds(IAcounter(firstInd + dupNumber), firstInd + dupNumber) = intersection;
          IAcounter(firstInd + dupNumber) += 1;
        }
      }
    }
    counts(Left(i)) = countTemp;
    y(k++) = Left(i++);
  }
  while (j < n2)
  {
    y(k++) = Right(j++);
  }
}



void mergeSort2_RM(arma::uvec& y, int left, int right, arma::uvec& counts,
                   arma::mat& sampleInds, arma::uvec& indicator, const arma::vec& X,
                   const arma::vec& Y, const arma::uvec& valToline,
                   const std::unordered_map<int, std::pair<int, int>>& lineKey,
                   arma::uvec& IAcounter, double theta_lo, double theta_hi)
{
  // assumes y is an vector of distinct integers 1, 2, ... n without gaps
  // 
  if (left < right)
  {
    int middle = left + (right - left) / 2;
    mergeSort2_RM(y, left, middle, counts, sampleInds, indicator, X, Y, valToline, lineKey, IAcounter, theta_lo, theta_hi);
    mergeSort2_RM(y, middle + 1, right, counts, sampleInds, indicator, X, Y, valToline, lineKey, IAcounter, theta_lo, theta_hi);
    merge2_RM(y, left, middle, right, counts, sampleInds, indicator, X, Y, valToline, lineKey, IAcounter, theta_lo, theta_hi);
  }
}


void sampleMedIA(arma::uvec y, arma::mat& sampleInds, const arma::vec& X,
                 const arma::vec& Y, const arma::uvec sampledObs,
                 const arma::uvec valToline, const arma::uvec lineToval,
                 const std::unordered_map<int, std::pair<int, int>>& lineKey,
                 double theta_lo, double theta_hi) {
  // sample median intersection abscissa
  const int n = y.n_elem;
  
  arma::uvec indicator(n,arma::fill::zeros);
  indicator(lineToval(sampledObs)).fill(1); // lookup by value
  
  arma::uvec IAcounter(sampleInds.n_cols, arma::fill::zeros); 
  arma::uvec counts(n, arma::fill::zeros);
  mergeSort2_RM(y, 0, n - 1, counts, sampleInds, indicator, X, Y,
                valToline, lineKey, IAcounter, theta_lo, theta_hi);
  
}


/////////////////////
// Repeated Median //
/////////////////////


// [[Rcpp::export]]
arma::vec rcpp_RepeatedMedian(const arma::vec X,
                              const arma::vec Y,
                              const bool verbose,
                              const arma::uword medind,
                              const arma::uvec medind2) {
  // Computes repeated median slope and intercept
  
  // global constants
  const arma::uword n = X.n_elem;
  const int c = 20; // stop contraction iterations when arma::sum(LCRi.col(1)) <= c * n
  const double beta = 0.5;
  const arma::uword r = std::ceil(std::pow(n, beta)); // number of dual lines to sample
  const arma::uword P = std::floor(std::pow(n, 1 - beta)); // number of IA to sample per line
  
  // Initialization
  int termination = 0; 
  double theta_lo = - arma::datum::inf;
  double theta_hi =  arma::datum::inf;
  double thetaP_lo, thetaP_hi;
  arma::uvec L;
  arma::uvec R;
  arma::uvec C = arma::regspace<arma::uvec>(0, n - 1);
  arma::umat LCRi(n, 3, arma::fill::zeros);
  LCRi.col(1).fill(n - 1); 
  arma::mat sampledinds(P, r, arma::fill::zeros);
  arma::uvec inds(r);
  arma::umat permutation(n, 3);
  arma::uvec ranksNeeded(r);
  arma::uvec IOcounter(2);
  IOcounter(1) = arma::sum(LCRi.col(1));
  IOcounter(0) = 2 * IOcounter(1) ;
  int iterationCounter = 0;
  
  if (verbose)
  {
    Rcpp::Rcout << "Initialization finished, starting interval contraction." << std::endl;
  }
  
  while (IOcounter(1) > c * n)
  {
    permutation = getInterPerm(X, Y, theta_lo, theta_hi, 0);
    
    inds = arma::sort(C(arma::randi<arma::uvec>(r, arma::distr_param(0, C.n_elem - 1))));
    
    std::unordered_map<int, std::pair<int, int>> lineKey;
    // map links a key to a pair
    // first element of a pair is the rownumber where the 
    // IA of the correspondling line should be inserted
    // second element is the number of times this line was sampled
    
    // first sample P indices for the first sampled line
    sampledinds.col(0) = arma::sort(arma::randi<arma::vec>(P,
                                    arma::distr_param(0, LCRi(inds(0), 1) - 1))); 
    lineKey.insert({inds(0), std::make_pair(0, 1)});
    
    // now iterate through the other sampled lines 
    // remember the lines were sampled with replacement
    // for each of them, we add a column of sampled 
    // intersection abscissa 
    for (arma::uword colnumber = 1; colnumber < r; colnumber++)
    {
      // sample P indices with replacement between 0 and LCRi(,1) - 1
      // i.e. between 0 and the number of IA of the line under consideration
      // in C.
      sampledinds.col(colnumber) = arma::sort(arma::randi<arma::vec>(P,
                                              arma::distr_param(0,
                                                                LCRi(inds(colnumber),
                                                                     1) - 1))); 
      if(inds(colnumber) != inds(colnumber - 1))
      {// if a new line: make new element in unordered map.
        // key is the line index. the pair denotes the start 
        // of the columns with sampled IA indices connected to this line
        // the second number denotes the number of columns connected
        // to this line
        lineKey.insert({inds(colnumber), std::make_pair(colnumber, 1)});
      }
      else
      {
        int counter = lineKey.at(inds(colnumber)).second + 1; 
        int startnumber = lineKey.at(inds(colnumber)).first;
        lineKey.at(inds(colnumber)) = std::make_pair(startnumber, counter);
      }
    }
    
    sampleMedIA(permutation.col(0), sampledinds, X, Y, inds,
                permutation.col(1), permutation.col(2),
                lineKey, theta_lo, theta_hi);
    arma::uvec leftcounts = LCRi.col(0);
    arma::uvec centercounts = LCRi.col(1);
    
    ranksNeeded = arma::conv_to< arma::uvec >::from(arma::ceil(P * (medind2(inds) -
      arma::conv_to< arma::vec >::from(leftcounts(inds))) / 
      arma::conv_to< arma::vec >::from(centercounts(inds))) - 1);
    arma::vec estimatedMedIAs(r);
    
    for (arma::uword i = 0; i < r; i++)
    {
      if (ranksNeeded(i) > P - 1)
      {
        estimatedMedIAs(i) = theta_hi;
      }
      else if (ranksNeeded(i) < 0)
      {
        estimatedMedIAs(i) = theta_lo;
      }
      else
      {
        arma::vec colvec = sampledinds.col(i);
        std::nth_element(colvec.begin(), colvec.begin() +
          ranksNeeded(i), colvec.end());
        estimatedMedIAs(i) =  colvec(ranksNeeded(i));
      }
    }
    
    // Part c
    arma::uword k = std::ceil(n / 2.0) - L.n_elem;
    arma::uword k_lo =  std::max(1.0, std::floor((r * k) /
      C.n_elem - (3.0 * std::sqrt(r)) / 2.0)) - 1;
    arma::uword k_hi = std::min(r + 0.0, std::floor((r * k) /
      C.n_elem + (3.0 * std::sqrt(r)) / 2.0)) - 1;
    
    std::nth_element(estimatedMedIAs.begin(),
                     estimatedMedIAs.begin() + k_lo,
                     estimatedMedIAs.end());
    thetaP_lo = estimatedMedIAs(k_lo);
    std::nth_element(estimatedMedIAs.begin(),
                     estimatedMedIAs.begin() + k_hi,
                     estimatedMedIAs.end());
    thetaP_hi = estimatedMedIAs(k_hi);
    
    if (thetaP_lo == thetaP_hi)
    {
      termination = 1;
      break;
    }
    
    //Part d
    
    //count the number of IA for every line
    // count1 = in (theta_lo, thetaP_lo] = I1
    // count2 = in (thetaP_lo, thetaP_hi) = I2
    // count3 = in [thetaP_hi, thetaP_hi] = I3
    // count4 = in (thetaP_hi, theta_hi] = I4
    permutation = getInterPerm(X, Y, theta_lo, thetaP_lo, 0);
    arma::uvec count1 = countInversions_RM(permutation.col(0), permutation.col(2));
    permutation = getInterPerm(X, Y, thetaP_lo, thetaP_hi, 1); // open interval
    arma::uvec count2 = countInversions_RM(permutation.col(0), permutation.col(2));
    arma::uvec count3 = permutation.col(3); 
    arma::uvec count4 = LCRi.col(1) - (count1 + count2 + count3);
    
    // Part e  
    // find lines in C for which the median IA lies in I1, I2, I3 and I4.
    
    arma::uvec Lowcounts = LCRi.col(0);
    arma::uvec inds_I1 = arma::find((count1(C) + Lowcounts(C)) >= medind2(C));
    arma::uvec inds_I2 = arma::intersect(arma::find((count1(C) + Lowcounts(C)) < medind2(C)),
                                         arma::find((count1(C) + count2(C) + Lowcounts(C))
                                                      >= medind2(C)));
    arma::uvec inds_I3 = arma::intersect(arma::find((count1(C) + count2(C) +
    Lowcounts(C)) < medind2(C)), arma::find((count1(C) + count2(C) + count3(C) + 
    Lowcounts(C)) >= medind2(C)));
    arma::uvec inds_I4 = arma::find((count1(C) + count2(C) + count3(C) +
      Lowcounts(C)) < medind2(C));
    
    // # interval contraction
    
    if (L.n_elem + inds_I1.n_elem >= medind)
    {
      theta_hi = thetaP_lo;
      // # L doesn't change here
      R = arma::join_cols(arma::join_cols(arma::join_cols(R,
                                                          C(inds_I2)),
                                                          C(inds_I3)),
                                                          C(inds_I4));
      C = C(inds_I1);
      // # update LRCi (now for every line,
      LCRi.col(1) = count1;
      LCRi.col(2) = LCRi.col(2) + count2 + count3 + count4;
    }
    else if (L.n_elem + inds_I1.n_elem + inds_I2.n_elem >= medind)
    { // most likely scenario
      theta_lo = thetaP_lo;
      theta_hi = thetaP_hi;
      L = arma::join_cols(L, C(inds_I1));
      R = arma::join_cols(R, C(inds_I4));
      C = arma::join_cols(C(inds_I2), C(inds_I3)); 
      LCRi.col(0)= LCRi.col(0) + count1;
      LCRi.col(1) = count2 + count3;
      LCRi.col(2) = LCRi.col(2) + count4;
    }
    else if (L.n_elem + inds_I1.n_elem + inds_I2.n_elem + inds_I3.n_elem >= medind)
    {
      termination = 2;
      break;
    }
    else
    {
      theta_lo = thetaP_hi;
      // # R doesn't change here
      L = arma::join_cols(arma::join_cols(arma::join_cols(L,
                                                          C(inds_I1)),
                                                          C(inds_I2)),
                                                          C(inds_I3));
      C = C(inds_I4);
      LCRi.col(0) = LCRi.col(0) + count1 + count2 + count3;
      LCRi.col(1) = count4;
    }
    
    // check if number of intersection ordinates in the interval 
    // went down over the last 2 iterations
    arma::uword IOinInterval = arma::sum(LCRi.col(1));
    if ((double) IOcounter(0) / (double)  IOinInterval < 1.1)
    {
      Rcpp::Rcout << "Warning: algorithm switched early" <<
        " to brute-force enumeration" << std::endl;
      break;
    }
    IOcounter(0) = IOcounter(1);
    IOcounter(1) = IOinInterval;   
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
    permutation = getInterPerm(X, Y, theta_lo, theta_hi, 0);
    arma::uword maxnbbIA = LCRi.col(1).max();
    arma::mat sampledindsfinal(maxnbbIA, C.n_elem, arma::fill::zeros);
    std::unordered_map<int, std::pair<int, int>> lineKey;
    
    for (arma::uword colnumber = 0; colnumber < C.n_elem; colnumber++)
    {
      arma::vec filler(maxnbbIA - LCRi(C(colnumber), 1));
      filler.fill(LCRi(C(colnumber), 1));
      sampledindsfinal.col(colnumber) = arma::join_cols(arma::regspace<arma::vec>(0,
                                                        LCRi(C(colnumber), 1) - 1), filler);
      lineKey.insert({C(colnumber), std::make_pair(colnumber, 1)});
    }
    
    sampleMedIA(permutation.col(0), sampledindsfinal, X, Y, C,
                permutation.col(1), permutation.col(2), lineKey, theta_lo, theta_hi);
    arma::vec medians(C.n_elem);
    
    for (arma::uword colnumber = 0; colnumber < C.n_elem; colnumber++)
    {
      arma::vec candidatemedians = sampledindsfinal.col(colnumber).head(LCRi(C(colnumber), 1));
      
      std::nth_element(candidatemedians.begin(),
                       candidatemedians.begin() +  (medind2(C(colnumber)) - LCRi(C(colnumber), 0) - 1),
                       candidatemedians.end());
      
      medians(colnumber) = candidatemedians(medind2(C(colnumber)) - LCRi(C(colnumber), 0) - 1);
    }
    std::nth_element(medians.begin(),
                     medians.begin() +  (medind - L.n_elem - 1),
                     medians.end());
    result(1) = medians(medind - L.n_elem - 1); // RM slope
  }
  
  // determine RM intercept
  arma::vec residuals = Y - result(1) * X;
  std::nth_element(residuals.begin(),
                   residuals.begin() + medind - 1,
                   residuals.end());
  result(0) = residuals(medind - 1); // rm intercept
  
  if (verbose) 
  {
    Rcpp::Rcout << "Algorithm finished" << std::endl;
  }
  return(result); 
}

