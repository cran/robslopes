#include "HelpFunctions.h"



//////////////////////////
// Small help functions //
//////////////////////////

orderRank rank(arma::vec& v)
{
  // "stable" rank vector: ties are "solved" by keeping original ordering
  // returns the ranks and order using 0-based indexing, i.e. in {0, 1, ... n-1}
  // result.orderVector contains order
  // result.rankVector contains ranks
  
  arma::uword n = v.size();
  orderRank result;
  result.orderVector = arma::regspace<arma::uvec>(0, n - 1);
  result.rankVector  = arma::uvec(n, arma::fill::zeros);
  
  std::stable_sort(result.orderVector.begin(), result.orderVector.end(), 
            [&v](int i, int j) { return v(i) < v(j);});
  result.rankVector(result.orderVector) = arma::regspace<arma::uvec>(0, n - 1);
  return result;
}


orderRank rankwTiebreak(arma::vec& v, arma::vec& tieBreaker)
{
  // rank vector with tiebreak: ties are solved by looking at
  // tieBreaker vector: the one with lowest tieBreaker value gets 
  // lowest rank
  // 
  // returns the ranks in {0, 1, ... n-1}! so not 1 to n
  arma::uword n = v.size();
  orderRank result;
  result.orderVector = arma::regspace<arma::uvec>(0, n - 1);
  result.rankVector  = arma::uvec(n, arma::fill::zeros);
  
  std::stable_sort(result.orderVector.begin(), result.orderVector.end(), 
            [&v](int i, int j) { return v(i) < v(j) ;});
  
  // dealing with ties
  for (arma::uword j, i = 0; i < n; i += j)
  {
    j = 1;
    while (i + j < n && v(result.orderVector(i)) ==
           v(result.orderVector(i + j)))
    {
      ++j; // count duplicates
    }
    if (j == 1)
    { // no ties
      result.rankVector(result.orderVector(i)) = i;
    } 
    else if (j > 1) { // Deal with ties using tieBreaker
      arma::uvec tieOrder = result.orderVector.subvec(i, i + j - 1);
      arma::vec subVector = tieBreaker(tieOrder);
      orderRank tieHandling = rank(subVector);
      //insert ranks of the tied obeservations
      result.rankVector(tieOrder) = i + tieHandling.rankVector;
      // permute order of the tied observations by the order of the tiebreaker
      result.orderVector.subvec(i, i + j - 1) = tieOrder.elem(tieHandling.orderVector);
    }
  }
  return result;
}

arma::umat getInterPerm(const arma::vec& X, const arma::vec& Y,
                        double theta_lo, double theta_hi, int open)
{
  // get the permutation (and inverse permutation) of the intersections
  // Handling of ties is done in the following way:
  // Ties at - infinity: by increasing y-intercept
  //         + infinity: by decreasing y-intercept
  //         finite lower bound: by increasing slope
  //         finite higher bound: by increasing slope
  // permutation describes the permutation of the intersections
  // of the lines with theta_lo to theta_hi. 
  // If open is 1, then tie breaking rules for Open interval are used
  // else, tie breaking rules for closed interval are used
  // For the open interval, column 3 of the return matrix 
  // contains the counts of the intersections at theta_hi
  
  // permutation itself takes values in {0, ..., n - 1}
  // valToline is such that given a certain Tag in {0, ..., n -1},
  // valToline[Tag] returns the linenumber in 0,..., n-1
  // line to val is such that given a certain lineNumber in {0, ...n - 1}
  // lineToval[linenumber] returns the Tag in {0, ..., n - 1}
  //   
  orderRank OR_lo, OR_hi;
  arma::umat result;
  arma::vec int_lo, int_hi, TB_lo, TB_hi;
  
  
  if (!arma::is_finite(theta_lo))
  {
    int_lo = -X;
    TB_lo = -Y;
  }
  else
  {
    int_lo = X * theta_lo - Y;
    TB_lo = X;
  }
  OR_lo = rankwTiebreak(int_lo, TB_lo);
  
  if (open == 0)
  {
    result = arma::umat(X.size(), 3);
    if (!arma::is_finite(theta_hi))
    {
      int_hi = X;
      TB_hi = Y;
    }
    else
    {
      int_hi = X * theta_hi - Y;
      TB_hi = X;
    }
    OR_hi = rankwTiebreak(int_hi, TB_hi);
  }
  else
  { //open == 1
    result = arma::umat(X.size(), 4);
    if (!arma::is_finite(theta_hi))
    { 
      int_hi = X;
      TB_hi = -Y;
    }
    else
    {
      int_hi = X * theta_hi - Y;
      TB_hi = -X;
    }
    
    OR_hi = rankwTiebreak(int_hi, TB_hi); 
    
    int_hi = int_hi(OR_hi.orderVector); // sorted intersections with theta_hi
    TB_hi = TB_hi(OR_hi.orderVector); // sorted intersections with theta_hi
    
    
    arma::uvec counts(int_hi.size(), arma::fill::zeros); 
    arma::ivec doubleDupsVec(int_hi.size(), arma::fill::zeros);
    for (arma::uword j, i = 0; i < int_hi.size(); i += j)
    {
      j = 1;
      int doubleDups = 0;
      while (i + j < int_hi.size() && int_hi(i) == int_hi(i + j))
      {
        // careful here: double duplicates (= x-and y-duplicates)
        // are assumed to meet at infinity. However, they 
        // cause equal values for both int_hi and the tie breaker.
        // therefore, they are counted as duplicates meeting at 
        // theta_hi. It also needs to work in a sequence like
        // 11112 with tie breaker 11223 
        // should we also check that theta_hi != Inf? 
        // in that case this correction is probably not needed.
        if (TB_hi(i + j - 1) == TB_hi(i + j)) { // make sure that (i,j) is not x-and y-duplicate
          ++doubleDups;
          
        } else {
          if (doubleDups > 0) {
            doubleDupsVec.subvec(i+j-1 - doubleDups,i+j-1).fill(doubleDups);
          }
          doubleDups = 0;
        }
        ++j; // count duplicates
      }
      if (doubleDups > 0) {
        doubleDupsVec.subvec(i+j-1 - doubleDups,i+j-1).fill(doubleDups);
      }
      if(j > 1) {
        arma::ivec replacement(j);
        replacement.fill(j-1);
        replacement = replacement - doubleDupsVec.subvec(i, i + j - 1);
        counts.subvec(i, i + j - 1) = arma::conv_to<arma::uvec>::from(replacement);
      }
    }
    
    result.col(3) = counts(OR_hi.rankVector); // put the intersection counts in the right place
  }
  result.col(0) = OR_lo.rankVector(OR_hi.orderVector); // permutation 
  result.col(1) = OR_lo.orderVector; // value-to-line key
  result.col(2) = OR_lo.rankVector; // line-to-val key
  return(result);
}
