#ifndef ARMA_USE_CXX11
#define ARMA_USE_CXX11
#endif

#ifndef ARMA_DONT_PRINT_ERRORS
#define ARMA_DONT_PRINT_ERRORS
#endif


#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]





// Calculate difference
double calcDiff(double x,double y) 
{
  double dRes = x-y;
  double EPS = 1e-12;
  if(dRes == 0) return 0;
  if(abs(dRes) < EPS * (abs(x) + abs(y)) / 2) return 0;
  return dRes;
}

//////////////////////////////////////////////////////
// Count Individual Inversions and concordant pairs //
//////////////////////////////////////////////////////

void merge_II(arma::vec& y,
              arma::uword left,
              arma::uword middle,
              arma::uword right,
              arma::vec& invCount,
              arma::vec& concCount,
              arma::vec& rankVector)
{
  // c is number of inversions
  arma::uword i, j, k;
  const arma::uword n1 = middle - left + 1;
  const arma::uword n2 = right - middle;
  
  arma::vec Left  = y.subvec(left, left + n1 - 1);
  arma::vec Right = y.subvec(middle + 1, middle + n2);
  
  arma::vec rankVector_Left  = rankVector.subvec(left, left + n1 - 1); 
  arma::vec rankVector_Right = rankVector.subvec(middle + 1, middle + n2); 
  
  arma::uvec tempidx;
  
  i = 0, j = 0, k = left;
  while (i < n1 && j < n2)
  {
    arma::uword ii = i;
    arma::uword jj = j;
    
    if (calcDiff(Left(i),Right(j)) == 0) {
      // count number of duplicates left and right
      while((ii < n1) && (calcDiff(Left(ii),Left(i)) == 0)) {ii++;}
      while((jj < n2) && (calcDiff(Right(jj),Right(j))== 0)) {jj++;}
      int nbLeftdups = ii-i; // number of dups including the first ocurence
      int nbRightdups = jj-j;
      
      invCount(arma::conv_to<arma::uvec>::from(rankVector_Left.subvec(i, i + nbLeftdups - 1)))+=j;
      concCount(arma::conv_to<arma::uvec>::from(rankVector_Right.subvec(j, j + nbRightdups - 1)))+= i;
      
      concCount(arma::conv_to<arma::uvec>::from(rankVector_Left.subvec(i, i + nbLeftdups - 1)))+= n2 - jj;
      invCount(arma::conv_to<arma::uvec>::from(rankVector_Right.subvec(j, j + nbRightdups - 1)))+= n1 - ii;
      
      rankVector.subvec(k, k + nbLeftdups - 1) = rankVector_Left.subvec(i,i + nbLeftdups - 1);
      rankVector.subvec(k + nbLeftdups, k+nbLeftdups+nbRightdups - 1) = rankVector_Right.subvec(j,j + nbRightdups - 1);
      
      arma::vec replacement(nbLeftdups + nbRightdups);
      replacement.fill(Left(i));
      y.subvec(k, k + nbLeftdups + nbRightdups - 1) = replacement;
      k = k + nbLeftdups + nbRightdups;
      i = i + nbLeftdups;
      j = j + nbRightdups;
    }
    else if (Left(i) < Right(j))
    {
      concCount(rankVector_Left(i))+= n2 - j;
      invCount(rankVector_Left(i))+=j;
      
      rankVector(k) = rankVector_Left(i);
      y(k++) = Left(i++);
    }
    else
    {//copy from the right an increase inversion count
      concCount(rankVector_Right(j))+= i;
      invCount(rankVector_Right(j))+= n1 - i;
      
      rankVector(k) = rankVector_Right(j);
      y(k++) = Right(j++);
    }
  }
  
  while (i < n1)
  {
    invCount(rankVector_Left(i))+=j;
    rankVector(k) = rankVector_Left(i);
    y(k++) = Left(i++);
  }
  
  while (j < n2)
  {
    concCount(rankVector_Right(j))+= i;
    rankVector(k) = rankVector_Right(j);
    y(k++) = Right(j++);
  }
}


void mergeSort_II(arma::vec& y, arma::uword left,
                  arma::uword right, arma::vec& invCount,
                  arma::vec& concCount,
                  arma::vec& rankVector)
{
  // 
  if (left < right)
  {
    arma::uword middle = left + (right - left) / 2;
    mergeSort_II(y, left, middle, invCount, concCount, rankVector);
    mergeSort_II(y, middle + 1, right, invCount, concCount, rankVector);
    merge_II(y, left, middle, right, invCount, concCount, rankVector);
  }
}


//[[Rcpp::export]]
arma::mat countIndividualInversions(arma::vec y)
{
  // calculate number of inversions/concordant pairs in vector y
  // for each element in y separately
  
  int n = y.size();
  arma::mat result(n, 2);
  arma::vec invCount(n,arma::fill::zeros);
  arma::vec concCount(n,arma::fill::zeros);
  arma::vec rankVector = arma::regspace(0, n - 1); 
  mergeSort_II(y, 0, n - 1, invCount, concCount, rankVector);
  result.col(0) = invCount;
  result.col(1) = concCount;
  return(result);
}


