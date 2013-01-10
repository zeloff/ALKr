#include "ALKr.h"
using namespace Rcpp;

NumericVector do_pj(NumericMatrix x) {
  return colSums(x) / sum(x);
}

NumericMatrix do_Nk(NumericMatrix Q, NumericMatrix Q_, NumericMatrix A,
           NumericVector fi, NumericVector pj) {
  
  int ni = fi.size();
  int nj = pj.size();
  
  NumericVector rsA = rowSums(A);
  NumericMatrix Qpj(ni, nj);
  for (int j = 0; j < nj; ++j)
    Qpj(_, j) = Q_(_, j) * pj[j];
  
  NumericVector rsQpj = rowSums(Qpj);
  NumericMatrix Nk(ni, nj);
  for (int i = 0; i < ni; ++i) {
    Nk(i, _) = rsA[i] * Q(i, _) + Qpj(i, _) * (fi[i] - rsA[i]) / rsQpj[i];
  }
  return Nk;
}

NumericMatrix do_Nz(NumericMatrix Q_, NumericVector fi, NumericVector pj) {
  int ni = fi.size();
  int nj = pj.size();

  NumericMatrix Qpj(ni, nj);
  for (int j = 0; j < nj; ++j) {
    Qpj(_, j) = Q_(_, j) * pj[j];
  }
  
  NumericVector rsQpj = rowSums(Qpj);
  NumericMatrix Nz(ni, nj);
  for (int i = 0; i < ni; ++i) {
      Nz(i, _) = Qpj(i, _) * fi[i] / rsQpj[i];
  }
  
  return Nz;
}

// [[Rcpp::export]]
List hoenigC(List AAk, List FFik, List FFiz, int threshold, int maxiter) {
  
  NumericMatrix Ak = as<NumericMatrix>(AAk[0]);
  
  int nr = Ak.nrow();
  int nc = Ak.ncol();
  int nlk = AAk.size();
  int nlz = FFiz.size();
  
  NumericMatrix NNk[nlk];
  NumericMatrix QQk[nlk];
  NumericMatrix Q_(nr, nc);
  NumericMatrix sumNk(nr, nc);
  
  for (int l = 0; l < nlk; ++l) {
    NumericMatrix Qk = calc_ALK(as<NumericMatrix>(AAk[l]));
    
    NumericMatrix Nk(nr, nc);
    
    for (int c = 0; c < nc; ++c) {
      Nk(_, c) = Qk(_, c) * as<NumericVector>(FFik[l]);
      sumNk(_, c) = sumNk(_, c) + Nk(_, c);
    }
    
    QQk[l] = Qk;
    NNk[l] = Nk;
  }
  
  NumericVector csNk = colSums(sumNk);
  for (int c = 0; c < nc; ++c) {
    if (csNk[c] == 0) csNk[c] = 1;
    Q_(_, c) = sumNk(_, c) / csNk[c];
  }
  
  NumericVector Pjk[nlk];
  for (int l = 0; l < nlk; ++l) {
    Pjk[l] = do_pj(NNk[l]);
  }
  
  NumericMatrix Nz[nlz];
  NumericMatrix Nz_old[nlz];
  NumericVector Pjz[nlz];
  double sumFiz;
  for (int l = 0; l < nlz; ++l) {
    Nz[l] = NumericMatrix(nr, nc);
    Nz_old[l] = NumericMatrix(nr, nc);
    sumFiz = sum(as<NumericVector>(FFiz[l]));
    for (int r= 0; r < nr; ++r) {
      for (int c = 0; c < nc; ++c) {
        Nz[l](r, c) = Q_(r, c) * (1.0/nc) * sumFiz;
        Nz_old[l](r, c) = Nz[l](r, c);
      }
    }
    Pjz[l] = do_pj(Nz[l]);
  }
  
  int iter = 0;
  double criterion = threshold * 2;
  
  while(criterion > threshold && iter < maxiter) {
    
    NumericMatrix N(nr, nc);
    
    for (int c = 0; c < nc; ++c) {
      for (int l = 0; l < nlk; ++l) {
        N(_, c) = N(_, c) + NNk[l](_, c);
      }
      for (int l = 0; l < nlz; ++l) {
        N(_, c) = N(_, c) + Nz[l](_, c);
      }
    }
   
    NumericVector cSumN = colSums(N);
    NumericMatrix Q_(nr, nc);
    for (int c = 0; c < nc; ++c) {
      if (cSumN[c] == 0) cSumN[c] = 1;
      Q_(_, c) = N(_, c) / cSumN[c];
    }
    for (int l = 0; l < nlk; ++l) {
      NNk[l] = do_Nk(QQk[l], Q_, as<NumericMatrix>(AAk[l]),
                    as<NumericVector>(FFik[l]), Pjk[l]);
      Pjk[l] = do_pj(NNk[l]);
    }
    
    for (int l = 0; l < nlz; ++l) {
      Nz[l] = do_Nz(Q_, as<NumericVector>(FFiz[l]), Pjz[l]);
      Pjz[l] = do_pj(Nz[l]);
      
      criterion = 0;
      for (int c = 0; c < nc; ++c) {
        double this_diff = max(abs(Nz[l](_, c) - Nz_old[l](_, c)));
        if (this_diff > criterion) criterion = this_diff;
        Nz_old[l](_, c) = Nz[l](_, c);
      }
    }
    
    iter += 1;
  }
  
  List result(nlz);
  for (int l = 0; l < nlz; ++l) {
    result[l] = Nz[l];
  }
  
  return List::create(_["Nz"] = result, _["iter"] = iter);
}