#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
LogicalVector get_unique_ind(DataFrame& df) {
  int n = df.nrow();
  LogicalVector idx(n);
  NumericVector s = df[0];
  NumericVector e = df[1];
  CharacterVector name = df[2];
  for(int i = 0;i<n-5;i++) {
    if(name[i] != name[i+2]) {       // no duplicates
      if(name[i] == name[i+1]) {
        if(s[i] != s[i+1] || e[i] != e[i+1]) {
          idx[i] = true;
          idx[i+1] = true;
        } 
        i++;
      }
    } else {                          // duplicates
      if( name[i] == name[i + 3]) {   // 4 duplicates
        if(s[i] != s[i+1] || e[i] != e[i+1]) {
          idx[i] = true;
          idx[i+1] = true;
        } else if (s[i] != s[i+2] || e[i] != e[i+2]) {
          idx[i] = true;
          idx[i+2] = true;
        } else if (s[i] != s[i+3] || e[i] != e[i+3]) {
          idx[i] = true;
          idx[i+3] = true;
        }
        i+= 3;
      } else {                        // 3 duplicates
        if(s[i] != s[i+1] || e[i] != e[i+1]) {
          idx[i] = true;
          idx[i+1] = true;
        } else if (s[i] != s[i+2] || e[i] != e[i+2]) {
          idx[i] = true;
          idx[i+2] = true;
        }
        i+= 2;
      }
    }
  }
  return idx;
}





// [[Rcpp::export]]
LogicalVector filter_exon_length(DataFrame& df, int max_length, int min_intron_size) {
  int n = df.nrow();
  LogicalVector idx(n);
  NumericVector s = df[0];
  NumericVector e = df[1];
  NumericVector rs = df[5];
  NumericVector re = df[6];
  int max_len = max_length;
  int itron_len = min_intron_size;
  for(int i = 0;i<n-1;i++) {
    if(s[i+1]-1 - (e[i]+1) < max_len & (rs[i+1] - re[i] + 1) <= itron_len) {
      idx[i] = true;
      idx[i+1] = true;
    } 
    i++;
  }
  return idx;
}



// [[Rcpp::export]]
DataFrame pred_exon_coord(DataFrame& df) {
  int n = df.nrow();
  int l = df.nrow()/2;
  CharacterVector seq(l);
  NumericVector le(l);
  NumericVector s(l);
  NumericVector e(l);
  NumericVector rs(l); 
  CharacterVector st(l);
  int j = 0;
  
  CharacterVector seqnames = df["seqnames"];
  NumericVector start = df["start"];
  NumericVector end = df["end"];
  CharacterVector strand = df["strand"];
  
  for(int i = 0;i<n-1;i++) {
    seq[j] = seqnames[i];
    le[j] = start[i] - 1;
    s[j] = end[i] + 1;
    e[j] = start[i+1] - 1;
    rs[j] = end[i+1] + 1;
    st[j] = strand[i];
    j++;
    i++;
  }
  return(Rcpp::DataFrame::create(Rcpp::Named("seqnames")=seq, 
                                 Rcpp::Named("lend")=le, Rcpp::Named("start")=s, 
                                 Rcpp::Named("end")=e, Rcpp::Named("rstart")=rs, 
                                 Rcpp::Named("strand")=st));
}
