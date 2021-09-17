

double max_offdiag_symmetric(arma::mat A, int& k, int &l){
  double max_val = 0;
  int mat_size = len(A);
  for(i=1; i<mat_size; i++){
    for(j=i+1; j<mat_size; j++){
      if(abs(A(i,j)) > max_val){
        max_val = abs(A(i,j);
        k = i;
        l = j;
      }
    }
  }
  return max_val
}

int main()
{
    arma::mat A (4,4, arma::fill::zeros);
    for(i=0;i<4;i++){
      A(i,i) = 1;
    }
    A(3,0) = 0.5;
    A(0,3) = 0.5;
    A(1,2) = -0.7;
    A(2,1) = -0.7;

    double max_A = max_offdiag_symmetric(A, int k, int l);

    std::cout << max_A;
}
