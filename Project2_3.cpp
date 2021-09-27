#include <armadillo>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <chrono>
#include <cmath>


int main()
{

  int n = 7; //Define n number of steps
  double h = 1./n; //Calculate the step length
  double a = -1/pow(h, 2); //Calculate the sup- and super-diagonal elements
  double d = 2/pow(h, 2); //Calculate the main-diagonal elements

  //Defining the matrix A
  int N = n-1; //Defining the size of the matrix
  arma::mat A (N,N, arma::fill::zeros);

  //Add the first row
  A(0,0) = d;
  A(0,1) = a;

  //Add the last row
  A(N-1, N-1) = d;
  A(N-1, N-2) = a;

  //Add all other rows
  for(int i = 1; i <= 4; i++){
    for(int j = 0; j <= 5; j++){
      if(i == j){
        A(i, j) = d; //Add main diagonal
      }
      else if(abs(i-j) == 1){
        A(i, j) = a;
      }
    }
  }
  std::cout << '\n' << "This is the original matrix: " << '\n' <<  A << '\n';


  //Analytical soultion of lambda and v:
  arma::vec lambda(N);
  arma::mat v(N, N);

  for(int i=1; i<=N; i++){
    lambda(i-1) = d + 2*a*cos((i*M_PI)/(N+1)); //Calculate the analytical eigenvalue values
    for(int j=1; j<=N; j++){
      v(j-1, i-1) = (sin((j*i*M_PI)/(N+1))); //Calculate the analytical eigenvector values
    }
    v.col(i-1) = normalise(v.col(i-1)); //Normalise the columns in v
  }


  std::cout << "This is the analytical solution of the eigenvalues: " << "\n" <<lambda << "\n";
  std::cout << "This is the analytical solution of the eigenvectors: " << "\n "<< v << "\n";

  //armadillos soultion
  arma::vec eigval;
  arma::mat eigvec;

  eig_sym(eigval, eigvec, A); //Solve it using armadillo

  std::cout << "This is the armadillo solution of the eigenvalues: " << "\n" <<eigval << "\n";
  std::cout << "This is the armadillo solution of the eigenvectors: " << "\n "<< eigvec << "\n";


  return 0;
}
