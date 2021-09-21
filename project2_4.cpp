#include <armadillo>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <chrono>
#include <cmath>


using namespace std;

// Determine the the max off-diagonal element of a symmetric matrix A
// - Saves the matrix element indicies to k and l
// - Returns absolute value of A(k,l) as the function return value
double max_offdiag_symmetric(const arma::mat& A, int& k, int& l, int size){
  double max_val = 0;
  std::cout << size;
  for(int i=0; i<size-1; i++){
    for(int j=i+1; j<size; j++){
      if(abs(A(i,j)) > max_val){
        max_val = abs(A(i,j));
        k = i;
        l = j;
      }
    }
  }
  return max_val;
}



// Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l){
  int size = sizeof A[0];
  double tau = 0;
  double t = 0;
  double theta = 0;
  tau = (A(l,l)-A(k,k))/(2*A(k,l));
  t = min(-tau + sqrt(1+pow(tau,2)), -tau - sqrt(1+pow(tau,2)));
  theta = atan(t);
  double s = sin(theta);
  double c = cos(theta);


  A(k,k) = A(k,k)*pow(c,2) - 2*A(k,l)*c*s + A(l,l)*pow(s,2);
  A(l,l) = A(l,l)*pow(c,2) + 2*A(k,l)*c*s + A(k,k)*pow(s,2);
  A(k,l) = 0;
  A(l,k) = 0;
  for(int i=0;i<size; i++){
    if(i != k && i != l){
      double prev_ik = A(i,k);
      A(i,k) = A(i,k)*c - A(i,l)*s;
      A(k,i) = A(i,k);
      A(i,l) = A(i,l)*c + prev_ik*s;
      A(l,i) = A(i,l);
    }
  double prev_ik_R = R(i,k);
  R(i,k) = R(i,k)*c - R(i,l)*s;
  R(i,l) = R(i,l)*c + prev_ik_R*s;

  }
}

// Jacobi method eigensolver:
// - Runs jacobo_rotate until max off-diagonal element < eps
// - Writes the eigenvalues as entries in the vector "eigenvalues"
// - Writes the eigenvectors as columns in the matrix "eigenvectors"
//   (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
// - Stops if it the number of iterations reaches "maxiter"
// - Writes the number of iterations to the integer "iterations"
// - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter
void jacobi_eigensolver(arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, const int maxiter, int& iterations, bool& converged){
  double max_off_squared = 1.;
  int size = sizeof A[0]-2;
  std::cout << size;
  int k = 0;
  int l = 0;
  arma::mat R(size,size, arma::fill::zeros);
  for(int i=0;i<size;i++){
    R(i,i) = 1.;
  }

  while(max_off_squared > eps && iterations <= maxiter){
    std::cout << "This is A before the iteration: " << '\n' << A << '\n';
    double max_off = max_offdiag_symmetric(A, k, l, size);
    std::cout << "This is max off element A: "<<max_off << '\n';

    jacobi_rotate(A, R, k, l);

    max_off_squared = pow(max_offdiag_symmetric(A, k, l, size), 2);
    iterations = iterations + 1;
  }
  if(max_off_squared < eps){
    converged = true;
  }else{
    converged = false;
  }
}

// A main function for task 4b
int main()
{
    arma::mat A_test (4,4, arma::fill::zeros);
    for(int i=0;i<4;i++){
      A_test(i,i) = 1;
    }
    A_test(3,0) = 0.5;
    A_test(0,3) = 0.5;
    A_test(1,2) = -0.7;
    A_test(2,1) = -0.7;

    int k;
    int l;
    double max_A_test = max_offdiag_symmetric(A_test, k, l, 4);

    std::cout << max_A_test << '\n';

  //Task 5:
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
  std::cout << "Original matrix A: " << '\n' << A << '\n';


  //Analytical soultion of lambda and v:
  arma::vec lambda(N);
  arma::vec v(N);
  for(int i=1; i<=N; i++){
    lambda(i-1) = d + 2*a*cos((i*M_PI)/(N+1));
    v(i-1) = (sin((i*i*M_PI)/(N+1)));
  }
  arma::vec v_norm = normalise(v);


  //Using Jacobi rotation algorithm:
  double eps = pow(10,-8);
  arma::vec eigenvalues;
  arma::mat eigenvectors;
  int maxiter = 2;
  int iterations = 0;
  bool converged = false;

  jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);
  std::cout << "After Jacobi rotation algorithm: " << A << '\n';

}
