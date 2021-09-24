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
  //std::cout << size;
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
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l, int size){
  //int size = sizeof A[0] - 2;
  double t = 0;
  double tau = (A(l,l)-A(k,k))/(2*A(k,l));
  //if(tau == 0){
    //std::cout << "TAU IS 0 /n";
  //}
  if(tau >= 0){
    t = -tau + sqrt(1+pow(tau,2));
  }
  else {
    t = -tau - sqrt(1+pow(tau,2));
  };
  double theta = atan(t);
  double s = sin(theta);
  double c = cos(theta);


  double prev_kk = A(k,k);
  A(k,k) = A(k,k)*pow(c,2) - 2*A(k,l)*c*s + A(l,l)*pow(s,2);
  A(l,l) = A(l,l)*pow(c,2) + 2*A(k,l)*c*s + prev_kk*pow(s,2);
  A(k,l) = 0;
  A(l,k) = 0;
  for(int i=0;i<size; i++){
    if(i != k && i != l){
      double prev_ik = A(i,k);
      //std::cout << "A(i,k)^m" <<prev_ik << "\n";
      A(i,k) = A(i,k)*c - A(i,l)*s;
      //std::cout << "A(i,k)^m after switch"  <<prev_ik << "\n";
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
void jacobi_eigensolver(arma::mat& A, arma::mat& R, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors,
  const int maxiter, int& iterations, bool& converged, int size){
  double max_off_squared = 1.;
  //int size = sizeof A[0]-2;
  std::cout << "Size is " << size << "\n";
  int k = 0;
  int l = 0;


  while(max_off_squared > eps && iterations <= maxiter){
    //std::cout << "This is A before the iteration: " << '\n' << A << '\n';
    double max_off = max_offdiag_symmetric(A, k, l, size);
    //std::cout << "This is max off element A: "<<max_off << '\n';
    //std::cout << "This is k: "<< k << " This is l: " << l << '\n';
    //std::cout << "Hei1";
    jacobi_rotate(A, R, k, l, size);

    max_off_squared = pow(max_offdiag_symmetric(A, k, l, size), 2);
    iterations = iterations + 1;
  }
  if(max_off_squared < eps){
    converged = true;
  }else{
    converged = false;
  }
}

void A_R(arma::mat& A, arma::mat& R, int N){

  double h = 1./(N+1); //Calculate the step length
  double a = -1/pow(h, 2); //Calculate the sup- and super-diagonal elements
  double d = 2/pow(h, 2); //Calculate the main-diagonal elements

  for(int i=0;i<N;i++){
    R(i,i) = 1.;
  }

  //Add the first row
  A(0,0) = d;
  A(0,1) = a;

  //Add the last row
  A(N-1, N-1) = d;
  A(N-1, N-2) = a;

  //Add all other rows
  for(int i = 1; i <= N-2; i++){
    for(int j = 0; j <= N-1; j++){
      if(i == j){
        A(i, j) = d; //Add main diagonal
      }
      else if(abs(i-j) == 1){
        A(i, j) = a;
      }
    }
  }
}

void analytical_solution(arma::mat& v, arma::vec& lambda, int N){
  double h = 1./(N+1); //Calculate the step length
  double a = -1/pow(h, 2); //Calculate the sup- and super-diagonal elements
  double d = 2/pow(h, 2); //Calculate the main-diagonal elements
  for(int i=1; i<=N; i++){
    lambda(i-1) = d + 2*a*cos((i*M_PI)/(N+1));
    for(int j=1; j<=N; j++){
      v(j-1, i-1) = (sin((i*i*M_PI)/(N+1)));
    }
    v.col(i-1) = normalise(v.col(i-1));
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

    //std::cout << max_A_test << '\n';

  //Task 5:
  int n = 7; //Define n number of steps
  double h = 1./n; //Calculate the step length
  double a = -1/pow(h, 2); //Calculate the sup- and super-diagonal elements
  double d = 2/pow(h, 2); //Calculate the main-diagonal elements

  //Defining the matrix A
  int N = n-1; //Defining the size of the matrix
  arma::mat A(N,N, arma::fill::zeros);
  arma::mat R(N,N, arma::fill::zeros);
  A_R(A, R, N);

  //std::cout << "Original matrix A: " << '\n' << A << '\n';


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
  int maxiter = 100000;
  int iterations = 0;
  bool converged = false;


  jacobi_eigensolver(A, R, eps, eigenvalues, eigenvectors, maxiter, iterations, converged, N);
  std::cout << "After Jacobi rotation algorithm A: \n" << A << '\n';
  std::cout << "After Jacobi rotation algorithm R: \n" << R << '\n';

/*
  arma::vec N_l = arma::linspace (5, 120, 24);

  std::ofstream myfile;
  myfile.open("iterations_per_N.txt");

  for(int i=0; i<24;i++){
    int N = N_l(i); //Defining the size of the matrix
    arma::mat A(N,N, arma::fill::zeros);
    arma::mat R(N,N, arma::fill::zeros);
    A_R(A, R, N);
    double eps = pow(10,-8);
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int maxiter = 100000;
    int iterations = 0;
    bool converged = false;
    jacobi_eigensolver(A, R, eps, eigenvalues, eigenvectors, maxiter, iterations, converged, N);
    myfile << std::scientific << iterations << " " << std::scientific << N << "\n";
  }
  */

  //Task 7:


  int N_task_7[2] = {9,99};

  for(int i=0; i<2; i++){
    std::ofstream myfile;
    myfile.open("vec_val_" + std::to_string(N_task_7[i]) + ".txt");
    int N = N_task_7[i]; //Defining the size of the matrix
    arma::mat A(N,N, arma::fill::zeros);
    arma::mat R(N,N, arma::fill::zeros);
    A_R(A, R, N);
    double eps = pow(10,-8);
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int maxiter = 100000;
    int iterations = 0;
    bool converged = false;
    jacobi_eigensolver(A, R, eps, eigenvalues, eigenvectors, maxiter, iterations, converged, N);
    //myfile << std::scientific << iterations << " " << std::scientific << N << "\n";
    arma::vec diagonal = A.diag();

    int idx1 = diagonal.index_min();
    diagonal[idx1] = diagonal.max()+2.2;


    int idx2 = diagonal.index_min();
    diagonal[idx2] = diagonal.max()+2.2;


    int idx3 = diagonal.index_min();
    diagonal[idx3] = diagonal.max()+2.2;

    arma::mat v;
    arma::vec lambda;
    analytical_solution(v, lambda, N);

    double h = 1./(N+1);
    double x_hat[N];
    x_hat[0] = h;
    for(int j=1;j<N;j++){
      x_hat[j] = (j+1) * h;
    }
    myfile << std::scientific << 0 << " " << std::scientific << 0 << " " << std::scientific << 0 <<
     " " << std::scientific << 0 << " " << std::scientific << 0 << " " << std::scientific << 0 << " " << std::scientific << 0 << "\n";

    for(int k=0; k<N; k++){
      myfile << std::scientific << x_hat[k] << " " << std::scientific << R.col(idx1)[k] <<
      " " << std::scientific << R.col(idx2)[k] << " " << std::scientific << R.col(idx3)[k] <<
      std::scientific << v.col(0)[k] << std::scientific << v.col(1)[k] << std::scientific << v.col(2)[k] << std::scientific <<"\n";
    }
    myfile << std::scientific << 1 << " " << std::scientific << 0 << " " << std::scientific << 0 <<
     " " << std::scientific << 0 << " " << std::scientific << 0 << " " << std::scientific << 0 << " " << std::scientific << 0 << "\n";


  }


}
