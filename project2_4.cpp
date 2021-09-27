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
// - Takes as its inputs a symmetric matrix A, two indices k and l, as wellas the number of rows/columns in A.
// - Saves the matrix element indicies to k and l
// - Returns absolute value of A(k,l) as the function return value
double max_offdiag_symmetric(const arma::mat& A, int& k, int& l, int size){
  double max_val = 0; //Storing the maximum value
  for(int i=0; i<size-1; i++){ //Loops through each row from start to the second last
    for(int j=i+1; j<size; j++){ //Loops through all elements to the right of the diagonal elements
      if(abs(A(i,j)) > max_val){ //If the element at index (i,j) is bigger then max_val replace max_val and update k and l
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
// - Takes as its inputs a symmetric matrix A, an identity matrix R with same size as A, two indices k and l, as wellas the number of rows/columns in A.
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l, int size){
  double t = 0; //Storing the t value
  double tau = (A(l,l)-A(k,k))/(2*A(k,l)); //Calculating the tau value
  if(tau >= 0){ //Choosing the t value that gives us the smallest theta.
    t = -tau + sqrt(1+pow(tau,2));
  }
  else {
    t = -tau - sqrt(1+pow(tau,2));
  };
  double theta = atan(t); //Calculate theta
  double s = sin(theta); //Calculate s value
  double c = cos(theta);//Calculate c value


  double prev_kk = A(k,k); //Storing acopy of the A_kk value
  A(k,k) = A(k,k)*pow(c,2) - 2*A(k,l)*c*s + A(l,l)*pow(s,2); //Updating the A_kk value
  A(l,l) = A(l,l)*pow(c,2) + 2*A(k,l)*c*s + prev_kk*pow(s,2); //Updating the A_ll value
  A(k,l) = 0;
  A(l,k) = 0;
  for(int i=0;i<size; i++){ //Looping through the rows of the matrix A
    if(i != k && i != l){//If we arenot at either row k or l, do thefollowing
      double prev_ik = A(i,k); //Store the value A_ik for use later
      A(i,k) = A(i,k)*c - A(i,l)*s; //Updating the A_ik value
      A(k,i) = A(i,k); //Updating the A_ki value
      A(i,l) = A(i,l)*c + prev_ik*s; //Updating the A_il value
      A(l,i) = A(i,l); //Updating the A_li value
    }
  //Rest of function updates the R matrix
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
  double max_off_squared = eps*10; //A variable for storing the maximum of diagonal element squared
  int k = 0; //Stores index k
  int l = 0; //Stores index l


  while(max_off_squared > eps && iterations <= maxiter){ // If the max_off_squared is bigger then our given eps and the number of iterations is less then maxiter, continue. Else stop
    double max_off = max_offdiag_symmetric(A, k, l, size); //Find the maxmimum of diagonal element
    jacobi_rotate(A, R, k, l, size); //Do the rotation part

    max_off_squared = pow(max_offdiag_symmetric(A, k, l, size), 2); //Calculate the new max_off_squared
    iterations = iterations + 1; //Add one iteration to the total number of iterations
  }
  if(max_off_squared < eps){ //Check if we have converged within the given maxiter
    converged = true;
  }else{
    converged = false;
  }
}

void A_R(arma::mat& A, arma::mat& R, int N){

  double h = 1./(N+1); //Calculate the step length
  double a = -1/pow(h, 2); //Calculate the sup- and super-diagonal elements
  double d = 2/pow(h, 2); //Calculate the main-diagonal elements

  for(int i=0;i<N;i++){ //Make the R matrix
    R(i,i) = 1.;
  }

  //Add the first row to A
  A(0,0) = d;
  A(0,1) = a;

  //Add the last row to A
  A(N-1, N-1) = d;
  A(N-1, N-2) = a;

  //Add all other rows to A
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
    lambda(i-1) = d + 2*a*cos((i*M_PI)/(N+1)); //Calculate the analytical eigenvalue values
    for(int j=1; j<=N; j++){
      v(j-1, i-1) = (sin((j*i*M_PI)/(N+1))); //Calculate the analytical eigenvector values
    }
    v.col(i-1) = normalise(v.col(i-1)); //Normalise the columns in v
  }
}



int main()
{
    //Task 4b:
    //Make the matrix we are given to test our function max_offdiag_symmetric
    arma::mat A_test (4,4, arma::fill::zeros);
    for(int i=0;i<4;i++){
      A_test(i,i) = 1; //Add the diagonal
    }
    //Add all other elements to A_test
    A_test(3,0) = 0.5;
    A_test(0,3) = 0.5;
    A_test(1,2) = -0.7;
    A_test(2,1) = -0.7;

    int k; //Initialize the k value
    int l; //Initialize the k value
    double max_A_test = max_offdiag_symmetric(A_test, k, l, 4); //Calculate the max of diagonal element of A_test

    std::cout << "The maximum absolute off diagonal element of the matrix in task 4b is: " <<max_A_test << '\n' << '\n';

  //Task 5:
  int n = 7; //Define n number of steps
  double h = 1./n; //Calculate the step length
  double a = -1/pow(h, 2); //Calculate the sup- and super-diagonal elements
  double d = 2/pow(h, 2); //Calculate the main-diagonal elements

  //Defining the matrix A and R
  int N = n-1; //Defining the size of the matrix
  arma::mat A(N,N, arma::fill::zeros);
  arma::mat R(N,N, arma::fill::zeros);
  A_R(A, R, N);

  //std::cout << "Original matrix A: " << '\n' << A << '\n';


  //Analytical soultion of lambda and v:
  arma::vec lambda(N);
  arma::mat v(N,N);
  analytical_solution(v, lambda, N);


  //Using Jacobi rotation algorithm:
  double eps = pow(10,-8); //Define epsilon
  arma::vec eigenvalues; //Define a vector for storing the eigenvalues
  arma::mat eigenvectors; //Define a matrix for storing the eigenvectors
  int maxiter = 100000; //Define the maximum number of iterations
  int iterations = 0; //Count the number of iterations used
  bool converged = false; //Store if the function has converged or not


  jacobi_eigensolver(A, R, eps, eigenvalues, eigenvectors, maxiter, iterations, converged, N);
  std::cout << "The eigenvalues using the Jacobi algorithm for task 5b is the diagonal elements of this matrix: \n" << A << '\n';
  std::cout << "The eigenvectors using the Jacobi algorithm for task 5b is the columns of this matrix: \n" << R << '\n';

  //Task 6a:
  arma::vec N_l = arma::linspace (5, 120, 24); //Make a vector of different N's to try

  std::ofstream myfile;
  myfile.open("iterations_per_N.txt");

  for(int i=0; i<24;i++){
    int N = N_l(i); //Defining the size of the matrix
    arma::mat A(N,N, arma::fill::zeros);
    arma::mat R(N,N, arma::fill::zeros);
    A_R(A, R, N); //Make the A and R matrix
    double eps = pow(10,-8);
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int maxiter = 100000;
    int iterations = 0;
    bool converged = false;
    jacobi_eigensolver(A, R, eps, eigenvalues, eigenvectors, maxiter, iterations, converged, N);
    myfile << std::scientific << iterations << " " << std::scientific << N << "\n"; //Write the number of iteration needed for convergence to file
  }


  //Task 7:
  int N_task_7[2] = {9,99};
  //Loops through the two values 9 and 99
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

    arma::vec diagonal = A.diag(); //Stores the eigenvalues

    int idx1 = diagonal.index_min(); //Find the index of the smallest eigenvalue
    diagonal[idx1] = diagonal.max()+2.2; //Update the index idx1 in diagonal such that it becomes the biggest element


    int idx2 = diagonal.index_min(); //Find the index of the second smallest eigenvalue
    diagonal[idx2] = diagonal.max()+2.2; //Update the index idx2 in diagonal such that it becomes the biggest element


    int idx3 = diagonal.index_min(); //Find the index of the third smallest eigenvalue
    diagonal[idx3] = diagonal.max()+2.2; //Update the index idx3 in diagonal such that it becomes the biggest element

    //Find the analytical solutions
    arma::mat v(N,N);
    arma::vec lambda(N);
    analytical_solution(v, lambda, N);

    double h = 1./(N+1);
    double x_hat[N]; //Define the x-values
    x_hat[0] = h;

    //Loops through the analytical and estimated solutions and write the values to a file
    for(int j=1;j<N;j++){
      x_hat[j] = (j+1) * h;
    }
    myfile << std::scientific << 0 << " " << std::scientific << 0 << " " << std::scientific << 0 <<
     " " << std::scientific << 0 << " " << std::scientific << 0 << " " << std::scientific << 0 << " "
     << std::scientific << 0 << "\n";

    for(int k=0; k<N; k++){
      myfile << std::scientific << x_hat[k] << " " << std::scientific << R.col(idx1)[k] <<
      " " << std::scientific << R.col(idx2)[k] << " " << std::scientific << R.col(idx3)[k]<< " " <<
      std::scientific << v.col(0)[k] << " " << std::scientific << v.col(1)[k] << " " <<  std::scientific
      << v.col(2)[k] <<"\n";
    }
    myfile << std::scientific << 1 << " " << std::scientific << 0 << " " << std::scientific << 0 <<
     " " << std::scientific << 0 << " " << std::scientific << 0 << " " << std::scientific << 0 << " "
     << std::scientific << 0 << "\n";


  }


}
