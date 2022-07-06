/*
build our linearsolver class for the problem
*/
#pragma once
#include <vector>
#include "Matrix.h"
#include "CSRMatrix.h"
using namespace std;
template <class T>
class LinearSolver
{
public:
    //linear solver variable
    Matrix <T> *A;
    CSRMatrix <T> *A_CSR;
    vector <double> x;
    vector <T> b;
    int rows;
    int cols;

    //constructor
    LinearSolver(Matrix<T>& A, vector<T>& b);

    //randomly construct the Matrix A and vector b
    LinearSolver(int Arow, int Acol);

    //constructor load Matrix A and vector b
    LinearSolver(string MatrixAFile, string VectorBfile);

    //destructor
    ~LinearSolver();

    //Print Linear System
    void PrintLinearSystem();

    //write solver answer, x
    void WriteSolution(string outfilename);

    // LU factorisation
    int LuFactorisation();

    // Gauss Elimination
    int GaussEli();

    //Jacobi
    int Jacobi(int iterationNum);

    //Gauss-Seidel
    int GaussSeidel(int iterationNum);

    //Jacobi sparse
    int JacobiCSR(int iterationNum);

    //LU sparse
    int LuFactorisationCSR();

    // ConjugateGradient
    int ConjugateGradient(long double tolerance);
};