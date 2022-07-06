# LinearSolver
this repo provides LinearSolver tool to solve Linear System

## run the code
```
g++ test_main.cpp -o test
```
```
./test
```

## test_main.cpp
the main function provides a user iteraction way to test our Solver
it provides two methods of testing:<br/>
**1.create Matrix randomly**<br/>
**2.load Matrix from files**<br/>

seven different solvers are implemented
**Jacobi, Gauss Elimination, GaussSeidel, LuFactorisation, Jacobi sparse solver, LuFactorisation sparse solver, Conjugate Gradient**

**it also provides the error_handler**,
when error occurs in the solver, it returns corresponding error message, 
which is an integer, and the error handler throws error information as follows:
error message number| what kind of error
---|---:
-1|Matrix A is not a square, rows not euqal cols
-2|Matrix A row number not equal size of the vector b
-3|Jacobi cannot converge
-4|LU decomposition first element is zero
-5|Jacobi sparse method cannot converge

you can follow the instructions in the command line,<br/>
enter 1 to create Matrix randomly<br/>
enter 2 to load Matrix from files

### 1.create Matrix randomly
follow the instructions to enter the row number and the column number of the Matrix
our LinearSolver will create Matrix with input row and input col randomly and create vector b
which size equals to the input row number;<br/>

enter 1 to print the Matrix A generated and corresponding CSR Matrix A and vector b<br/>

follow the instructions, you can test each solver (Jacobi, Gauss Elimination, GaussSeidel, LuFactorisation, Jacobi sparse version)<br/>

the x will be written to txt file after<br/>

enter 0 to exit the program.

### 2.load Matrix from files
follow the instructions
press any key to load the default Matrix A and vector b we provide(MatrixA.txt VectorB.txt)
press 0 to input the filename of Matrix A data and vector b data;<br/>

the following steps are the same as create Matrix randomly<br/>

### Write the x to the file
when a solver solves the linear system sucessfully, a file named "x_SolverType.txt"<br/>
will be generated.  e.g if LuFactorisation function is called, no error is raised, when the vector x<br/>
is solved, the "x_LuFactorisation.txt" will be generated and writes the vector x to it<br/><br/><br/>

## LinearSolver class<br/>
### constructor with given Matrix A and vector b
```cpp
LinearSolver(Matrix<T>& A, vector<T>& b);
```
### construct by randomly generate the Matrix A (Arow*Acol) and vector b(Arow)
```cpp
LinearSolver(int Arow, int Acol);
```
### construct by loading Matrix A and vector b
```cpp
LinearSolver(string MatrixAFile, string VectorBfile);
```
### Print Linear System, Print Matrix A, CSR Matrix A, vector b
```cpp
void PrintLinearSystem();
```
### write solver answer, x
```cpp
void WriteSolution(string outfilename);
```
### Gauss Elimination solver, return 1 if solver works successfully, return error code if error occuers
```cpp
int GaussEli();
```
### LU factorisation solver, return 1 if solver works successfully, return error code if error occuers
```cpp
int LuFactorisation();
```
### Jacobi solver, return 1 if solver works successfully, return error code if error occuers
```cpp
  int Jacobi(double tolerance);
```
### Gauss-Seidel solver, return 1 if solver works successfully, return error code if error occuers
```cpp
int GaussSeidel(int iterationNum);
```
### Jacobi sparse version, return 1 if solver works successfully, return error code if error occuers
```cpp
int JacobiCSR(int iterationNum);
```
### LuFactorisationCSR sparse version, return 1 if solver works successfully, return error code if error occuers
```cpp
int LuFactorisationCSR();
```
### Conjugate Gradient, return 1 if solver works successfully, return error code if error occuers
```cpp
int ConjugateGradient(long double tolerance);
```




[![Open in Visual Studio Code](https://classroom.github.com/assets/open-in-vscode-f059dc9a6f8d3a56e377f745f24479a46679e63a5d9fe6f495e02850cd0d8118.svg)](https://classroom.github.com/online_ide?assignment_repo_id=6702442&assignment_repo_type=AssignmentRepo)
