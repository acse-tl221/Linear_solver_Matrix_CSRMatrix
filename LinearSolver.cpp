/*
build our linearsolver class for the problem
*/
#include <iostream>
#include "LinearSolver.h"
#include<algorithm>
#include<math.h>
#include <fstream>
#include <sstream>
#include <string>

//construct our Solver with Matrix A and vector b
template <class T>
LinearSolver<T>::LinearSolver(Matrix<T>& A, vector<T>& b){
    cout<<"call constructor"<<endl;
    this->rows = A.rows;
    this->cols = A.cols;
    this->A = new Matrix<T>(A);
    this->A_CSR = new CSRMatrix<T>(A);
    for(auto m:b){
        this->b.push_back(m);
    }
}

//construct our Solver, Randomly generated Matrix A and Vector b
template <class T>
LinearSolver<T>::LinearSolver(int Arow, int Acol){
    cout<<"***create matrix A and vector b randomly***"<<endl;
	this->rows = Arow;
	this->cols = Acol;
    srand((int)time(0));
    T *A_val = new T[Arow * Acol];
    for(int i =0; i < Arow*Acol; i++){
        A_val[i] = rand()%100;
    }
	Matrix<T> A(Arow,Acol,A_val);
    this->A = new Matrix<T>(A);
	this->A_CSR = new CSRMatrix<T>(A);
    for(int i=0; i<Arow; i++){
        this->b.push_back(rand()%100);
    }
}

//construct our Solver, load Matrix A and vector b from files
template <class T>
LinearSolver<T>::LinearSolver(string MatrixAFile, string VectorBfile){
	cout<<"\n\n*** load matrix A, vector b from files***\n\n";
    ifstream infile(MatrixAFile);
    ifstream infile2(VectorBfile);
    //start loading matrix
    vector <T> mat_val; 
    string temp;
    int row=0, col=0;
    getline(infile,temp);
    //get first line
    stringstream line_input(temp);
    //get the columns
    while(line_input>>temp){
        mat_val.push_back(stof(temp));
        col++;
    }
    row++;
    //iterate from the second line to the end
    while(getline(infile,temp)) 
    {   
        stringstream line_input(temp);
        while(line_input>>temp){
        mat_val.push_back(stof(temp));
        }
        row++;//get the rows
    }
	this->rows = row;
	this->cols = col;
    T *A_val = new T[col*row];
    for (int i=0; i<row*col; i++){
        A_val[i] = mat_val[i];
    }
	Matrix<T> A(this->rows,this->cols,A_val);
    this->A = new Matrix<T>(A);
    this->A_CSR = new CSRMatrix<T>(A);
    vector <T> b;
    getline(infile2,temp);
    stringstream line_input2(temp);
    while(line_input2>>temp){
        b.push_back(stof(temp));
    }
	this->b = b;
}

template <class T>
LinearSolver<T>::~LinearSolver(){
    cout<<"call destructor"<<endl;
    delete A;
    delete A_CSR;
}

template<class T> 
void LinearSolver<T>::PrintLinearSystem(){
	cout<<"print Matrix A"<<endl;
	(this->A)->printMatrix();
    cout<<endl<<"print b"<<endl;
	for(auto it: b){
		cout<<it<<",";
	}
    cout<<"\n\nprint CSRMatrix A"<<endl;
    (this->A_CSR)->printMatrix();
}

template<class T>
void LinearSolver<T>::WriteSolution(string outfilename){
    ofstream outfile(outfilename);
    for (auto it:this->x){
        outfile<<it<<" ";
    }
}

template <class T>
int LinearSolver<T>::Jacobi(int iterationNum){
	if(this->rows!=this->cols){
		return -1;
	}
    else if (this->rows != this->b.size()) {
        return -2;
    }
	else{
		int nrow = this->rows;
		double x_new[nrow];
		double x_old[nrow];
		T *a = (this->A)->values;
        double err;
        double err_previous = 1e10;
		int cnt = 0;
		for (int i=0;i<nrow;i++){
			x_new[i]=0;
			x_old[i]=0;
		}
		while (cnt<iterationNum)
		{
            err = 0;
			for (int i=0;i<nrow;i++)
			{
				x_old[i]=(1/a[i*nrow+i])*(this->b[i]);
				for (int j=0;j<nrow;j++)
				{
					if (j!=i)
					x_old[i]=x_old[i]-(1/a[i*nrow+i])*(a[i*nrow+j] * x_new[j]);
				}
			}
			for (int i=0;i<nrow;i++){
                err += fabs(x_new[i]-x_old[i]);
				x_new[i]=x_old[i];
            }
            if(err>err_previous) return -3;
            err_previous = err;
			cnt++;
		}
        vector <double> x_ans(x_new,x_new+nrow);
		this->x = x_ans;
        return 1;
	}
}

template <class T>
int LinearSolver<T>::GaussSeidel(int iterationNum){
    if(this->rows!=this->cols){
		return -1;
	}
    else if (this->rows != this->b.size()) {
        return -2;
    }
    else{
        int n = this->rows;
        T *a = (this->A)->values;
        double *x_ = new double[this->rows];
        double *y = new double[this->rows];
        int cnt;
        for (cnt=0; cnt<iterationNum; cnt++)
        {
            for (int i = 0; i < n; i++)
                {
                    y[i] = (this->b[i] / a[i*n+i]);
                    for (int j = 0; j < n; j++)
                    {
                        if (j == i) continue;
                        y[i] = y[i] - ((a[i*n+j] / a[i*n+i]) * x_[j]);
                        x_[i] = y[i];
                    }
                }
        }
        vector <double> x_ans(x_,x_+n);
        this->x = x_ans;
        return 1;
    }
}

template <class T>
int LinearSolver<T>::LuFactorisation(){
	int n = this->rows;
    double *l = new double[n*n];//lower triangular
    double *u = new double[n*n];//upper triangle
    double *y = new double[n];//first solve LY=B,then solve UX=Y, y used to store Y
    double *x_ =new double[n];
    T *a = (this->A)->values;
    //decompose matrix into upper triangle matrix and lower triangle matrix
	for(int i =0;i<n*n;i++){
		l[i]=0;
		u[i]=0;
	}
	if(a[0]==0) return -4;
    int i, r, k;
    for (i = 0; i<n; i++)
    {
        u[i] = a[i];
    }
    for (i = 1; i<n; i++)
    {
        l[i*n] = a[i*n] / u[0];
    }

    for (r = 1; r<n; r++)
    {
        for (i = r; i <n; i++)
        {
            double sum = 0;
            for (k = 0; k < r; k++)
            {
                sum += l[r*n+k] * u[k*n+i];// sum of l(r, k) * u(k, i)
            }
            u[r*n+i] = a[r*n+i] - sum;//calculate u(r,i)
        }

        if(r!=n)
            for(i=r+1;i<n;i++)
            {
                double sum = 0;
                for (k = 0; k<r; k++)
                {
                    sum += l[i*n+k] * u[k*n+r];//sum of l(i,k) * u(k,r)
                }
                    l[i*n+r] = (a[i*n+r] - sum) / u[r*n+r];//calculate l(i,r)
            }
    }
    //first solve LY=B
	for(int i_ =0; i_ <n; i_++){
		y[i_] = 0;
	}
    y[0] = b[0];
    for (i = 1; i<n; i++)
    {
        double sum = 0;
        for (k = 0; k<i; k++)
            sum += l[i*n+k] * y[k];
        y[i] = b[i] - sum;
    }
    //then solve UX=Y
    x_[n - 1] = y[n - 1] / u[(n - 1)*n+n - 1];
    for (i = n - 2; i >= 0; i--)
    {
        double sum = 0;
        for (k = i + 1; k<n; k++)
        sum += u[i*n+k] * x_[k];
        x_[i] = (y[i] - sum) / u[i*n+i];
    }
    vector <double> x_ans(x_,x_+n);
	this->x = x_ans;
	return 1;
}

template<class T>
int LinearSolver<T>::GaussEli(){

    // Gauss elimination method
    // Summary: Firstly, transfer the matrix A into upper triangle format
    //          Secondly, back substitute it and get the results

    // upper_triangle function
    // Summary: A function to covert A into upper triangluar form through row operations.
    //            The same row operations are performed on the vector b.
    //            Note that this implementation does not use partial pivoting which is introduced below.

    int n;
    int rows, cols;
    n = this->b.size();
    rows = this->rows;
    cols = this->cols;

    // check whether A is a square(tianchen modified, standard the error message handler)
    if (rows != cols) {
        return -1;
    }
    //and check A has the same number of rows as the size of the vector b
    if (rows != n) {
        cout << "Matrix A does not have the same number as the size of the vector b"<<endl;
        return -2;
    }
    cout<<"Hello A"<<endl;
    // converting operations
    int k, i, j;
    double s;
    double* Avalues_ = new double [this->rows * this->cols];

    vector<double> b_ = this->b;

    // copy the values of matrix A into a new one
    int item;
    for(item = 0; item < rows*cols; item ++){
        Avalues_[item] = this->A->values[item];
    }

    // Loop over each pivot row - all but the last row which we will never need to use as a pivot
    for(k = 0; k < n -1; k++){

    //Loop over each row below the pivot row, including the last row which we do need to update
        for(i = k+1; i < n; i++){
    // Define the scaling factor for this row outside the innermost loop otherwise
    // its value gets changed as you over-write A!!
    // There's also a performance saving from not recomputing things when not strictly necessary
            s = (Avalues_[i * cols + k] / Avalues_[k + k * rows]);
    //Update the current row of A by looping over the column j
    //start the loop from k as we can assume the entries before this are already zero
            for (j = k; j < n; j ++){
                Avalues_[i*cols + j] = Avalues_[i*cols+j] - s*Avalues_[k*cols + j];
            }
            // and update the corresponding entry of b
            b_[i] = b_[i] - s*b_[k];
        }
    }

    int k2,j2;
    double *x_ =new double[n];
    double s2;
    for(k2 = n-1; k2 > -1; k2-- ){
    // note that we could do this update in a single vectorised line
        s2 = 0.0;
        for(j2= k2+1; j2<n; j2++){

            s2 = s2 + Avalues_[k2*cols + j2] * x_[j2];
        }
        x_[k2] = (b_[k2] - s2)/Avalues_[k2*cols + k2];
    }

    vector <double> x_ans(x_,x_+n);
    this->x = x_ans;
    return 1;
}

template<class T>
int LinearSolver<T>::JacobiCSR(int iterationNum){
	if(this->rows!=this->cols){
		return -1;
	}
    else if (this->rows != this->b.size()) {
        return -2;
    }
	else{
		int nrow = this->rows;
		double x_new[nrow];
		double x_old[nrow];
		T *val = (this->A_CSR)->values;
        int *row_pos = (this->A_CSR)->row_position;
        int *col_idx = (this->A_CSR)->col_index;
        double err;
        double err_previous = 1e10;
		int cnt = 0;
		for (int i=0;i<nrow;i++){
			x_new[i]=0;
			x_old[i]=0;
		}
		while (cnt<iterationNum)
		{
            err = 0;
			for (int i=0;i<nrow;i++)
			{
                double tmp1=0;
                for(int k=row_pos[i]; k<row_pos[i+1];k++){
                    if(col_idx[k]==i){
                        tmp1 = val[k];
                    }
                }
				x_old[i]=(1/tmp1)*(this->b[i]);
				for (int j=0;j<nrow;j++)
				{
                    double tmp2 = 0;
                    for(int k=row_pos[i]; k<row_pos[i+1];k++){
                        if(col_idx[k]==j){
                            tmp2 = val[k];
                        }
                    }
					if (j!=i)
					x_old[i]=x_old[i]-(1/tmp1)*(tmp2 * x_new[j]);
				}
			}
			for (int i=0;i<nrow;i++){
                err += fabs(x_new[i]-x_old[i]);
				x_new[i]=x_old[i];
            }
            if(err>err_previous) return -5;
            err_previous = err;
			cnt++;
		}
        vector <double> x_ans(x_new,x_new+nrow);
		this->x = x_ans;
        return 1;
	}
}

template<class T>
int LinearSolver<T>::LuFactorisationCSR(){
	if(this->rows!=this->cols){
		return -1;
	}
    else if (this->rows != this->b.size()) {
        return -2;
    }
    else{
        int n = this->rows;
        double *l = new double[n*n];//lower triangular
        double *u = new double[n*n];//upper triangle
        double *y = new double[n];//first solve LY=B,then solve UX=Y, y used to store Y
        double *x_ =new double[n];
        T *val = (this->A_CSR)->values;//using sparse matrix instead of dense matrix
        int *row_pos = (this->A_CSR)->row_position;
        int *col_idx = (this->A_CSR)->col_index;
        //decompose matrix into upper triangle matrix and lower triangle matrix
        for (int i = 0; i < n; i++)
        {
            // upper triangular
            for (int k = i; k < n; k++)
            {
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (l[i*n+j] * u[j*n+k]);// sum of l(i, j) * u(j, k)
                // calculate u(i,k), get a(i,k) using sparse matrix
                double tmp = 0;
                for(int it=row_pos[i]; it<row_pos[i+1];it++){
                        if(col_idx[it]==k){
                            tmp = val[it];
                        }
                    }
                u[i*n+k] = tmp - sum;
            }
            // lower triangular
            for (int k = i; k < n; k++)
            {
                if (i == k)
                    l[i*n+i] = 1; // set the diagonal elements to 1
                else
                {
                    double sum = 0;
                    for (int j = 0; j < i; j++)
                        sum += (l[k*n+j] * u[j*n+i]);//sum of l(k, j) * u(j, i)
    
                    // calculate l(k, i), get a(k,i) using sparse matrix
                    double tmp2 = 0;
                    for(int it=row_pos[k]; it<row_pos[k+1];it++){
                        if(col_idx[it]==i){
                            tmp2 = val[it];
                        }
                    }
                    l[k*n+i] = (tmp2 - sum) / u[i*n+i];
                }
            }
        }

        for(int i_ =0; i_ <n; i_++){
            y[i_] = 0;
        }
        y[0] = b[0];
        for (int i = 1; i<n; i++)
        {
            double sum = 0;
            for (int k = 0; k<i; k++)
                sum += l[i*n+k] * y[k];
            y[i] = b[i] - sum;
        }
        x_[n - 1] = y[n - 1] / u[(n - 1)*n+n - 1];
        for (int i = n - 2; i >= 0; i--)
        {
            double sum = 0;
            for (int k = i + 1; k<n; k++)
            sum += u[i*n+k] * x_[k];
            x_[i] = (y[i] - sum) / u[i*n+i];
        }
        vector <double> x_ans(x_,x_+n);
        this->x = x_ans;
        return 1;
    }
}

template<class T>
int LinearSolver<T>::ConjugateGradient(long double tolerance){
//Summary: an iterative method which is suitable for sparse matrix.
//Parameter: tolerance indicating the minimum value of residual, which is the stop condition of iteration.

    // initialize size parameters
    int n;
    int rows, cols;
    n = this->b.size();
    rows = this->rows;
    cols = this->cols;
    int* row_p = this->A_CSR->row_position;
    int* cols_i = this->A_CSR->col_index;

    // check whether A is a square(tianchen modified, standard the error message handler)
    if (rows != cols) {
        return -1;
    }
    //and check A has the same number of rows as the size of the vector b
    if (rows != n) {
        cout << "Matrix A does not have the same number as the size of the vector b"<<endl;
        return -2;
    }

    // make a copy of matrix A and b
    double* Avalues_ = new double [this->A_CSR->nnzs];
    int item;
    for(item = 0; item < this->A_CSR->nnzs; item ++){
        Avalues_[item] = this->A_CSR->values[item];
    }
    vector<double> b_ = this->b;

    // 'redefinition' of zero
    const T isNearZero = 1.0e-5;

    // instantiation of iterable x and initially guessed 0
    double* x_tem = new double [n];
    for(int i = 0; i < n; i++)
        x_tem[i] = 0;

    T* r = new T [n];// r
    T* AxPdt = new T[n];// product of A * b
    T AvaluItem;
    int i;
    int j;
    // initialize r
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            AvaluItem = 0;
            for(int ii = row_p[i]; ii < row_p[i+1]; ii ++){
                if(cols_i[ii] == j){
                    // attain the value in ith row jth col
                    AvaluItem = Avalues_[ii];
                }
            AxPdt[i] += AvaluItem * x_tem[j];
            }
        }
        r[i] = b_[i] - AxPdt[i];
    }

    // initialize p = r
    T* p = new T[n];
    for(int i = 0; i<n; i++){
        p[i] = r[i];
    }
    // inner product of transposed r and r
    T rrPdt = innerProduct(r, r, n);
    // product of A and p
    T * ApPdt = new T[n];

    int k = 0; // iteration index
    while (k < n){

        // get the ApPdt in each iteration
        for(int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int ii = row_p[i]; ii < row_p[i + 1]; ii++) {
                    if (cols_i[ii] == j) {
                        AvaluItem = Avalues_[ii];
                    }
                }
                ApPdt[i] += AvaluItem * p[j];
            }
        }

        // get the coefficient alpha
//        T alpha = rrPdt / max(innerProduct(p, ApPdt, n), isNearZero);
        T alpha = rrPdt / innerProduct(p, ApPdt, n);

        // next generation of  x = x + alpha * p
        for(int i = 0; i < n; i++)
            x_tem[i] += alpha * p[i];
        // next generation of residual r = r - alpha * ap;
        for(int i = 0; i < n; i++)
            r[i] -= alpha * ApPdt[i];

        // update the product of r and transposed r
        // convergence condition check with tolerance (argument)
        T new_rrPdt;
        new_rrPdt = innerProduct(r, r, n);
        if(sqrt(new_rrPdt) < tolerance) break;

        // attain the coefficient beta
//        T beta = new_rrPdt / max(rrPdt, isNearZero);
        T beta = new_rrPdt / rrPdt;
        // next generation of p
        for(int i = 0; i < n; i++)
            p[i] = r[i] + beta * p[i];
        //p = r + beta * p;
        rrPdt = new_rrPdt;
        k++;
    }

    vector <double> x_ans(x_tem,x_tem+n);
    this->x = x_ans;
    return 1;
}