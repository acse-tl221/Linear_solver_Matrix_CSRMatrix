#include <iostream>
#include "CSRMatrix.h"

// constructor for csrmatrix
// we are calling the base contructor for matrix here with a false for the preallocate
// that way we don't accidently allocate rows*cols space for the values array
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, bool preallocate): Matrix<T>(rows, cols, false), nnzs(nnzs)
{
   // If we don't pass false in the initialisation list base constructor, it would allocate values to be of size
   // rows * cols in our base matrix class
   // So then we need to set it to the real value we had passed in
   this->preallocated = preallocate;
   // this is set to be rows*cols in our parent class Matrix
   // which is not correct for our CSRMatrix
   this->size_of_values = nnzs;

   // If we want to handle memory ourselves
   if (this->preallocated)
   {
      // Must remember to delete this in the destructor
      this->values = new T[this->nnzs];
      this->row_position = new int[this->rows+1];
      this->col_index = new int[this->nnzs];
   }
}

// Constructor - now just setting the value of our T pointer
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, T *values_ptr, int *row_position, int *col_index): Matrix<T>(rows, cols, values_ptr), nnzs(nnzs), row_position(row_position), col_index(col_index)
{
   // this is set to be rows*cols in our parent class Matrix
   // which is not correct for our CSRMatrix
   this->size_of_values = nnzs;  
}

//constructor - convert a dense matrix to a sparse matrix
template<class T>
CSRMatrix<T>::CSRMatrix(Matrix<T>& dense_mat)
{
   T *arr = dense_mat.values;
   int nnzs = 0;
   T *value; int *row_position; int *col_position;
   int dense_row = dense_mat.rows;
   int dense_cols = dense_mat.cols;
   for (int i = 0; i < dense_row; i++)
        for (int j = 0; j < dense_cols; j++){
           if(arr[i*dense_cols+j]!=0) nnzs++;
        }
   //std::cout<<nnzs<<std::endl;
   value = new T[nnzs];
   row_position = new int[dense_row+1];
   col_index = new int[nnzs];
   int cnt = 0;
   row_position[0] = 0;
   for (int i = 0; i < dense_row; i++){
        for (int j = 0; j < dense_cols; j++){
           if(arr[i*dense_cols+j]!=0) {
              value[cnt] = arr[i*dense_cols+j];
              col_index[cnt] = j;
              cnt++;
           }
        }
      row_position[i+1] = cnt;
   }

   //Matrix<T>(dense_row, dense_cols, value);
   this->rows = dense_row;
   this->cols =dense_cols;
   this->values =value;
   this->nnzs = nnzs;
   this->row_position = row_position;
   this->col_index =col_index;
}

// destructor
template <class T>
CSRMatrix<T>::~CSRMatrix()
{
   // Delete the values array
   if (this->preallocated){
      delete[] this->row_position;
      delete[] this->col_index;
   }
   // The super destructor is called after we finish here
   // This will delete this->values if preallocated is true
}

// Explicitly print out the values in values array as if they are a matrix
template <class T>
void CSRMatrix<T>::printMatrix() 
{ 
   std::cout << "Printing matrix" << std::endl;
   std::cout << "Values: ";
   for (int j = 0; j< this->nnzs; j++)
   {  
      std::cout << this->values[j] << " ";      
   }
   std::cout << std::endl;
   std::cout << "row_position: ";
   for (int j = 0; j< this->rows+1; j++)
   {  
      std::cout << this->row_position[j] << " ";      
   }
   std::cout << std::endl;   
   std::cout << "col_index: ";
   for (int j = 0; j< this->nnzs; j++)
   {  
      std::cout << this->col_index[j] << " ";      
   }
   std::cout << std::endl;   
}

// Do a matrix-vector product
// output = this * input
template<class T>
void CSRMatrix<T>::matVecMult(T *input, T *output)
{
   if (input == nullptr || output == nullptr)
   {
      std::cerr << "Input or output haven't been created" << std::endl;
      return;
   }

   // Set the output to zero
   for (int i = 0; i < this->rows; i++)
   {
      output[i] = 0.0;
   }

   // Loop over each row
   for (int i = 0; i < this->rows; i++)
   {
      // Loop over all the entries in this col
      for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
      {
         // This is an example of indirect addressing
         // Can make it harder for the compiler to vectorise!
         output[i] += this->values[val_index] * input[this->col_index[val_index]];

      }
   }
}

// Blank matmatmult
template<class T>
void CSRMatrix<T>::matMatMult(CSRMatrix<T>& mat_left, CSRMatrix<T>& output)
{
   std::cout << "Inside sparse matmatmult" << std::endl;

}