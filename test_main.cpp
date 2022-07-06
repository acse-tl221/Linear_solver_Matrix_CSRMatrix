#include <iostream>
#include "LinearSolver.h"
#include "LinearSolver.cpp"
#include "Matrix.cpp"
#include "CSRMatrix.cpp"
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

/*
provide two methods of testing:
create Matrix randomly, load Matrix from files
three solvers are implemented
Jacobi solver, Gauss Elimination, LuFactorisation solver
*/
void error_handler(int err_msg){
    switch (err_msg)
    {
    case -1:
        throw "\nerror:Matrix A is not a square, rows not euqal cols\n";
        break;
    case -2:
        throw "\nerror:Matrix A row number not equal size of the vector b\n";
        break;
    case -3:
        throw "\nerror:Jacobi cannot converge\n";
        break;
    case -4:
        throw "\nerror:LU decomposition first element is zero\n";
        break;
    case -5:
        throw "\nerror:Jacobi sparse method cannot converge\n";
        break;
    default:
        break;
    }
}

int main(){
    int err_message;
    cout<<"welcome to use the linear solver tool\n";
    cout<<"Create Matrix A and Vector B randomly, press 1\n";
    cout<<"Load Matrix A and Vector B from files, press 2\n";
    int mode;
    while(cin>>mode){
        if(mode == 1 || mode ==2) break;
        cout<<"press 1 or press 2\n";
    }

    if(mode == 1 ){
        int Arow, Acol;
        cout<<"please input the Arow and the Acol"<<endl;
        cin>>Arow>>Acol;
        LinearSolver<double> test_solver(Arow,Acol);//initialize class
        cout<<"\ndo you want to show the value of Matrix A and Vector b?\n";
        cout<<"press 1 to show, press other number to skip\n";
        int ifShow;
        cin>>ifShow;
        if(ifShow == 1) test_solver.PrintLinearSystem();
        while(1){
            int solverChoice;
            cout<<"\nchoose one solver\n";
            cout<<"Press 1, choose LuFactorisation solver\n";
            cout<<"Press 2, choose Gauss Elimination solver\n";
            cout<<"Press 3, choose GaussSeidel solver\n";
            cout<<"Press 4, choose Jacobi solver\n";
            cout<<"Press 5, choose LuFactorisation sparse method solver\n";
            cout<<"Press 6, choose Jacobi sparse method solver\n";
            while(cin>>solverChoice){
                if(solverChoice == 1 || solverChoice == 2 || solverChoice ==3 || solverChoice == 4 || solverChoice == 5 || solverChoice == 6) break;
                cout<<"press 1 or 2 or 3 or 4 or 5 or 6 to select the solver\n";
            }
            switch (solverChoice)
            {
            case 1:
                err_message = test_solver.LuFactorisation();
                try{
                    error_handler(err_message);
                }
                catch (const char* msg) {
                    cout << msg << endl;
                }
                if(err_message==1){
                    cout<<endl<<"LuFactorisation solver answer"<<endl<<"printing x"<<endl;
                    test_solver.WriteSolution("x_LuFactorisation.txt");
                    for(auto m:test_solver.x){cout<<m<<",";}
                    cout<<endl;
                }
                break;
            case 2:
                err_message = test_solver.GaussEli();
                try{
                    error_handler(err_message);
                }
                catch (const char* msg) {
                    cout << msg << endl;
                }
                if(err_message==1){
                    cout<<endl<<"Gauss Elimination solver answer"<<endl<<"printing x"<<endl;
                    test_solver.WriteSolution("x_GaussEli.txt");
                    for(auto m:test_solver.x){cout<<m<<",";}cout<<endl;
                }
                break;
            case 3:
                err_message = test_solver.GaussSeidel(1e4);
                try{
                    error_handler(err_message);
                }
                catch (const char* msg) {
                    cout << msg << endl;
                }
                if(err_message==1){
                    cout<<endl<<"GaussSeidel solver answer"<<endl<<"printing x"<<endl;
                    test_solver.WriteSolution("x_GaussSeidel.txt");
                    for(auto m:test_solver.x){cout<<m<<",";}cout<<endl;
                }
                break;
            case 4:
                err_message = test_solver.Jacobi(1e4);
                try {
                    error_handler(err_message);
                }
                catch (const char *msg) {
                    cout << msg << endl;
                }
                if (err_message == 1) {
                    cout << endl << "Jacobi solver answer" << endl << "printing x" << endl;
                    test_solver.WriteSolution("x_Jacobi.txt");
                    for (auto m: test_solver.x) { cout << m << ","; }
                    cout << endl;
                }
                break;
            case 5:
                err_message = test_solver.LuFactorisationCSR();
                try {
                    error_handler(err_message);
                }
                catch (const char *msg) {
                    cout << msg << endl;
                }
                if (err_message == 1) {
                    cout << endl << "LuFactorisation sparse method solver answer" << endl << "printing x" << endl;
                    test_solver.WriteSolution("x_LuFactorisation_sparse.txt");
                    for (auto m: test_solver.x) { cout << m << ","; }
                    cout << endl;
                }
                break;
            case 6:
                err_message = test_solver.JacobiCSR(1e4);
                try {
                    error_handler(err_message);
                }
                catch (const char *msg) {
                    cout << msg << endl;
                }
                if (err_message == 1) {
                    cout << endl << "Jacobi sparse method solver answer" << endl << "printing x" << endl;
                    test_solver.WriteSolution("x_Jacobi_sparse.txt");
                    for (auto m: test_solver.x) { cout << m << ","; }
                    cout << endl;
                }
                break;            
            default:
                break;
            }
            int ifContinue;
            cout<<"\ndo you want to try aother solver, press 0 to exit, press any other number to continue\n";
            cin>>ifContinue;
            if(ifContinue==0) break;
        }
    }

    else{
        string MatrixFile = "MatrixA.txt", VectorFile = "VectorB.txt";
        int ifDefault = 1;
        cout<<"do you want to use the default loading file?\n";
        cout<<"press 0 to input your own input file name, press any other number to use default files\n";
        cin>>ifDefault;
        if(!ifDefault){
            cout<<"\ninput Matrix A loading filename\n";
            cin>>MatrixFile;
            cout<<"input Vector B loading filename\n";
            cin>>VectorFile;
        }
        LinearSolver<double> test_solver(MatrixFile,VectorFile);
        cout<<"\ndo you want to show the value of Matrix A and Vector b?\n";
        cout<<"press 1 to show, press any key to skip\n";
        int ifShow;
        cin>>ifShow;
        if(ifShow == 1) test_solver.PrintLinearSystem();
        while(1){
            int solverChoice;
            cout<<"\nchoose one solver\n";
            cout<<"Press 1, choose LuFactorisation solver\n";
            cout<<"Press 2, choose Gauss Elimination solver\n";
            cout<<"Press 3, choose GaussSeidel solver\n";
            cout<<"Press 4, choose Jacobi solver\n";
            cout<<"Press 5, choose LuFactorisation sparse method solver\n";
            cout<<"Press 6, choose Jacobi sparse method solver\n";
            while(cin>>solverChoice){
                if(solverChoice == 1 || solverChoice == 2 || solverChoice ==3 || solverChoice == 4 || solverChoice == 5 || solverChoice == 6) break;
                cout<<"press 1 or 2 or 3 or 4 or 5 or 6 to select the solver\n";
            }
            switch (solverChoice)
            {
            case 1:
                err_message = test_solver.LuFactorisation();
                try{
                    error_handler(err_message);
                }
                catch (const char* msg) {
                    cout << msg << endl;
                }
                if(err_message==1){
                    cout<<endl<<"LuFactorisation solver answerJacobi solver answer"<<endl<<"printing x"<<endl;
                    test_solver.WriteSolution("x_LuFactorisation.txt");
                    for(auto m:test_solver.x){cout<<m<<",";}
                    cout<<endl;
                }
                break;
            case 2:
                err_message = test_solver.GaussEli();
                try{
                    error_handler(err_message);
                }
                catch (const char* msg) {
                    cout << msg << endl;
                }
                if(err_message==1){
                    cout<<endl<<"Gauss Elimination solver answer"<<endl<<"printing x"<<endl;
                    test_solver.WriteSolution("x_GaussEli.txt");
                    for(auto m:test_solver.x){cout<<m<<",";}cout<<endl;
                }
                break;
            case 3:
                err_message = test_solver.GaussSeidel(1e4);
                try{
                    error_handler(err_message);
                }
                catch (const char* msg) {
                    cout << msg << endl;
                }
                if(err_message==1){
                    cout<<endl<<"GaussSeidel solver answer"<<endl<<"printing x"<<endl;
                    test_solver.WriteSolution("x_GaussSeidel.txt");
                    for(auto m:test_solver.x){cout<<m<<",";}cout<<endl;
                }
                break;
            case 4:
                err_message = test_solver.Jacobi(1e4);
                try {
                    error_handler(err_message);
                }
                catch (const char *msg) {
                    cout << msg << endl;
                }
                if (err_message == 1) {
                    cout << endl << "Jacobi solver answer" << endl << "printing x" << endl;
                    test_solver.WriteSolution("x_Jacobi.txt");
                    for (auto m: test_solver.x) { cout << m << ","; }
                    cout << endl;
                }
                break;
            case 5:
                err_message = test_solver.LuFactorisationCSR();
                try {
                    error_handler(err_message);
                }
                catch (const char *msg) {
                    cout << msg << endl;
                }
                if (err_message == 1) {
                    cout << endl << "LuFactorisation sparse method solver answer" << endl << "printing x" << endl;
                    test_solver.WriteSolution("x_LuFactorisation_sparse.txt");
                    for (auto m: test_solver.x) { cout << m << ","; }
                    cout << endl;
                }
                break;
            case 6:
                err_message = test_solver.JacobiCSR(1e4);
                try {
                    error_handler(err_message);
                }
                catch (const char *msg) {
                    cout << msg << endl;
                }
                if (err_message == 1) {
                    cout << endl << "Jacobi sparse method solver answer" << endl << "printing x" << endl;
                    test_solver.WriteSolution("x_Jacobi_sparse.txt");
                    for (auto m: test_solver.x) { cout << m << ","; }
                    cout << endl;
                }
                break;  
            default:
                break;
            }
            int ifContinue;
            cout<<"\ndo you want to try aother solver, press 0 to exit, press any other number to continue\n";
            cin>>ifContinue;
            if(ifContinue==0) break;
        }
    }
}
