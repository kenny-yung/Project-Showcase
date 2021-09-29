//
// Created by Kenny Yung on 5/20/21.
//

#ifndef PROJECT1_P1SOLVERS_H
#define PROJECT1_P1SOLVERS_H
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <utility>
#include "vect.h"
#include "matrix_csr.h"
#include "mat.h"
using namespace std;

vector<vector<double>> read_COO(string filename);
vector<vector<double>> COO2CSR(vector<vector<double>> input);
matrix_csr read_mtx(string filename);
void print_COO(vector<vector<double>> matrix);
void print_CSR(matrix_csr m);
vect GMRES(matrix_csr &A, vect x0, vect b, int m, int &iter, ofstream &file, double restart_r0norm);
vect rGMRES(matrix_csr &A, vect x0, vect b, int m, double tol);
matrix_csr precond_jacobi(matrix_csr &A);
matrix_csr precond_gs(matrix_csr &A);
void print_mat(mat m);
void print_vect(vect v);
vect Cong_Grad(matrix_csr &A, vect &x0, vect b, double tol);

#endif //PROJECT1_P1SOLVERS_H
