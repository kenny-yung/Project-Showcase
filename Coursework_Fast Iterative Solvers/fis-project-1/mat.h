//
// Created by Kenny Yung on 6/6/21.
//

#ifndef PROJECT1_MAT_H
#define PROJECT1_MAT_H
#include <iostream>
#include <vector>
#include <cmath>
#include "vect.h"
using namespace std;

class mat {
public:
    mat(int nrow, int ncol, double val);
    mat(vector<vector<double>> v);
    int m_nrow, m_ncol;
    vector<vector<double>> m_matrix;
    double& operator()(int i, int j);
//    double operator()(int i, int j) const;
    mat transpose();
    vect get_col(int col);
    void set_col(int col, vect v);
    vect dot_vector(vect x);
    mat dot_matrix(mat m);
    void delete_row(int nrow);
    void delete_col(int ncol);
    vect back_sub(vect y);
};

double inner_product(vector<double> v1, vector<double> v2);
vector<double> matrix_vector_product(vector<vector<double>> &A, vector<double> &x);

#endif //PROJECT1_MAT_H
