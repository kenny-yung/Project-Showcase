//
// Created by Kenny Yung on 6/8/21.
//

#ifndef PROJECT1_MATRIX_CSR_H
#define PROJECT1_MATRIX_CSR_H
#include <iostream>
#include <cmath>
#include "vect.h"

class matrix_csr {
public:
    matrix_csr();
    matrix_csr(vector<double> &val, vector<int> &col, vector<int> &row, bool sym);
    vector<double> m_val;
    vector<int> m_col, m_row;
    bool m_sym;
    int m_nrow;
    int m_ncol;
    vect dot_vector(vect x);
    vect forward_sub(vect vect);
    vect precond_jacobi(vect v);
    vect precond_gs(vect v);
};


#endif //PROJECT1_MATRIX_CSR_H
