//
// Created by Kenny Yung on 6/8/21.
//

#include "matrix_csr.h"
matrix_csr::matrix_csr() {}
matrix_csr::matrix_csr(vector<double> &val, vector<int> &col, vector<int> &row, bool sym){
    m_val = val;
    m_col = col;
    m_row = row;
    m_sym = sym;
    m_nrow = row.size()-1;
    m_ncol = col.at(distance(col.begin(),max_element(col.begin(),col.end())))+1; // max(col)
}
vect matrix_csr::dot_vector(vect x) {
    int n = x.size();
    if (n != m_ncol)
        throw invalid_argument("Vector size not match");
    vect result(n,0);
    for (int r = 0; r < m_row.size()-1; r++) {
        for (int j = m_row[r]; j < m_row[r+1]; j++) {
            int c = m_col[j];
            result[r] += m_val[j] * x[c];
            if (m_sym && r != c)
                result[c] += m_val[j] * x[r];
        }
    }
    return result;
}
vect matrix_csr::forward_sub(vect y) {
    int n = y.size();
    if (n != m_ncol)
        throw invalid_argument("vector size and column not match");
    vect x (n, 0);
    for (int i = 0; i < n; i++) {
        double alpha = y[i];
//        for (int j = i+1; j < n; j++) {
        double Lii;
        for (int j = m_row[i]; j < m_row[i+1]; j++) {
            int c = m_col[j];
            double Lij = m_val[j];
            if (c < i)
                alpha = alpha - Lij * x[c];
            if (c == i)
                Lii = Lij;
        }
        x[i] = alpha / Lii;
    }
    return x;
}
vect matrix_csr::precond_jacobi(vect v) {
    // return M (diagonal of A) ^-1
    vector<double> val;
    vector<int> row,col;
    for (int r = 0; r < m_row.size()-1; r++) {
        for (int j = m_row[r]; j < m_row[r+1]; j++) {
            int c = m_col[j];
            if (c == r) {
                double temp = m_val[j];
                val.push_back(1/temp);
                row.push_back(r);
                col.push_back(c);
            }
        }
    }
    row.push_back(col.size());
    matrix_csr M_inv (val,col,row,false);
    return M_inv.dot_vector(v);
}
vect matrix_csr::precond_gs(vect v) {
    // return M (D - L) lower triangular
    vector<double> val;
    vector<int> row,col;
    row.push_back(0);
    int count = 0;
    for (int r = 0; r < m_row.size()-1; r++) {
        for (int j = m_row[r]; j < m_row[r+1]; j++) {
            int c = m_col[j];
            if (c <= r) {
                double temp = m_val[j];
                val.push_back(-temp);
                col.push_back(c);
                count++;
            }
        }
        row.push_back(count);
    }
//    row.push_back(val.size());
    matrix_csr M (val,col,row,false);
    return M.forward_sub(v);
}