//
// Created by Kenny Yung on 6/6/21.
//

#include "mat.h"

mat::mat(int nrow, int ncol, double val) {
    m_matrix = vector<vector<double>>(nrow,vector<double>(ncol, val));
    m_nrow = nrow;
    m_ncol = ncol;
}
mat::mat(vector<vector<double>> v) {
    m_matrix = v;
    m_nrow = v.size();
    m_ncol = v[0].size();
}
double& mat::operator()(int i, int j){
    return m_matrix[i][j];
}
//double mat::operator()(int i, int j) const{
//    return m_matrix[i][j];
//}
vect mat::get_col(int col) {
    vect out;
    for (auto & i : m_matrix)
        out.push_back(i[col]);
    return out;
}
void mat::set_col(int col, vect v) {
    if (col > m_ncol)
        throw invalid_argument("Column specified exceeds matrix column");
    for (int i = 0; i < m_nrow; i++)
        m_matrix[i][col] = v[i];
}
void mat::delete_row(int irow) {
    if (irow >= m_nrow)
        throw invalid_argument("row index exceeds size");
    m_matrix.erase(m_matrix.begin() + irow);
    m_nrow--;
}
void mat::delete_col(int icol) {
    for (auto & i : m_matrix)
        i.erase(i.begin() + icol);
    m_ncol--;
}
vect mat::back_sub(vect y) {
    int n = m_ncol;
    vect x (n, 0);
    for (int i = n-1; i >= 0; i--) {
        double xt = y[i];
        for (int j = i+1; j < n; j++) {
            xt = xt - m_matrix[i][j]*x[j];
        }
        x[i] = xt / m_matrix[i][i];
    }
    return x;
}
mat mat::transpose() {
    mat result(m_ncol,m_nrow,0);
    for (int i = 0; i < m_nrow; i++){
        for (int j = 0; j < m_ncol; j++){
            result(j,i) = m_matrix[i][j];
        }
    }
    return result;
}

vect mat::dot_vector(vect x) {
    if (m_ncol != x.size()) {
        throw invalid_argument("vector size not match");
    }
    vect result(m_nrow,0);
    for (int i = 0; i < m_nrow; i++) {
        for (int j = 0; j < m_ncol; j++) {
            result[i] += m_matrix[i][j] * x[j];
        }
    }
    return result;
}
mat mat::dot_matrix(mat m) {
    if (m_ncol != m.m_nrow) {
        throw invalid_argument("matrix size not match");
    }
    mat result(m_nrow,m.m_nrow,0);
    for (int i = 0; i < m_nrow; i++) {
        for (int j = 0; j < m.m_ncol; j++) {
            double sum = 0;
            for (int k = 0; k < m_ncol; k++) {
                sum += m_matrix[i][k] * m(k,j);
            }
            result(i,j) = sum;
        }
    }
    return result;
}
//vector<vector<double>> transpose(vector<vector<double>> &a) {
//    vector<vector<double>> out (a[0].size(),vector<double>(a.size(),0));
//    for (int i = 0; i < a.size(); i++)
//        for (int j = 0; i < a[0].size(); i++)
//            out[j][i] = a[i][j];
//    return out;
//}
//vector<double> matrix_vector_product(vector<vector<double>> &A, vector<double> &x) {
//    int n = x.size();
//    vector<double> result(n,0);
//    for (int i = 0; i < A.size(); i++) {
//        for (int j = 0; j < n; j++) {
//            result[i] += A[i][j] * x[j];
//        }
//    }
//    return result;
//}
//void array_add(vector<double> &a, double b) {
//    for (auto & i : a)
//        i += b;
//}
//void array_multiply(vector<double> &a, double b) {
//    for (auto & i : a)
//        i = i * b;
//}
//double inner_product(vector<double> v1, vector<double> v2){
//    double result = 0;
//    if (v1.size() != v2.size()){
//        cout << "vector length not match" << endl;
//        return -1;
//    }
//    for (int i = 0; i < v1.size(); i++)
//        result += v1[i] * v2[i];
//    return result;

