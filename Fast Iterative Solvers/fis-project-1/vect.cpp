//
// Created by Kenny Yung on 6/8/21.
//

#include "vect.h"

vect::vect(int n, double val){
    for (int i = 0; i < n; i++)
        this->push_back(val);
//    m_vect = vector<double>(n,val);
}
vect::vect()= default;
vect::vect(vector<double> v){
    this->clear();
    for (auto & i : v)
        this->push_back(i);
}
vect& vect::operator=(vector<double> &v) {
    this->clear();
    for (auto & i : v){
        this->push_back(i);
    }
//    m_n = this->size();
    return *this;
}
vect vect::operator+(vect v){
    if (v.size() != this->size())
        throw invalid_argument("vector size not match");
    vect out;
    for (int i = 0; i < this->size(); i++)
        out.push_back(v[i] + this->at(i));
    return out;
}
vect vect::operator-(vect v){
    if (v.size() != this->size())
        throw invalid_argument("vector size not match");
    vect out;
    for (int i = 0; i < this->size(); i++)
        out.push_back(this->at(i) - v[i]);
    return out;
}
vect vect::operator*(double x){
    vect out;
    for (int i = 0; i < this->size(); i++)
        out.push_back(this->at(i) * x);
    return out;
}
double vect::operator*(vect v){
    if (v.size() != this->size())
        throw invalid_argument("vector size not match");
    double out = 0;
    for (int i = 0; i < v.size(); i++)
        out += this->at(i) * v[i];
    return out;
}
vect vect::operator/(double x){
    vect out(this->size(),0);
    for (int i = 0; i < this->size(); i++)
        out[i] = (this->at(i) / x);
    return out;
}
double vect::norm2() {
    double sum = 0;
    for (auto & i : *this)
        sum += pow(i, 2);
    return sqrt(sum);
}
//double& vect::operator[](int i) {
//    return m_vect[i];
//}

//int vect::length() {
//    return m_vect.size();
//}
//void vect::push_back(double x) {
//    m_vect.push_back(x);
//}
