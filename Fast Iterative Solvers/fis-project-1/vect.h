//
// Created by Kenny Yung on 6/8/21.
//

#ifndef PROJECT1_VECT_H
#define PROJECT1_VECT_H
#include <vector>
#include <iostream>
#include <cmath>
using namespace std;

class vect : public vector<double> {
public:
    vect();
    vect(int n, double val);
    vect(vector<double> v);
//    vector<double> m_vect;
//    int m_n;
    vect operator+(vect v);
    vect operator-(vect v);
    vect operator*(double x);
    double operator*(vect v);
    vect operator/(double x);
    vect& operator=(vector<double> &v);
//    double& operator[](int i);
    double norm2();
};

#endif //PROJECT1_VECT_H
