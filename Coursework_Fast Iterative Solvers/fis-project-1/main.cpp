#include "p1solvers.h"
#include "mat.h"
#include <chrono>

void test();
void full_gmres();
void rgmres();
void conj_grad();

using namespace std;
using namespace chrono;

int main() {

//    test();
//    full_gmres();
    rgmres();
//    conj_grad();

    return 0;
}

void full_gmres() {
    string input;
    matrix_csr A;
    int m,n;
    vect x0,b,sol,xm;
    double tol = 1e-8;

    input = "../orsirr_1.mtx";
    A = read_mtx(input);
    n = A.m_ncol;
    cout << "n = " << n << endl;
    sol  = vect(n,1);     // prescribe solution as {1,1,...1}
    b = A.dot_vector(sol);
    x0 = vect(n,0);     // initial guess {0,0,...0}

    m = 165;
    double conv;
    double r0norm = (b - A.dot_vector(x0)).norm2();
    do {
        cout << "Full GMRES with m = " << m << endl;
        int i = 0;
        ofstream file;
        auto t1 = high_resolution_clock::now();
        xm = GMRES(A, x0, b, m, i, file,1);
        auto t2 = high_resolution_clock::now();
        auto ms_int = duration_cast<milliseconds>(t2 - t1);
        conv = (b - A.dot_vector(xm)).norm2() / r0norm;
        cout << "xm: ";
        print_vect(xm);
        cout << "r = " << conv << endl;
        cout << "Time duration: " << ms_int.count() << "ms" << endl;
        m = m + 1;
    } while (conv > tol);
}

void rgmres() {
    string input;
    matrix_csr A;
    int m,n;
    vect x0,b,sol,xm;
    double tol = 1e-8;
    input = "../orsirr_1.mtx";
    A = read_mtx(input);
    n = A.m_ncol;
    cout << "n = " << n << endl;
    sol  = vect(n,1);     // prescribe solution as {1,1,...1}
    b = A.dot_vector(sol);
    x0 = vect(n,0);     // initial guess {0,0,...0}
    // Restart GMRES
    m = 10;
    cout << "Restarted GMRES with m = " << m << endl;
    auto t1 = high_resolution_clock::now();
    xm = rGMRES(A,x0,b,m,tol);
    auto t2 = high_resolution_clock::now();
    auto ms_int = duration_cast<milliseconds>(t2 - t1);
    cout << "Time duration: " << ms_int.count() << "ms" << endl;
}

void conj_grad() {
    // conjugate gradient
    string input;
    matrix_csr A;
    int m,n;
    vect x0,b,sol,xm;
    double tol = 1e-8;

    cout << "Conjugate Gradient Method:" << endl;
    input = "../s3rmt3m3.mtx";
    A = read_mtx(input);
    n = A.m_ncol;
    sol = vect(n,1);
    b = A.dot_vector(sol);
    x0 = vect(n,0);
    auto t1 = high_resolution_clock::now();
    xm = Cong_Grad(A,x0,b,tol);
    auto t2 = high_resolution_clock::now();
    auto ms_int_cg = duration_cast<milliseconds>(t2 - t1);
    cout << "xm: ";
    print_vect(xm);
    cout << "Time duration: " << ms_int_cg.count() << "ms" << endl;
}

void test() {
    //test 3x3 matrix
//    input = "../test.txt";
//    x0 = vect({1,1,1});
//    b  = vect({5,4,3});
//    A = read_mtx(input);
//matrix = read_COO(input);
    //matrix = COO2CSR(test);
    //print_COO(matrix);
    //    print_CSR(A);
//    vect x ({1,2,3});
//    vect y = A.dot_vector(x);
//    cout << "result:" << endl;
//    for (auto & i : y){
//        cout << i << endl;
//    }
//    b  = vect({1,3,5});
//    matrix_csr matrix = read_mtx("../test2.txt");
//    print_vect(matrix.forward_sub(vect({1,3,5})));
}
