#include "p1solvers.h"

// restart GMRES
vect rGMRES(matrix_csr &A, vect x0, vect b, int m, double tol) {
    ofstream file ("../rgmres.csv");             // output results
    file << "Iteration #,Relative Residual\n";
    vect r0 = b - A.dot_vector(x0);
    double r0norm = r0.norm2();
    vect rk, xm;
    double r;
    int n_iter = 0;
    xm = x0;
    int k = 1;
    do {
        xm = GMRES(A,xm,b,m,k,file,r0norm);
        cout << "xm: ";
        print_vect(xm);
        rk = b - A.dot_vector(xm);
//        rk = A.precond_jacobi(rk);        // precond jacobi
//        rk = A.precond_gs(rk);          // precond gs
        r = rk.norm2()/r0norm;
        n_iter++;
        cout << "iter " << n_iter << " |r|: " << r << endl;
//        file << n_iter << "," << r << endl;     // record residual result to csv
    } while (r > tol);
    file.close();
    return xm;
}
vect GMRES(matrix_csr &A, vect x0, vect b, int m, int &iter, ofstream &restart, double r0n) {
    // GramSchmidt
    int n = A.m_ncol;
    if (m > n)
        cout << "Warning: m exceeds matrix dimension" << endl;
//    ofstream file ("../gmres_gs.csv");
//    file << "iteration k,v1*vk,relative residual\n";
    mat H (m+1, m,0);
    mat v (n, m+1, 0);
    vect Ax = A.dot_vector(x0);
    vect r0 = b - Ax;
//    r0 = A.precond_jacobi(r0);        // Jacobi preconditioner
//    r0 = A.precond_gs(r0);            // Gauss Seidel preconditioner
    double r0norm = r0.norm2();
    v.set_col(0, r0 / r0norm);      // set v1
    vect g (m+1, 0);
    g[0] = r0norm;
    vect w,vj,vj1;
    vect c (m,0);
    vect s (m,0);
    for (int j = 0; j < m; j++) {
        // GetKrylov
        vj = v.get_col(j);
        w = A.dot_vector(vj);
//        w = A.precond_jacobi(w);        // Jacobi precond
//        w = A.precond_gs(w);            // gauss siedel precond
        for (int i = 0; i <= j; i++) {
            vect vi = v.get_col(i);
            H(i,j) = w * vi;
            vect temp = w - vi * H(i,j);
            w = temp;
        }
        H(j+1,j) = w.norm2();
        vj1 = w / H(j+1,j);
        v.set_col(j+1, vj1);
        // GivensRotation
        double temp1,temp2,mag;
        for (int k = 1; k <= j; k++) {
            temp1 = c[k-1]*H(k-1,j) + s[k-1]*H(k,j);
            temp2 = -s[k-1]*H(k-1,j) + c[k-1]*H(k,j);
            H(k-1,j) = temp1;
            H(k,j) = temp2;
        }
        mag = hypot(H(j,j),H(j+1,j));
        c[j] = H(j,j)/mag;
        s[j] = H(j+1,j)/mag;
        temp1 = c[j]*H(j,j) + s[j]*H(j+1,j);
        H(j,j) = temp1;
        g[j+1] = -s[j]*g[j];
        g[j] = c[j]*g[j];
//        file << j << "," << vj*vj1 << "," << abs(g[j+1])/r0norm << endl;
        restart << iter << "," << abs(g[j+1])/r0n << endl;        // for restart gmres
        iter++;
    }
//    file.close();
    mat Rm = H;
    Rm.delete_row(m);
    g.pop_back();
    v.delete_col(m);
    vect Rmg = Rm.back_sub(g);
    vect xm = x0 + v.dot_vector(Rmg);
    return xm;
}
vect Cong_Grad(matrix_csr &A, vect &x0, vect b, double tol) {
    ofstream file;
    file.open("../cg.csv");
    file << "Iteration #,Residual (2-norm),Residual (A-norm)\n";
    vect r0 = b - A.dot_vector(x0);
    double r0norm = r0.norm2();
    vect xs (A.m_ncol,1);
    vect pm = r0;
    vect rm = r0;
    vect xm = x0;
    vect xm1, rm1, pm1, err;
    int n_iter = 0;
    double rk, anorm;
    do {
        vect Apm = A.dot_vector(pm);
        double alpha = rm * rm / (Apm * pm);
        xm1 = xm + pm * alpha;
        rm1 = rm - Apm * alpha;
        double beta = rm1*rm1 / (rm*rm);
        pm1 = rm1 + pm * beta;
//        conv = (xm1-xm).norm2();
        rk = rm1.norm2();
        xm = xm1;
        rm = rm1;
        pm = pm1;
        n_iter++;
        err = xs - xm;
        anorm = sqrt(A.dot_vector(err)*err);
        if (n_iter % 200 == 0) {
            file << n_iter << "," << rk << "," << anorm << endl;
            cout << "iter" << n_iter << " |r|2 = " << rk << " |e|A = " << anorm <<endl;
        }
    } while (rk/r0norm > tol);
    cout << "# of iterations: " << n_iter << endl;
    return xm;
}
matrix_csr read_mtx(string filename) {
    vector<vector<double>> COO;
    bool sym = false;
    ifstream input(filename);
    string line;
    int m,n,entries;
    if (input.is_open()) {
        double row, col, val;
        getline(input, line);
        istringstream iss(line);
        string word;
        while (iss >> word) {
            if (word == "symmetric")
                sym = true;
        }
        while (line.substr(0, 1) == "%") {
            getline(input, line);
        }
        iss = istringstream (line);
        iss >> m;
        iss >> n;
        iss >> entries;
        while (getline(input, line)) {
            istringstream iss(line);
            iss >> row;
            iss >> col;
            iss >> val;
            COO.insert(COO.end(), {row, col, val});
        }
        input.close();
    } else {
        cout << "unable to read file";
    }
    sort(COO.begin(),COO.end());
    vector<double> val;
    vector<int> col;
    vector<int> row;
    int r = 0;
    int r_curr;
    int N = 0;
    row.insert(row.begin(),0);
    for (auto & l : COO) {
        val.insert(val.end(),l[2]);
        col.insert(col.end(),l[1]-1);
        r_curr = l[0]; // current row index
        if (r_curr != r+1) {
            row.insert(row.end(),N);
            r++;
        }
        N++;
    }
    row.insert(row.end(),N);
    return matrix_csr(val,col,row,sym);
}
void print_mat(mat m) {
    for (auto & i : m.m_matrix) {
        for (auto & j : i) {
            cout << j << ", ";
        }
        cout << endl;
    }
}
void print_vect(vect v) {
    for (auto & i : v)
        cout << i << " ";
    cout << endl;
}

void print_COO(vector<vector<double>> matrix) {
    for (auto & i : matrix) {
        for (auto & j : i) {
            cout << j << " ";
        }
        cout << endl;
    }
}
void print_CSR(matrix_csr m) {
    cout << "values: ";
    for (auto & v : m.m_val)
        cout << v << " ";
    cout << endl;
    cout << "# vals " << m.m_val.size() << endl;
    cout << "columns: ";
    for (auto & c : m.m_col)
        cout << c << " ";
    cout << endl;
    cout << "# cols " << m.m_col.size() << endl;
    cout << "row offset: ";
    for (auto & r : m.m_row)
        cout << r << " ";
    cout << endl;
    cout << "# rows " << m.m_row.size() << endl;
}
vector<vector<double>> read_COO(string filename) {
    vector<vector<double>> output;
    bool sym = false;
    ifstream input(filename);
    string line;
    if (input.is_open()) {
        double row, col, val;
        getline(input, line);
        istringstream iss(line);
        string word;
        while (iss >> word) {
            if (word == "symmetric")
                sym = true;
        }
        while (getline(input, line)) {
            if (line.substr(0, 1) != "%") {
                istringstream iss(line);
                iss >> row;
                iss >> col;
                iss >> val;
                output.insert(output.end(), {row, col, val});
//                if (sym && row != col)
//                    output.insert(output.end(),{col,row,val});
            }
        }
        input.close();
    } else {
        cout << "unable to read file";
    }
    return output;
}

vector<vector<double>> COO2CSR(vector<vector<double>> input) {
    vector<double> val;
    vector<double> col;
    vector<double> row;
    sort(input.begin(),input.end());
    int r = 0;
    int r_curr;
    int N = 0;
    row.insert(row.begin(),0);
    for (auto & line : input) {
        val.insert(val.end(),line[2]);
        col.insert(col.end(),line[1]-1);
        r_curr = line[0]; // current row index
        if (r_curr != r+1) {
            row.insert(row.end(),N);
            r++;
        }
        N++;
    }
    row.insert(row.end(),N+1);
    vector<vector<double>> output = {val,col,row};
    return output;
}