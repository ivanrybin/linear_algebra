/* Ivan Rybin 2019.
   Programm finds the best solution for unsolvable system of linear equations, where A (n > m). */
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

vector<vector<double>> InputAugmentedMatrix(int m, int n) {
    // m - number of rows, n - number of columns
    vector<vector<double>> A(m, vector<double>(n, 0));
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j) {
            int element;
            cin >> element;
            A[i][j] = element;
        }
    return A;
}

vector<vector<double>> TransposeMatrix(const vector<vector<double>>& A) {
    vector<vector<double>> At(A[0].size(), vector<double> (A.size()));
    for (int i = 0; i < A.size(); ++i)
        for (int j = 0; j < A[0].size(); ++j)
            At[j][i] = A[i][j];
    return At;
}

bool ZeroRowColTest(const vector<vector<double>>& A) {
    // zero row (column) test 
    for (int i = 0; i < A.size(); ++i) {
        int r_count = 0;
        int c_count = 0;
        for (int j = 0; j < A.size(); ++j) {
            if (A[i][j] == 0)
                r_count++;
            if (A[j][i] == 0)
                c_count++;
        }
        if (r_count == A.size() || c_count == A.size())
            return true; // det(A) = 0
    }
    return false;
}

void Gauss(vector<vector<double>>& A, double& det) {
    for (int column = 0; column < A.size() - 1; ++column) {
        // search of max element (number of row) in current column
        double max_element = A[column][column];
        int max_row = column;
        for (int i = column; i < A.size(); ++i)
            if (abs(A[i][column]) > abs(max_element)) {
                max_element = A[i][column];
                max_row = i;
            }
        // permutation max_row and current row
        if (max_element != A[column][column]) {
            vector<double> tmp = A[max_row];
            A[max_row] = A[column];
            A[column] = tmp;
            det *= -1;
        }
        // Gauss down in current column
        if (A[column][column] != 0)
            for (int row = column + 1; row < A.size(); ++row) {
                double p = -A[row][column]/A[column][column];
                for (int j = column; j < A[0].size(); ++j)
                    A[row][j] += A[column][j] * p; 
            }
    }
    // Gauss up (deletes linear dependet rows)
    for (int i = A.size() - 1; i > 0; --i) {
        if (abs(A[i][A.size() - 1]) > 0.000000001) // troubles with double precision
            for (int k = i - 1; k >= 0; --k) {
                double p = -A[k][A.size() - 1]/A[i][A.size() - 1];
                for (int j = 0; j < A[0].size(); ++j)
                    A[k][j] += A[i][j] * p;
            }
    }
}

int Rank(vector<vector<double>>& A) {
    int rank = A.size();
    for (int i = A.size() - 1 ; i >= 0; --i) {
        int zero_count = 0;
        for (int j = 0; j < A[0].size(); ++j)
            if(abs(A[i][j]) <= 0.000000001)
                zero_count++;
        if (zero_count == A[0].size())
            rank--;
    }
    return rank;
}

void DetAfterGauss(const vector<vector<double>>& A, double& det) {
    for (int i = 0; i < A.size(); ++i)
        det *= A[i][i];
}

void ZeroRowKiller(vector<vector<double>>& A) {
    vector<vector<double>> B;
    for (int i = 0; i < A.size(); ++i) {
        int zero_count = 0;
        for (int j = 0; j < A[0].size(); ++j)
            if(abs(A[i][j]) <= 0.000000001)
                zero_count++;
        if (zero_count != A[0].size())
            B.push_back(A[i]);
    }
    A = B;
}

vector<double> MatrixVectorMult(const vector<vector<double>>& A,
                                         const vector<double>& b) {
    vector<double> D(A.size());
    for (int i = 0; i < A.size(); ++i) {
        double sum = 0;
        for (int j = 0; j < b.size(); ++j)
            sum += A[i][j] * b[j];
        D[i] = sum;
    }
    return D;
}

vector<vector<double>> MatrixMatrixMult(const vector<vector<double>>& A, 
                                            const vector<vector<double>>& B) {
    vector<vector<double>> C(A.size(), vector<double> (B[0].size()));
    for (int k = 0; k < C[0].size(); ++k)
        for (int i = 0; i < A.size(); ++i) {
            double sum = 0;
            for (int j = 0; j < B.size(); ++j)
                sum += A[i][j] * B[j][k];
            C[i][k] = sum;
        }
    return C;
}

void SeparateAugmentedMatrix(const vector<vector<double>>& Ab, 
                            vector<vector<double>>& A, vector<double>& b) {
    for (int i = 0; i < Ab.size(); ++i)
        for (int j = 0; j < Ab[0].size() - 1; ++j)
            A[i][j] = Ab[i][j];
    for (int i = 0; i < Ab.size(); ++i)
        b[i] = Ab[i][Ab[0].size() - 1];
}

void ShowMatrix(const vector<vector<double>>& A) {
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[0].size(); ++j) 
            cout << fixed << setprecision(3) << right << setw(8) << A[i][j];
        cout << endl;
    }
}

int main() {
    int m, n; // m - number of rows, n - number of columns
    cout << "Input equations number: ";
    cin >> m;
    cout << "Input variables number (< equations number): ";
    cin >> n;
    cout << "\n" << "Example of input 4x2 augmented matrix (A|b):" << "\n" << "a1 a2 b1" <<
                                 "\n" << "a3 a4 b2" << "\n" << "a5 a6 b3" << "\n" << "a7 a8 b4" << "\n" << "\n";
    cout << "Your input:" << "\n";
    vector<vector<double>> Ab = InputAugmentedMatrix(m, n + 1);
    vector<vector<double>> A(m, vector<double> (n, 0));
    vector<double> b(m, 0);
    SeparateAugmentedMatrix(Ab, A, b);

    vector<vector<double>> At = TransposeMatrix(A);
    vector<vector<double>> C = MatrixMatrixMult(At, A);
    vector<double> D = MatrixVectorMult(At, b);

    // CX = D
    vector<vector<double>> CD (C.size(), vector<double> (C[0].size() + 1));
    for (int i = 0; i < CD.size(); ++i)
        for (int j = 0; j < CD[0].size() - 1; ++j)
            CD[i][j] = C[i][j];
    
    for (int i = 0; i < CD.size(); ++i)
        CD[i][CD[0].size() - 1] = D[i];

    // X_hat
    vector<vector<double>> C2 = C;
    double det_C = 1;
    Gauss(C, det_C);
    DetAfterGauss(C, det_C);
    vector<double> cramers_dets(C.size(), 1);
    vector<double> x_hat;

    for (int j = 0; j < C.size(); ++j) {
        vector<vector<double>> bX = C2; // Cramer
        for (int i = 0; i < bX.size(); ++i)
            bX[i][j] = D[i];
        Gauss(bX, cramers_dets[j]);
        DetAfterGauss(bX, cramers_dets[j]);
    }

    for (const auto& det: cramers_dets)
        x_hat.push_back(det/det_C);
    
    cout << "\n" << "The best solution for Ax = b with Ordinary least squares (OLS):" << "\n";
    for (int i = 0; i < x_hat.size(); ++i)
        cout << "X" << i + 1 << " = " << setw(12) << right << setprecision(6) << fixed << x_hat[i] << "\n";

    return 0;
}
