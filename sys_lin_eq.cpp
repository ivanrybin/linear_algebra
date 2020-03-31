/* 
   Ivan Rybin 2019. 
   Programm solves system of linear equations. 
*/
#include <iostream>
#include <iomanip>
#include <string>
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
        if (abs(A[i][A.size() - 1]) > 0.000000001) // troubles with double's precision
            for (int k = i - 1; k >= 0; --k) {
                double p = -A[k][A.size() - 1]/A[i][A.size() - 1];
                for (int j = 0; j < A[0].size(); ++j)
                    A[k][j] += A[i][j] * p;
            }
    }
}

int Rank(const vector<vector<double>>& A) {
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

void SeparateAugmentedMatrix (const vector<vector<double>>& Ab, 
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
        cout << "\n";
    }
}

int main() {
    int m, n; // m - number of rows, n - number of columns
    cout << "Input equations number: ";
    cin >> m;
    cout << "Input variables number: ";
    cin >> n;
    cout << "\n" << "Example of input 3x4 augmented matrix (A|b):" << "\n" << "a1 a2 a3 b1" <<
                                 "\n" << "a4 a5 a6 b2" << "\n" << "a7 a8 a9 b3" << "\n" << "\n";
    cout << "Your input:" << "\n";
    vector<vector<double>> Ab = InputAugmentedMatrix(m, n + 1); // augmented matrix
    cout << "\n";

    vector<vector<double>> A(m, vector<double> (n, 0)); // matrix A
    vector<double> b(m, 0); // matrix B
    SeparateAugmentedMatrix(Ab, A, b);

    double det_A = 1; 
    Gauss(Ab, det_A);
    ZeroRowKiller(Ab);
    DetAfterGauss(Ab, det_A);

    vector<vector<double>> A2(m, vector<double> (n, 0));
    vector<double> b2(m, 0);
    SeparateAugmentedMatrix(Ab, A2, b2);
    int rank_Ab = Rank(Ab);
    int rank_A = Rank(A2);

    if (rank_A < rank_Ab) { // no solution
        cout << "The system has NO solution <= rank(A) = " << rank_A <<
                    " < " << "rank(A|b) = " << rank_Ab << " (Kronecker Capelli theorem)" << "\n" << "\n";
    } else if (rank_Ab == rank_A && rank_Ab < n) { // inf many solutions
        cout << "The system has INFINITELY many solutions <= rank(A|b) = rank(A) = " << 
                                                        rank_A << " < n = " << n << "\n" << "\n";
    } else if (rank_A == rank_Ab && rank_A == n) { // one solution (Cramer's rule)
        cout << "The system has a single unique solution" <<
                         " <= rank(A|b) = rank(A) = n = " << n << "\n" << "\n";
        vector<double> cramers_dets(n, 1);
        for (int j = 0; j < n; ++j) {
            vector<vector<double>> bX = A; // Cramer
            for (int i = 0; i < bX.size(); ++i)
                bX[i][j] = b[i];
            Gauss(bX, cramers_dets[j]);
            ZeroRowKiller(bX);
            DetAfterGauss(bX, cramers_dets[j]);
        }
        cout << "Solution:" << "\n";
        for (int i = 0; i < cramers_dets.size(); ++i)
            cout << "X" << i + 1 << " = " << setw(12) << right << setprecision(6) << fixed <<(cramers_dets[i] / det_A) << "\n";
        cout << "\n" << "Cramer's rule: Xi = det_i / det(A)" << "\n";
        cout << "\n" << "Det(A) = " << setw(12) << right << setprecision(3) << det_A << "\n";
        cout << "\n" << "Cramer's determinants:" << "\n";
        for (int i = 0; i < cramers_dets.size(); ++i)
            cout << "det #" << i + 1 << " = " << setw(12) << right << setprecision(3) << fixed << cramers_dets[i] << "\n";
        cout << "\n";
    }
    cout << "Augmented matrix (A|b) after Gaussian elimination:" << "\n";
    ShowMatrix(Ab);
    return 0;
}
