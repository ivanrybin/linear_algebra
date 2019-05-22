/* 
Ivan Rybin 2019.
Programm calculates determinant of quadratic matrix NxN with Gaussian elimination (rref).
*/
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

vector<vector<double>> InputMatrix(int n) {
    vector<vector<double>> A(n, vector<double>(n, 0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            int element;
            cin >> element;
            A[i][j] = element;
        }
    return A;
}

bool ZeroRowColTest(const vector<vector<double>>& A) {
    // zero row (column) test 
    for (int i = 0; i < A.size(); ++i) {
        int count = 0;
        for (int j = 0; j < A.size(); ++j)
            if (A[i][j] == 0)
                count++;
        if (count == A.size())
            return true; // det(A) = 0
    }
    for (int j = 0; j < A.size(); ++j) {
        int count = 0;
        for (int i = 0; i < A.size(); ++i) 
            if (A[i][j] == 0)
                count++;
        if (count == A.size())
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

void ShowMatrix(const vector<vector<double>>& A) {
    int max = 0;
    int width = 8;
    for (int i = 0; i < A.size(); ++i)
        for (int j = 0; j < A.size(); ++j) 
            if (abs(A[i][j] > max))
                max = abs(A[i][j]);
    while(max != 0) { max /= 10; width++;}
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A.size(); ++j) 
            cout << setprecision(3) << right << setw(width) << A[i][j];
        cout << endl;
    }
}

int main() {
    cout << "Input matrix size n: ";
    int n;
    cin >> n;
    cout << "Input elements of " << n << "x" << n << " matrix:" << endl;
    vector<vector<double>> A = InputMatrix(n);
    cout << endl;
    if (ZeroRowColTest(A)) {
        cout << "Det(A) = 0. Matrix has zero row (column).";
    } else {
        double det = 1;
        Gauss(A, det);
        for (int i = 0; i < n; ++i)
            det *= A[i][i];
        cout << setprecision(3) << fixed;
        cout << "Det(A) = " << round(det) << endl;
        cout << endl;
        cout << "Matrix after Gauss elimination:" << endl;
        ShowMatrix(A);
    }
    return 0;
}
