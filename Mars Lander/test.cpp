
#include<bits/stdc++.h>

using namespace std;

class Matrix {
public:
    int row, column;
    float **data;
    Matrix() : row(0), column(0) {}
    Matrix(int row, int column) : row(row), column(column) {
        data = new float*[row];
        for (int i = 0; i < row; i++) {
            data[i] = new float[column];
        }
    }
    void set(int i, int j, float value) {
        data[i][j] = value;
    }
    float get(int i, int j) {
        return data[i][j];
    }
    void print() {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                cout << data[i][j] << " ";
            }
            cout << endl;
        }
    }
    Matrix operator*(Matrix m) {
        Matrix result(row, m.column);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < m.column; j++) {
                float sum = 0;
                for (int k = 0; k < column; k++) {
                    sum += data[i][k] * m.data[k][j];
                }
                result.set(i, j, sum);
            }
        }
        return result;
    }
    Matrix operator+(Matrix m) {
        Matrix result(row, column);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                result.set(i, j, data[i][j] + m.data[i][j]);
            }
        }
        return result;
    }
    Matrix operator-(Matrix m) {
        Matrix result(row, column);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                result.set(i, j, data[i][j] - m.data[i][j]);
            }
        }
        return result;
    }
    Matrix transpose() {
        Matrix result(column, row);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                result.set(j, i, data[i][j]);
            }
        }
        return result;
    }

};

int main() {
    Matrix A;
    A = Matrix(3, 3);
    A.print();
    cerr << 3 << endl;

}