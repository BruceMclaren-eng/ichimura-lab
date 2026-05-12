#include <iostream>
#include <vector>

struct CRSMatrix {
    int n_rows;
    int n_cols;
    int nnz;
    std::vector<double> values;
    std::vector<int> col_idx;
    std::vector<int> row_ptr;
};

void init_matrix(CRSMatrix& A, int n_rows, int n_cols, int nnz){
    A.n_rows = n_rows;
    A.n_cols = n_cols;
    A.nnz = nnz;
    A.values.resize(nnz, 0.0);
    A.col_idx.resize(nnz, 0);
    A.row_ptr.resize(n_rows + 1,0);
}

void print_matrix(const CRSMatrix& A) {
    std::cout << "n_rows: " << A.n_rows << "\n";
    std::cout << "n_cols: " << A.n_cols << "\n";
    std::cout << "nnz: " << A.nnz << "\n";
}

int main() {
    CRSMatrix A;
    init_matrix(A,4,4,8);
    print_matrix(A);

    return 0;
}