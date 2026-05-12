#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <numeric>

// ============================================================
// fem_types モジュール相当
// Fortran: type :: pair / real(4) :: x, y
// ============================================================
struct Pair {
    float x, y;
};

// CRS行列構造体
// Fortran: row_ptr, col_idx, values をバラバラに管理していたものをまとめる
struct CRSMatrix {
    int n_rows;
    int nnz;
    std::vector<float> values;
    std::vector<int>   col_idx;
    std::vector<int>   row_ptr;
};

// ============================================================
// グローバル変数（Fortranのprogram変数に相当）
// ============================================================
int nx, ny;
int n_nodes, n_elems, n_edges;
int nnz_global;
std::vector<Pair>  nodes;
std::vector<std::vector<int>> elems; // elems[e][0..2]

// ============================================================
// generate_nodes
// Fortran: subroutine generate_nodes(nx,ny,x_max,x_min,y_max,y_min,nodes)
// ============================================================
void generate_nodes(int nx, int ny,
                    float x_max, float x_min,
                    float y_max, float y_min,
                    std::vector<Pair>& nodes)
{
    float dx = (x_max - x_min) / nx;
    float dy = (y_max - y_min) / ny;

    int node_id = 0;
    for (int j = 0; j <= ny; j++) {
        for (int i = 0; i <= nx; i++) {
            nodes[node_id].x = x_min + i * dx;
            nodes[node_id].y = y_min + j * dy;
            node_id++;
        }
    }
}

// ============================================================
// generate_elems
// Fortran: subroutine generate_elems(nx,ny,elems)
// ============================================================
void generate_elems(int nx, int ny,
                    std::vector<std::vector<int>>& elems)
{
    int elem_id = 0;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int n1 = j * (nx + 1) + i;         // 0始まりに変換
            int n2 = n1 + 1;
            int n3 = n1 + nx + 1;
            int n4 = n3 + 1;

            // triangle 1: (n1, n2, n4)
            elems[elem_id][0] = n1;
            elems[elem_id][1] = n2;
            elems[elem_id][2] = n4;
            elem_id++;

            // triangle 2: (n1, n4, n3)
            elems[elem_id][0] = n1;
            elems[elem_id][1] = n4;
            elems[elem_id][2] = n3;
            elem_id++;
        }
    }
}

// ============================================================
// calc_elem_stiffness
// Fortran: subroutine calc_elem_stiffness(p1,p2,p3,k_mat)
// ============================================================
void calc_elem_stiffness(const Pair& p1, const Pair& p2, const Pair& p3,
                         float k_mat[3][3])
{
    float x[3] = {p1.x, p2.x, p3.x};
    float y[3] = {p1.y, p2.y, p3.y};

    // ヤコビアン
    float J[2][2];
    J[0][0] = -x[0] + x[1];
    J[0][1] = -y[0] + y[1];
    J[1][0] = -x[0] + x[2];
    J[1][1] = -y[0] + y[2];

    float detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
    float val_m = 1.0f / detJ;
    float area  = std::abs(detJ) * 0.5f;

    // ヤコビアン逆行列
    float L[2][2];
    L[0][0] =  val_m * J[1][1];
    L[0][1] = -val_m * J[0][1];
    L[1][0] = -val_m * J[1][0];
    L[1][1] =  val_m * J[0][0];

    // 局所微分
    float v[3][2] = {{-1.0f, -1.0f}, {1.0f, 0.0f}, {0.0f, 1.0f}};

    // 全体座標微分
    float res[3][2];
    for (int i = 0; i < 3; i++) {
        res[i][0] = J[0][0]*v[i][0] + L[1][0]*v[i][1];
        res[i][1] = J[0][1]*v[i][0] + L[1][1]*v[i][1];
    }

    // 要素剛性行列
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            k_mat[i][j] = (res[i][0]*res[j][0] + res[i][1]*res[j][1]) * area;
        }
    }
}

// ============================================================
// calc_global_stiffness
// Fortran: subroutine calc_global_stiffness(...)
// ============================================================
void calc_global_stiffness(const std::vector<std::vector<std::vector<float>>>& K_semilocal,
                           CRSMatrix& A)
{
    const int max_nnz = 10;
    std::vector<int>              nnz_row(n_nodes, 0);
    std::vector<std::vector<int>> col_idx_temp(n_nodes, std::vector<int>(max_nnz, 0));

    // スパースパターンの構築
    for (int e = 0; e < n_elems; e++) {
        for (int i = 0; i < n_edges; i++) {
            int gi = elems[e][i];
            for (int j = 0; j < n_edges; j++) {
                int gj = elems[e][j];
                bool already_exists = false;
                for (int k = 0; k < nnz_row[gi]; k++) {
                    if (col_idx_temp[gi][k] == gj) {
                        already_exists = true;
                        break;
                    }
                }
                if (!already_exists) {
                    if (nnz_row[gi] >= max_nnz) {
                        std::cerr << "Error: Exceeded maximum non-zeros per row\n";
                        std::exit(1);
                    }
                    col_idx_temp[gi][nnz_row[gi]++] = gj;
                }
            }
        }
    }

    // row_ptr の構築
    A.row_ptr.resize(n_nodes + 1);
    A.row_ptr[0] = 0;  // C++は0始まり（Fortranは1始まりだった）
    for (int gi = 0; gi < n_nodes; gi++) {
        A.row_ptr[gi+1] = A.row_ptr[gi] + nnz_row[gi];
    }
    nnz_global = A.row_ptr[n_nodes];
    A.nnz = nnz_global;

    // col_idx の構築
    A.col_idx.resize(nnz_global);
    int crs_pos = 0;
    for (int gi = 0; gi < n_nodes; gi++) {
        for (int ptr = 0; ptr < nnz_row[gi]; ptr++) {
            if (col_idx_temp[gi][ptr] != 0 || ptr == 0) {
                A.col_idx[crs_pos++] = col_idx_temp[gi][ptr];
            }
        }
    }

    // values の構築
    A.values.resize(nnz_global, 0.0f);
    for (int e = 0; e < n_elems; e++) {
        for (int i = 0; i < n_edges; i++) {
            int gi = elems[e][i];
            for (int j = 0; j < n_edges; j++) {
                int gj = elems[e][j];
                for (int k = A.row_ptr[gi]; k < A.row_ptr[gi+1]; k++) {
                    if (A.col_idx[k] == gj) {
                        A.values[k] += K_semilocal[e][i][j];
                        break;
                    }
                }
            }
        }
    }
}

// ============================================================
// dirichlet_row
// Fortran: subroutine dirichlet_row(node, n, values, col_idx, row_ptr, F_vec, bc_val)
// ============================================================
void dirichlet_row(int node, CRSMatrix& A,
                   std::vector<float>& F_vec, float bc_val)
{
    for (int k = A.row_ptr[node]; k < A.row_ptr[node+1]; k++) {
        A.values[k] = (A.col_idx[k] == node) ? 1.0f : 0.0f;
    }
    F_vec[node] = bc_val;
}

// ============================================================
// matvec（SpMV）
// Fortran: subroutine matvec2(n, values, col_idx, row_ptr, x, Ap)
// ============================================================
void matvec(const CRSMatrix& A,
            const std::vector<float>& x,
            std::vector<float>& Ap)
{
    for (int i = 0; i < A.n_rows; i++) {
        Ap[i] = 0.0f;
        for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; k++) {
            Ap[i] += A.values[k] * x[A.col_idx[k]];
        }
    }
}

// ============================================================
// pcg_solver
// Fortran: subroutine pcg_solver(n, values, col_idx, row_ptr, b, x, diag)
// ============================================================
void pcg_solver(const CRSMatrix& A,
                const std::vector<float>& b,
                std::vector<float>& x,
                const std::vector<float>& diag)
{
    int n = A.n_rows;
    std::vector<float> r(n), p(n), Ap(n), z(n);

    // 初期残差
    matvec(A, x, Ap);
    float r0_sq = 0.0f;
    for (int i = 0; i < n; i++) {
        r[i]  = b[i] - Ap[i];
        z[i]  = r[i] / diag[i];
        p[i]  = z[i];
        r0_sq += r[i] * r[i];
    }

    float rz = 0.0f;
    for (int i = 0; i < n; i++) rz += r[i] * z[i];

    double t_spmv = 0, t_pap = 0, t_update = 0, t_pvec = 0;
    int iter = 0;

    for (iter = 0; iter < 100000; iter++) {

        auto a1 = std::chrono::high_resolution_clock::now();
        matvec(A, p, Ap);
        auto a2 = std::chrono::high_resolution_clock::now();
        t_spmv += std::chrono::duration<double>(a2-a1).count();

        auto b1 = std::chrono::high_resolution_clock::now();
        float pAp = 0.0f;
        for (int i = 0; i < n; i++) pAp += p[i] * Ap[i];
        auto b2 = std::chrono::high_resolution_clock::now();
        t_pap += std::chrono::duration<double>(b2-b1).count();

        float alpha = rz / pAp;

        auto c1 = std::chrono::high_resolution_clock::now();
        float rz_new = 0.0f, r_norm_sq = 0.0f;
        for (int i = 0; i < n; i++) {
            x[i]     += alpha * p[i];
            r[i]     -= alpha * Ap[i];
            z[i]      = r[i] / diag[i];
            rz_new   += r[i] * z[i];
            r_norm_sq += r[i] * r[i];
        }
        auto c2 = std::chrono::high_resolution_clock::now();
        t_update += std::chrono::duration<double>(c2-c1).count();

        float beta = rz_new / rz;

        auto d1 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < n; i++) p[i] = z[i] + beta * p[i];
        auto d2 = std::chrono::high_resolution_clock::now();
        t_pvec += std::chrono::duration<double>(d2-d1).count();

        rz = rz_new;
        if (std::sqrt(r_norm_sq) / std::sqrt(r0_sq) < 1.0e-8f) break;
    }

    std::cout << "Number of iterations: " << iter << "\n";
    std::cout << "t_spmv:  " << t_spmv  << " seconds\n";
    std::cout << "t_pap:   " << t_pap   << " seconds\n";
    std::cout << "t_update:" << t_update << " seconds\n";
    std::cout << "t_pvec:  " << t_pvec  << " seconds\n";
}

// ============================================================
// main
// ============================================================
int main() {
    auto t_total_start = std::chrono::high_resolution_clock::now();

    // Step1: メッシュ定義
    float x_min = 0.0f, x_max = 1.0f;
    float y_min = 0.0f, y_max = 2.0f;
    nx = 15000;
    ny = 15000;
    n_nodes = (nx+1) * (ny+1);
    n_elems = nx * ny * 2;
    n_edges = 3;

    // Step2: 配列初期化
    nodes.resize(n_nodes);
    elems.resize(n_elems, std::vector<int>(3));
    std::vector<float> F_vec(n_nodes, 0.0f);
    std::vector<float> u_vec(n_nodes, 0.0f);

    // Step3: メッシュ生成
    auto a1 = std::chrono::high_resolution_clock::now();
    generate_nodes(nx, ny, x_max, x_min, y_max, y_min, nodes);
    generate_elems(nx, ny, elems);
    auto a2 = std::chrono::high_resolution_clock::now();
    std::cout << "Mesh Generated\n";
    std::cout << "Generated Nodes:    " << n_nodes << "\n";
    std::cout << "Generated Elements: " << n_elems << "\n";

    // Step4: 要素剛性行列
    // Fortranのreal(4) K_semilocal(n_elems, n_edges, n_edges) に相当
    std::vector<std::vector<std::vector<float>>>
        K_semilocal(n_elems, std::vector<std::vector<float>>(3, std::vector<float>(3)));

    auto c1 = std::chrono::high_resolution_clock::now();
    for (int e = 0; e < n_elems; e++) {
        float k_mat[3][3];
        calc_elem_stiffness(nodes[elems[e][0]],
                            nodes[elems[e][1]],
                            nodes[elems[e][2]], k_mat);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                K_semilocal[e][i][j] = k_mat[i][j];
    }
    auto c2 = std::chrono::high_resolution_clock::now();

    // Step4続き: 全体剛性行列（CRS変換）
    CRSMatrix A;
    A.n_rows = n_nodes;

    auto d1 = std::chrono::high_resolution_clock::now();
    calc_global_stiffness(K_semilocal, A);
    auto d2 = std::chrono::high_resolution_clock::now();
    std::cout << "CRS Conversion Done: nnz = " << nnz_global << "\n";

    // Step5: 境界条件
    std::vector<bool> is_fixed(n_nodes, false);
    auto e1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n_nodes; i++) {
        if (std::abs(nodes[i].y - 0.0f) < 1.0e-6f) {
            is_fixed[i] = true;
            u_vec[i] = 0.0f;
            dirichlet_row(i, A, F_vec, 0.0f);
        } else if (std::abs(nodes[i].y - 2.0f) < 1.0e-6f) {
            is_fixed[i] = true;
            u_vec[i] = 1.0f;
            dirichlet_row(i, A, F_vec, 1.0f);
        }
    }
    auto e2 = std::chrono::high_resolution_clock::now();

    // Jacobi前処理用対角成分
    std::vector<float> diag(n_nodes, 0.0f);
    for (int i = 0; i < n_nodes; i++) {
        for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; k++) {
            if (A.col_idx[k] == i) {
                diag[i] = A.values[k];
                break;
            }
        }
    }

    // Step6: PCGソルバー
    auto f1 = std::chrono::high_resolution_clock::now();
    pcg_solver(A, F_vec, u_vec, diag);
    auto f2 = std::chrono::high_resolution_clock::now();

    // Step7: 計算時間の出力
    auto t_total_end = std::chrono::high_resolution_clock::now();

    std::cout << "generation of nodes and mesh time:              "
              << std::chrono::duration<double>(a2-a1).count() << " seconds\n";
    std::cout << "Calculation of Element Stiffness Matrices time: "
              << std::chrono::duration<double>(c2-c1).count() << " seconds\n";
    std::cout << "Calculation of Global Stiffness Matrix time:    "
              << std::chrono::duration<double>(d2-d1).count() << " seconds\n";
    std::cout << "Application of Boundary Conditions time:        "
              << std::chrono::duration<double>(e2-e1).count() << " seconds\n";
    std::cout << "CG Method Wall Time:                            "
              << std::chrono::duration<double>(f2-f1).count() << " seconds\n";
    std::cout << "Total Time:                                     "
              << std::chrono::duration<double>(t_total_end-t_total_start).count() << " seconds\n";

    return 0;
}