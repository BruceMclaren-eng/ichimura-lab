#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <numeric>

// x,y のペアを表す構造体
struct Pair {
    float x, y;
};

// CRS形式の構造体
struct CRSMatrix {
    int n_rows;
    int nnz;
    std::vector<float> values;
    std::vector<int>   col_idx;
    std::vector<int>   row_ptr;
};

int nx, ny;
int n_nodes, n_elems, n_edges;
int nnz_global;
std::vector<Pair>  nodes;
std::vector<std::vector<int>> elems; 


// ============================================================
// generate_nodes, elems
// ============================================================
void generate_nodes(int nx, int ny, // 分割数
                    float x_max, float y_max, // 領域の最大値
                    float x_min, float y_min, // 領域の最小値
                    std::vector<Pair>& nodes)
{
    float dx = (x_max - x_min) / nx;
    float dy = (y_max - y_min) / ny;

    int node_id = 0;
    for (int j = 0; j <= ny; j++) {
        for (int i = 0; i <= nx; i++)
    }
}
