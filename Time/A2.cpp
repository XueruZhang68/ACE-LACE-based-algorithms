#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>

int gcd(int a, int b) {
    while (b != 0) {
        int temp = b;
        b = a % b;
        a = temp;
    }
    return std::abs(a);
}

std::vector<int> Euler(int N) {
    std::vector<int> t;
    for (int j = 1; j < N; ++j) {
        if (gcd(j, N) == 1) {
            t.push_back(j);
        }
    }
    return t;
}

std::vector<std::vector<int>> GLP(int N, const std::vector<int>& hm) {
    std::vector<int> h(N);
    for (int i = 0; i < N; ++i) {
        h[i] = i + 1;
    }
    std::vector<std::vector<int>> udt(N, std::vector<int>(hm.size(), 0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < hm.size(); ++j) {
            udt[i][j] = (h[i] * hm[j]) % N;
            if (udt[i][j] == 0) {
                udt[i][j] = N;
            }
        }
    }
    return udt;
}

std::vector<std::vector<int>> level(const std::vector<std::vector<int>>& D, int l) {
    int nrow = D.size();
    if (nrow == 0) return {};
    int ncol = D[0].size();
    std::vector<std::vector<int>> y(nrow, std::vector<int>(ncol, 0));
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            y[i][j] = (D[i][j] + l) % nrow;
            if (y[i][j] == 0) {
                y[i][j] = nrow;
            }
        }
    }
    return y;
}

// level_leave_one_out_com 函数
std::vector<std::vector<int>> level_leave_one_out_com(const std::vector<std::vector<int>>& D, const std::vector<int>& com) {
    int n_rows = D.size();
    if (n_rows == 0) return {};
    int n_cols_D = D[0].size(); // D 的列数
    int n_cols_com = com.size(); // com 的长度
    
    // 初始化 D1，每列来自 level(D, com[j]) 的结果，列数为 n_cols_D * n_cols_com
    std::vector<std::vector<int>> D1(n_rows, std::vector<int>(n_cols_D * n_cols_com));
    
    // 遍历 com 的每个元素
    for (int j = 0; j < n_cols_com; ++j) {
        // 对整个矩阵 D 应用 level 函数
        std::vector<std::vector<int>> level_result = level(D, com[j]);
        
        // 将 level_result 中等于 n_rows 的元素替换为 com[j]
        for (int i = 0; i < n_rows; ++i) {
            for (int k = 0; k < n_cols_D; ++k) {
                if (level_result[i][k] == n_rows) {
                    level_result[i][k] = com[j]; // 替换为 com[j]
                }
            }
        }
        
        // 将修改后的 level_result 填入 D1 的对应位置
        for (int i = 0; i < n_rows; ++i) {
            for (int k = 0; k < n_cols_D; ++k) {
                D1[i][j * n_cols_D + k] = level_result[i][k];
            }
        }
    }
    
    // 移除最后一行
    D1.pop_back();
    
    return D1;
}

// L_2
double L_2(const std::vector<std::vector<int>>& D, int p) {
    int N = D.size();
    std::vector<double> min_values;
    for (int x = 0; x < N - 2; ++x) {
        double min_dist = -1;
        for (int i = x + 1; i < N; ++i) {
            double sum = 0;
            for (int j = 0; j < D[0].size(); ++j) {
                sum += std::pow(std::abs(D[i][j] - D[x][j]), p);
            }
            if (min_dist < 0 || sum < min_dist) {
                min_dist = sum;
            }
        }
        min_values.push_back(min_dist);
    }
    double last_dist = 0;
    for (int j = 0; j < D[0].size(); ++j) {
        last_dist += std::pow(std::abs(D[N-1][j] - D[N-2][j]), p);
    }
    min_values.push_back(last_dist);

    double min_val = min_values[0];
    for (double val : min_values) {
        if (val < min_val) {
            min_val = val;
        }
    }
    return min_val;
}

// SA_leave_one_out_half_naive 函数
std::vector<std::vector<int>> SA_leave_one_out_half_naive(int N, int S, double T, int p, int total, double r) {
    // 生成 com = 1:((N+1)/2)-1
    std::vector<int> com((N + 1) / 2 - 1);
    for (int i = 0; i < com.size(); ++i) {
        com[i] = i; // 0-based 索引，从 0 到 (N+1)/2-2
    }

    // 生成 D
    std::vector<std::vector<int>> D = level_leave_one_out_com(GLP(N, Euler(N)), com);
    int n_rows = D.size();
    int n_cols = D[0].size();

    // 随机数生成器
    std::random_device rd;
    std::mt19937 gen(rd());

    // 随机选择初始 column_optimal
    std::vector<int> indices(n_cols);
    for (int i = 0; i < n_cols; ++i) indices[i] = i;
    std::shuffle(indices.begin(), indices.end(), gen);
    std::vector<int> column_optimal(S);
    for (int i = 0; i < S; ++i) column_optimal[i] = indices[i];

    // 计算初始 value_optimal
    std::vector<std::vector<int>> D_optimal(n_rows, std::vector<int>(S));
    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < S; ++j) {
            D_optimal[i][j] = D[i][column_optimal[j]];
        }
    }
    double value_optimal = L_2(D_optimal, p) / std::floor(std::pow((N - 1) / 2.0, p - 1) * N * S / 3.0);

    // 模拟退火迭代
    for (int i = 0; i < total; ++i) {
        // 获取剩余列
        std::vector<int> remaining_cols;
        for (int col = 0; col < n_cols; ++col) {
            if (std::find(column_optimal.begin(), column_optimal.end(), col) == column_optimal.end()) {
                remaining_cols.push_back(col);
            }
        }

        // 随机选择一个剩余列
        std::uniform_int_distribution<> dis_remain(0, remaining_cols.size() - 1);
        int try_column = remaining_cols[dis_remain(gen)];

        // 随机替换 column_optimal 中的一列
        std::vector<int> column_try = column_optimal;
        std::uniform_int_distribution<> dis_optimal(0, S - 1);
        int replace_idx = dis_optimal(gen);
        column_try[replace_idx] = try_column;

        // 计算 value_try
        std::vector<std::vector<int>> D_try(n_rows, std::vector<int>(S));
        for (int r = 0; r < n_rows; ++r) {
            for (int c = 0; c < S; ++c) {
                D_try[r][c] = D[r][column_try[c]];
            }
        }
        double value_try = L_2(D_try, p) / std::floor(std::pow((N - 1) / 2.0, p - 1) * N * S / 3.0);

        // 模拟退火决策
        if (value_try > value_optimal) {
            column_optimal = column_try;
            value_optimal = value_try;
        } else {
            std::uniform_real_distribution<> dis_prob(0, 1);
            double prob = std::exp((value_try - value_optimal) / T);
            if (prob > dis_prob(gen)) {
                column_optimal = column_try;
                value_optimal = value_try;
            }
        }

        // 更新温度
        T *= r;
    }

    // 返回优化后的子矩阵
    std::vector<std::vector<int>> result(n_rows, std::vector<int>(S));
    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < S; ++j) {
            result[i][j] = D[i][column_optimal[j]];
        }
    }
    return result;
}
