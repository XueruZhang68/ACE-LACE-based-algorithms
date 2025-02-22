#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <chrono>

#define Min(a, b) ((a < b) ? a : b)
#define Max(a, b) ((a > b) ? a : b)

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

// 列合并函数（替代 R 的 cbind）
std::vector<std::vector<int>> column_bind(const std::vector<std::vector<int>>& D1, const std::vector<std::vector<int>>& D2) {
    int n_rows = D1.size();
    int n_cols1 = D1[0].size();
    int n_cols2 = D2[0].size();
    std::vector<std::vector<int>> result(n_rows, std::vector<int>(n_cols1 + n_cols2));
    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols1; ++j) {
            result[i][j] = D1[i][j];
        }
        for (int j = 0; j < n_cols2; ++j) {
            result[i][n_cols1 + j] = D2[i][j];
        }
    }
    return result;
}

// Translated construction_leave_one_out function
std::vector<std::vector<int>> construction_leave_one_out(int SetN) {
    std::vector<int> aa = Euler(SetN); // aa的大小为φ(SetN)
    std::vector<std::vector<int>> D = GLP(SetN, aa); // D为SetN × φ(SetN)
    int ncol = aa.size(); // D的列数
    // D1的列数变为SetN × ncol
    std::vector<std::vector<int>> D1(SetN, std::vector<int>(SetN * ncol, 0));
    
    for (int j = 0; j < SetN; ++j) {
        std::vector<std::vector<int>> leveled = level(D, j); // leveled为SetN × ncol
         // 将 level_result 中等于 n_rows 的元素替换为 com[j]
        for (int i = 0; i < SetN; ++i) {
            for (int k = 0; k < ncol; ++k) {
                if (leveled[i][k] == SetN) {
                    leveled[i][k] = j; // 替换为 com[j]
                }
            }
        }
        for (int i = 0; i < SetN; ++i) {
            for (int k = 0; k < ncol; ++k) {
                D1[i][j * ncol + k] = leveled[i][k]; // 取全部列
            }
        }
    }
    // Remove the last row (index SetN-1 in 0-based indexing)
    D1.pop_back();
    return D1;
}

std::vector<std::vector<int>> subset_columns(const std::vector<std::vector<int>>& D, const std::vector<int>& columns) {
    int n_rows = D.size();
    int n_cols = columns.size();
    std::vector<std::vector<int>> result(n_rows, std::vector<int>(n_cols));
    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            result[i][j] = D[i][columns[j]];
        }
    }
    return result;
}

// SA_leave_one_out_naive function
std::vector<std::vector<int>> SA_leave_one_out_naive(int N, int S, double T, int p, int total, double r) {
    std::vector<int> aa = Euler(N);
    int total_columns = N * aa.size();

    if (S < total_columns) {
        // Compute D
        std::vector<std::vector<int>> D = construction_leave_one_out(N);

        // Random number generation setup
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        // Initial column selection
        std::vector<int> all_columns(total_columns);
        for (int i = 0; i < total_columns; ++i) {
            all_columns[i] = i;
        }
        std::shuffle(all_columns.begin(), all_columns.end(), gen);
        std::vector<int> column_optimal(all_columns.begin(), all_columns.begin() + S);

        // Initial value
        std::vector<std::vector<int>> D_optimal = subset_columns(D, column_optimal);
        double normalization = std::floor(std::pow((N - 1) / 2.0, p - 1) * N * S / 3.0);
        double value_optimal = L_2(D_optimal, p) / normalization;

        // Simulated annealing loop
        for (int i = 0; i < total; ++i) {
            // Find remaining columns
            std::vector<int> remaining_columns;
            for (int j = 0; j < total_columns; ++j) {
                if (std::find(column_optimal.begin(), column_optimal.end(), j) == column_optimal.end()) {
                    remaining_columns.push_back(j);
                }
            }

            // Sample one column from remaining
            std::uniform_int_distribution<> dis_remaining(0, remaining_columns.size() - 1);
            int try_column = remaining_columns[dis_remaining(gen)];

            // Create column_try by replacing one column
            std::vector<int> column_try = column_optimal;
            std::uniform_int_distribution<> dis_optimal(0, S - 1);
            column_try[dis_optimal(gen)] = try_column;

            // Compute new value
            std::vector<std::vector<int>> D_try = subset_columns(D, column_try);
            double value_try = L_2(D_try, p) / normalization;

            // Acceptance criterion
            if (value_try > value_optimal) {
                column_optimal = column_try;
                value_optimal = value_try;
            } else {
                double prob = std::exp((value_try - value_optimal) / T);
                if (prob > dis(gen)) {
                    column_optimal = column_try;
                    value_optimal = value_try;
                }
            }

            // Update temperature
            T *= r;
        }

        return subset_columns(D, column_optimal);
    } else {
        return construction_leave_one_out(N);
    }
}


// New function: SA_leave_one_out_naive_combined
std::vector<std::vector<int>> SA_leave_one_out_naive_combined(int N, int S, double T, int p, int total, double r, int nn, std::vector<int> column_optimal) {
    std::vector<int> aa = Euler(N);
    int total_columns = N * aa.size();

    if (S < total_columns) {
        // Compute D
        std::vector<std::vector<int>> D = construction_leave_one_out(N);

        // Random number generation setup
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        // Initial value using provided column_optimal
        std::vector<std::vector<int>> D_optimal = subset_columns(D, column_optimal);
        double normalization = std::floor(std::pow((N - 1) / 2.0, p - 1) * N * S / 3.0);
        double value_optimal = L_2(D_optimal, p) / normalization;

        // Simulated annealing loop starting from nn
        for (int i = nn; i <= total; ++i) {
            // Find remaining columns
            std::vector<int> remaining_columns;
            for (int j = 0; j < total_columns; ++j) {
                if (std::find(column_optimal.begin(), column_optimal.end(), j) == column_optimal.end()) {
                    remaining_columns.push_back(j);
                }
            }

            // Sample one column from remaining
            std::uniform_int_distribution<> dis_remaining(0, remaining_columns.size() - 1);
            int try_column = remaining_columns[dis_remaining(gen)];

            // Create column_try by replacing one column
            std::vector<int> column_try = column_optimal;
            std::uniform_int_distribution<> dis_optimal(0, S - 1);
            column_try[dis_optimal(gen)] = try_column;

            // Compute new value
            std::vector<std::vector<int>> D_try = subset_columns(D, column_try);
            double value_try = L_2(D_try, p) / normalization;

            // Acceptance criterion
            if (value_try > value_optimal) {
                column_optimal = column_try;
                value_optimal = value_try;
            } else {
                double prob = std::exp((value_try - value_optimal) / T);
                if (prob > dis(gen)) {
                    column_optimal = column_try;
                    value_optimal = value_try;
                }
            }

            // Update temperature
            T *= r;
        }

        return subset_columns(D, column_optimal);
    } else {
        return construction_leave_one_out(N);
    }
}

// Helper function to compute variance
double variance(const std::vector<double>& values) {
    double mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
    double sum_sq_diff = 0.0;
    for (double value : values) {
        sum_sq_diff += (value - mean) * (value - mean);
    }
    return sum_sq_diff / (values.size() - 1);
}

// New function: SA_leave_one_out_ACE_combined
std::vector<std::vector<int>> SA_leave_one_out_ACE_combined(int N, int S, double T, int p, int total, double r, int step, double delta) {
    std::vector<int> aa = Euler(N);
    int phi_N = aa.size();

    if (S < phi_N) {
        return SA_leave_one_out_naive(N, S, T, p, total, r);
    } else {
        int m = S / phi_N;
        int K = S % phi_N;
        std::vector<std::vector<int>> D = GLP(N, aa);

        if (m < N) {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(0.0, 1.0);

            // Step 2: Initial U_optimal
            std::vector<int> U_0(N);
            for (int i = 0; i < N; ++i) {
                U_0[i] = i + 1; // 1-based indices
            }
            std::shuffle(U_0.begin(), U_0.end(), gen);
            std::vector<int> U_optimal(U_0.begin(), U_0.begin() + m);

            // Initial column_optimal
            std::vector<int> all_columns((N - m) * phi_N);
            for (int i = 0; i < (N - m) * phi_N; ++i) {
                all_columns[i] = i;
            }
            std::shuffle(all_columns.begin(), all_columns.end(), gen);
            std::vector<int> column_optimal(all_columns.begin(), all_columns.begin() + K);

            // Compute remaining indices
            std::vector<int> U_remaining;
            for (int i = 0; i < N; ++i) {
                if (std::find(U_optimal.begin(), U_optimal.end(), U_0[i]) == U_optimal.end()) {
                    U_remaining.push_back(U_0[i]);
                }
            }

            // Initial value_optimal
            std::vector<int> U_optimal_minus_1 = U_optimal;
            for (int& u : U_optimal_minus_1) u -= 1; // Convert to 0-based
            std::vector<int> U_remaining_minus_1 = U_remaining;
            for (int& u : U_remaining_minus_1) u -= 1;
            std::vector<std::vector<int>> D1 = level_leave_one_out_com(D, U_optimal_minus_1);
            std::vector<std::vector<int>> D2 = subset_columns(level_leave_one_out_com(D, U_remaining_minus_1), column_optimal);
            std::vector<std::vector<int>> Target = column_bind(D1, D2);
            double normalization = std::floor(std::pow((N - 1) / 2.0, p - 1) * N * S / 3.0);
            double value_optimal = L_2(Target, p) / normalization;

            // Value history
            std::vector<double> value(step, value_optimal);
            int iiii = 0;

            // SA loop
            while (variance(value) > delta && iiii <= total) {
                // Perturb U_try
                std::vector<int> U_try = U_optimal;
                std::uniform_int_distribution<> dis_optimal(0, m - 1);
                std::uniform_int_distribution<> dis_remaining(0, U_remaining.size() - 1);
                int idx_opt = dis_optimal(gen);
                int idx_rem = dis_remaining(gen);
                U_try[idx_opt] = U_remaining[idx_rem];

                // New column_try
                std::shuffle(all_columns.begin(), all_columns.end(), gen);
                std::vector<int> column_try(all_columns.begin(), all_columns.begin() + K);

                // Compute value_try
                std::vector<int> U_try_minus_1 = U_try;
                for (int& u : U_try_minus_1) u -= 1;
                std::vector<int> U_remaining_try;
                for (int i = 0; i < N; ++i) {
                    if (std::find(U_try.begin(), U_try.end(), U_0[i]) == U_try.end()) {
                        U_remaining_try.push_back(U_0[i]);
                    }
                }
                std::vector<int> U_remaining_try_minus_1 = U_remaining_try;
                for (int& u : U_remaining_try_minus_1) u -= 1;
                D1 = level_leave_one_out_com(D, U_try_minus_1);
                D2 = subset_columns(level_leave_one_out_com(D, U_remaining_try_minus_1), column_try);
                std::vector<std::vector<int>> Target_try = column_bind(D1, D2);
                double value_try = L_2(Target_try, p) / normalization;

                // Acceptance criterion
                if (value_try > value_optimal) {
                    column_optimal = column_try;
                    U_optimal = U_try;
                    U_remaining = U_remaining_try;
                    value_optimal = value_try;
                    Target = Target_try;
                } else {
                    double prob = std::exp((value_try - value_optimal) / T);
                    if (prob > dis(gen)) {
                        column_optimal = column_try;
                        U_optimal = U_try;
                        U_remaining = U_remaining_try;
                        value_optimal = value_try;
                        Target = Target_try;
                    }
                }

                // Update value history
                value.erase(value.begin());
                value.push_back(value_optimal);
                T *= r;
                iiii++;
            }

            // Early termination handling
            if (iiii < total) {
                std::vector<int> new_column_optimal;
                for (int x : U_optimal) {
                    int start = (x - 1) * phi_N;
                    for (int j = start; j < start + phi_N; ++j) {
                        new_column_optimal.push_back(j);
                    }
                }
                if (!column_optimal.empty()) {
                    for (int x : column_optimal) {
                        int col = x % phi_N;
                        if (col == 0) col = phi_N;
                        int row_idx = (x - 1) / phi_N + 1;
                        int h = col + (U_remaining[row_idx] - 1) * phi_N - 1;
                        new_column_optimal.push_back(h);
                    }
                }
                return SA_leave_one_out_naive_combined(N, S, T, p, total, r, iiii + 1, new_column_optimal);
            }
            return Target;
        } else {
            return construction_leave_one_out(N);
        }
    }
}

