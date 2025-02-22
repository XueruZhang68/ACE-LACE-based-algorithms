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


std::vector<std::vector<int>> construction(int SetN) {
    std::vector<int> aa = Euler(SetN); // aa的大小为φ(SetN)
    std::vector<std::vector<int>> D = GLP(SetN, aa); // D为SetN × φ(SetN)
    int ncol = aa.size(); // D的列数
    // D1的列数变为SetN × ncol
    std::vector<std::vector<int>> D1(SetN, std::vector<int>(SetN * ncol, 0));
    
    for (int j = 0; j < SetN; ++j) {
        std::vector<std::vector<int>> leveled = level(D, j); // leveled为SetN × ncol
        for (int i = 0; i < SetN; ++i) {
            for (int k = 0; k < ncol; ++k) {
                D1[i][j * ncol + k] = leveled[i][k]; // 取全部列
            }
        }
    }
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
