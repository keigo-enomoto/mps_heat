#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
// #include <omp.h>

// ベクトルの内積を計算
double dot_product(int n, const std::vector<double>& x, const std::vector<double>& y) {
    double result = 0.0;
    #pragma omp parallel for reduction(+:result)
    for (int i = 0; i < n; i++) {
        result += x[i] * y[i];
    }
    return result;
}

// ベクトルの線形結合 y = y + alpha * x
void axpy(int n, double alpha, const std::vector<double>& x, std::vector<double>& y) {
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        y[i] += alpha * x[i];
    }
}

// ベクトルのスカラー倍 x = alpha * x
void scale(int n, double alpha, std::vector<double>& x) {
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        x[i] *= alpha;
    }
}

// Incomplete Cholesky分解を実行
void incomplete_cholesky(int n, const std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& L) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = A[i][j];
            for (int k = 0; k < j; k++) {
                sum -= L[i][k] * L[j][k];
            }
            if (i == j) {
                L[i][j] = std::sqrt(sum);
            } else {
                L[i][j] = sum / L[j][j];
            }
        }
    }
}

// ICCG法を実行
void iccg(int n, const std::vector<std::vector<double>>& A, const std::vector<double>& b, std::vector<double>& x, int max_iter, double tol) {
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));

    incomplete_cholesky(n, A, L);

    std::vector<double> r(n), z(n), p(n), q(n);

    // 初期残差計算 r = b - Ax
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        r[i] = b[i] - sum;
    }

    // z = L^(-1) * r
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        double sum = r[i];
        for (int j = 0; j < i; j++) {
            sum -= L[i][j] * z[j];
        }
        z[i] = sum / L[i][i];
    }

    #pragma omp parallel for
    for (int i = n - 1; i >= 0; i--) {
        double sum = z[i];
        for (int j = i + 1; j < n; j++) {
            sum -= L[j][i] * z[j];
        }
        z[i] = sum / L[i][i];
    }

    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        p[i] = z[i];
    }

    double rz_old = dot_product(n, r, z);

    for (int k = 0; k < max_iter; k++) {
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += A[i][j] * p[j];
            }
            q[i] = sum;
        }

        double alpha = rz_old / dot_product(n, p, q);
        axpy(n, alpha, p, x);
        axpy(n, -alpha, q, r);

        double norm_r = std::sqrt(dot_product(n, r, r));
        if (norm_r < tol) {
            break;
        }

        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            double sum = r[i];
            for (int j = 0; j < i; j++) {
                sum -= L[i][j] * z[j];
            }
            z[i] = sum / L[i][i];
        }

        #pragma omp parallel for
        for (int i = n - 1; i >= 0; i--) {
            double sum = z[i];
            for (int j = i + 1; j < n; j++) {
                sum -= L[j][i] * z[j];
            }
            z[i] = sum / L[i][i];
        }

        double rz_new = dot_product(n, r, z);
        double beta = rz_new / rz_old;
        rz_old = rz_new;

        scale(n, beta, p);
        axpy(n, 1.0, z, p);
    }
}

int main() {
        // 計算開始時刻を記録
    auto start = std::chrono::high_resolution_clock::now();

    int n = 4; // 行列のサイズ
    std::vector<double> b = {1, 2, 3, 4};
    std::vector<double> x(n, 0.0); // 初期解はゼロベクトル
    std::vector<std::vector<double>> A = {
        {4, 1, 0, 0},
        {1, 4, 1, 0},
        {0, 1, 4, 1},
        {0, 0, 1, 4}
    };

    int max_iter = 1000;
    double tol = 1e-6;

    iccg(n, A, b, x, max_iter, tol);

    // 計算終了時刻を記録
    auto end = std::chrono::high_resolution_clock::now();
    // 経過時間を計算（ミリ秒単位）
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "解ベクトル x:\n";
    for (int i = 0; i < n; i++) {
        std::cout << x[i] << std::endl;
    }
    std::cout << "計算時間: " << duration.count() << " ミリ秒" << std::endl;

    return 0;
}
