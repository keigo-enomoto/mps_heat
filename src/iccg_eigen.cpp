
// これでいけそう！！
// https://eigen.tuxfamily.org/dox/TopicMultiThreading.html
// ICCG is parallelized in default settings

#include <iostream>
#include <chrono>

// #include <Eigen/Sparse>
// #include <Eigen/IterativeLinearSolvers>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/IterativeLinearSolvers>

using namespace Eigen;
using namespace std::chrono;

int main() {

    // 計算開始時刻を記録
    auto start = high_resolution_clock::now();


    // 行列Aのサイズ
    int n = 4;

    // 正定値対称行列Aの初期化
    SparseMatrix<double> A(n, n);
    std::vector<Triplet<double>> tripletList;

    // 対角要素
    tripletList.push_back(Triplet<double>(0, 0, 4));
    tripletList.push_back(Triplet<double>(1, 1, 4));
    tripletList.push_back(Triplet<double>(2, 2, 4));
    tripletList.push_back(Triplet<double>(3, 3, 4));

    // 非対角要素 (対称行列なので両側に同じ値を入れる)
    tripletList.push_back(Triplet<double>(0, 1, 1));
    tripletList.push_back(Triplet<double>(1, 0, 1));

    tripletList.push_back(Triplet<double>(1, 2, 1));
    tripletList.push_back(Triplet<double>(2, 1, 1));

    tripletList.push_back(Triplet<double>(2, 3, 1));
    tripletList.push_back(Triplet<double>(3, 2, 1));

    A.setFromTriplets(tripletList.begin(), tripletList.end());

    // ベクトルbの初期化
    VectorXd b(n);
    b << 1, 2, 3, 4;

    // 解ベクトルxの初期化
    VectorXd x = VectorXd::Zero(n);

    // ICCG法の実行
    // ConjugateGradient<SparseMatrix<double>, Lower | Upper, IncompleteCholesky<double>> solver;
    ConjugateGradient<SparseMatrix<double>, Lower | Upper, SimplicialCholesky<SparseMatrix<double>>> solver;
    solver.compute(A);

    if (solver.info() != Success) {
        std::cerr << "Failed to decompose matrix A with Incomplete Cholesky." << std::endl;
        return -1;
    }

    x = solver.solve(b);

    if (solver.info() != Success) {
        std::cerr << "Failed to solve the system with ICCG method." << std::endl;
        return -1;
    }

    // 計算終了時刻を記録
    auto end = high_resolution_clock::now();
    // 経過時間を計算（ミリ秒単位）
    auto duration = duration_cast<milliseconds>(end - start);


    // 結果の出力
    std::cout << "解ベクトル x:\n" << x << std::endl;
    std::cout << "計算時間: " << duration.count() << " ミリ秒" << std::endl;

    VectorXd b_check(n);
    b_check = A * x;
    // 結果の出力
    std::cout << "右辺の確認 b_check = A * x:\n" << b_check << std::endl;

    return 0;
}

