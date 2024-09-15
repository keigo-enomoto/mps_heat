


# Strategy

## 入力ファイル
### パラメータのインプット
[数値計算向けパラメータファイルを読み込むシングルファイル C++ ライブラリ](https://qiita.com/kaityo256/items/88a30fe1ecd9abc8449a) を利用し、`param` クラスを作成することにする。
ちなみに、このクラスは後からも変数を追加できることに注意する。


## Eigen
行列計算は[Eigen](https://kawakamik.com/post/eigen/)に一任する
Eigenの内部でOpenMP並列化を実装しているので、for文の中にEigenの操作を入れるとややこしくなってしまうことに注意

- 圧力の陰解法
- 速度の陰解法

# Reference
- [SPH Codes](https://www.spheric-sph.org/sph-projects-and-codes)
- [SPH simulation](https://github.com/takagi/blog-codes/blob/master/20091101/cpp/sph.cpp)
    `vec3` を上手に使えるといいかもしれない（各次元の計算を一気に終わらせられる）