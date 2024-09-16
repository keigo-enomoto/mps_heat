


# Strategy

## 流れ
- `param.hpp` の`param` クラスで諸々のパラメータを引き受ける
```
L=32
p=2.0
```
のような形式。これは最初だけ使うようにしたい。

- `Particle` 等のセッティング


## クラス

### `InputParams`
[kaityoさんの実装したクラス](https://qiita.com/kaityo256/items/88a30fe1ecd9abc8449a) をそのまま利用。

### `Parameters`
インプットパラメータの管理を行う。`InputParams`を使用して、ファイルから変数を読み込む

### `Particle`
個々の粒子の物理特性（位置、速度、加速度）を保持し、更新する

- メンバー変数 : 位置や速度、加速度などその粒子の物性および性質
- メンバー関数 : 各メンバー変数の更新

### `Simulation`
シミュレーション全体の管理をするクラス
- メンバー変数 : `Particle`のvector

### `PressureSolver`
  行列Aを入力としてもらい、EigenでICCG法を用いて行列の解を計算する。その後値を返す

### `CellList`
セルリスト用のパラメータの管理

- メンバー変数：インデックス周りの配列、メッシュの大きさ等の変数
- メンバー関数 : 上記の変数を受け取るための関数

### `Output`
出力用のクラス

### `Util`
その他の雑関数を集めたクラス



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
- [C++ Eigenの使い方](https://kawakamik.com/post/eigen/)