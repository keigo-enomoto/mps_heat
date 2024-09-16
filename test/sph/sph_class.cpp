// https://chatgpt.com/share/a9268c69-8a5f-440f-bf35-0ef4c4f4d62f

#include <vector>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <tuple>
#include <fstream>
#include <iomanip>

// ヘルパークラス（座標の整数化）
struct Int3 {
    int x, y, z;
    bool operator==(const Int3 &other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

// Int3のハッシュ関数
namespace std {
    template <>
    struct hash<Int3> {
        std::size_t operator()(const Int3& k) const {
            return ((std::hash<int>()(k.x)
                    ^ (std::hash<int>()(k.y) << 1)) >> 1)
                    ^ (std::hash<int>()(k.z) << 1);
        }
    };
}

// 定数
const double h = 0.1;                // スムージング長
const double cell_size = 2.0 * h;    // セルサイズ（2倍のスムージング長）
const double rho0 = 1000.0;          // 基準密度
const double k = 1000.0;             // 圧力係数
const double alpha = 0.01;           // 熱拡散係数
const double Tm = 0.5;               // 融点温度
const double Lf = 0.1;               // 潜熱
const double cp = 1.0;               // 比熱
const double epsilon = 0.01;         // 粘性項の補正定数
const double c = 10.0;               // 音速
const double alpha_viscosity = 0.1;  // 人工粘性係数
const double k_solid = 10000.0;      // 固体粒子間の弾性係数

// 粒子の状態を表す列挙型
enum class ParticleState {
    Liquid,
    Solid
};

// 粒子クラス
class Particle {
public:
    ParticleState state;
    std::vector<double> r;  // 位置
    std::vector<double> v;  // 速度
    double mass;
    double rho;
    double p;
    double T;
    double H;

    // コンストラクタ (Initialize)
    Particle() : mass(1.0), rho(rho0), p(0), T(0), H(0), state(ParticleState::Liquid), v(3, 0.0), r(3, 0.0) {}

    // 粒子の圧力を計算
    void computePressure() {
        p = k * (rho - rho0);
    }

    // 粒子のエンタルピーと状態を更新
    void updateEnthalpyAndState() {
        if (T <= Tm) {
            H = cp * T;
            if (state == ParticleState::Liquid) {
                state = ParticleState::Solid;  // 固化
            }
        } else {
            H = cp * Tm + Lf;
        }
    }
};

// SPHシミュレーションクラス
class SPHSimulation {
public:
    std::vector<Particle> particles;

    // コンストラクタ：長方形の初期配置を作成
    SPHSimulation(int num_particles_x, int num_particles_y, int num_particles_z, double spacing) {
        // 初期化
        particles.reserve(num_particles_x * num_particles_y * num_particles_z);

        // 格子状に粒子を配置
        for (int i = 0; i < num_particles_x; ++i) {
            for (int j = 0; j < num_particles_y; ++j) {
                for (int k = 0; k < num_particles_z; ++k) {
                    Particle p;
                    p.r = {i * spacing, j * spacing, k * spacing}; // 位置を格子状に設定
                    particles.push_back(p);
                }
            }
        }
        // 追加: 粒子数が正しいか確認するための出力
        std::cout << "Initialized " << particles.size() << " particles." << std::endl;

    }

    // シミュレーションを実行
    void simulate(double dt, int steps) {
        for (int step = 0; step < steps; ++step) {
            // 近傍リストの構築
            auto neighbor_list = buildNeighborList();

            // 各種物理量の更新
            updateDensity(neighbor_list);
            updatePressure();
            updateVelocity(neighbor_list, dt);
            updateTemperature(neighbor_list, dt);
            updateEnthalpyAndState();

            // 位置の更新
            updatePosition(dt);

            if(step%100 == 0){
                outputPhysicalQuantities(step);
            }
        }
    }

private:
    // セルに基づく近傍リストの構築
    std::vector<std::vector<int>> buildNeighborList() {
        std::unordered_map<Int3, std::vector<int>> cell_map;

        // 粒子をセルに分配
        for (int i = 0; i < particles.size(); ++i) {
            Int3 cell = getCell(particles[i]);
            cell_map[cell].push_back(i);
        }

        // 近傍リストの構築
        std::vector<std::vector<int>> neighbor_list(particles.size());

        for (int i = 0; i < particles.size(); ++i) {
            Int3 cell = getCell(particles[i]);

            // 近傍セルを探索
            for (int x = -1; x <= 1; ++x) {
                for (int y = -1; y <= 1; ++y) {
                    for (int z = -1; z <= 1; ++z) {
                        Int3 neighbor_cell = { cell.x + x, cell.y + y, cell.z + z };
                        if (cell_map.find(neighbor_cell) != cell_map.end()) {
                            for (int j : cell_map[neighbor_cell]) {
                                if (i != j) {
                                    neighbor_list[i].push_back(j);
                                }
                            }
                        }
                    }
                }
            }
        }

        return neighbor_list;
    }

    // 質量保存則
    void updateDensity(const std::vector<std::vector<int>>& neighbor_list) {
        for (int i = 0; i < particles.size(); ++i) {
            Particle& p_i = particles[i];
            p_i.rho = 0.0;
            for (int j : neighbor_list[i]) {
                const Particle& p_j = particles[j];
                auto r_ij = distance(p_i.r, p_j.r);
                p_i.rho += p_j.mass * W(r_ij);
            }
        }
    }

    // 粒子の圧力を更新
    void updatePressure() {
        for (auto& p : particles) {
            p.computePressure();
        }
    }

    // 運動量保存則
    void updateVelocity(const std::vector<std::vector<int>>& neighbor_list, double dt) {
        for (int i = 0; i < particles.size(); ++i) {
            Particle& p_i = particles[i];
            std::vector<double> dv = {0.0, 0.0, 0.0};
            for (int j : neighbor_list[i]) {
                const Particle& p_j = particles[j];
                if (&p_i == &p_j) continue;

                auto r_ij = distance(p_i.r, p_j.r);
                auto gradW_ij = gradW(r_ij);

                double pressure_term = (p_i.p / (p_i.rho * p_i.rho)) + (p_j.p / (p_j.rho * p_j.rho));
                double viscosity_term = artificialViscosity(p_i, p_j, r_ij);

                for (int d = 0; d < 3; ++d) {
                    dv[d] += -p_j.mass * (pressure_term + viscosity_term) * gradW_ij[d];
                }

                // 固体粒子同士の弾性相互作用
                if (p_i.state == ParticleState::Solid && p_j.state == ParticleState::Solid) {
                    auto elastic_force = computeElasticForce(p_i, p_j, r_ij);
                    for (int d = 0; d < 3; ++d) {
                        dv[d] += elastic_force[d];
                    }
                }
            }
            for (int d = 0; d < 3; ++d) {
                p_i.v[d] += dv[d] * dt;
            }
        }
    }

    // 熱伝導方程式の更新
    void updateTemperature(const std::vector<std::vector<int>>& neighbor_list, double dt) {
        for (int i = 0; i < particles.size(); ++i) {
            Particle& p_i = particles[i];
            double dT = 0.0;
            for (int j : neighbor_list[i]) {
                const Particle& p_j = particles[j];
                auto r_ij = distance(p_i.r, p_j.r);
                dT += alpha * (p_j.T - p_i.T) * W(r_ij) / p_j.rho;
            }
            p_i.T += dT * dt;
        }
    }

    // エンタルピーと粒子の状態を更新
    void updateEnthalpyAndState() {
        for (auto& p : particles) {
            p.updateEnthalpyAndState();
        }
    }

    // 粒子の位置を更新
    void updatePosition(double dt) {
        for (auto& p : particles) {
            for (int d = 0; d < 3; ++d) {
                p.r[d] += p.v[d] * dt;
            }
        }
    }

    // セルを取得
    Int3 getCell(const Particle& p) {
        return {
            static_cast<int>(std::floor(p.r[0] / cell_size)),
            static_cast<int>(std::floor(p.r[1] / cell_size)),
            static_cast<int>(std::floor(p.r[2] / cell_size))
        };
    }

    // カーネル関数（Poly6カーネル）
    double W(const std::vector<double>& r) {
        double r_len = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        if (r_len >= 0 && r_len <= h) {
            double q = r_len / h;
            return (315.0 / (64.0 * M_PI * pow(h, 9))) * pow(h*h - r_len*r_len, 3);
        } else {
            return 0.0;
        }
    }

    // カーネルの勾配
    std::vector<double> gradW(const std::vector<double>& r) {
        double r_len = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        if (r_len > 0 && r_len <= h) {
            double coef = -945.0 / (32.0 * M_PI * pow(h, 9)) * pow(h*h - r_len*r_len, 2);
            return { coef * r[0] / r_len, coef * r[1] / r_len, coef * r[2] / r_len };
        } else {
            return {0.0, 0.0, 0.0};
        }
    }

    // 粘性項の計算
    double artificialViscosity(const Particle& p_i, const Particle& p_j, const std::vector<double>& r_ij) {
        double r_len = sqrt(r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2]);
        if (r_len == 0.0) return 0.0;

        std::vector<double> v_ij = { p_i.v[0] - p_j.v[0], p_i.v[1] - p_j.v[1], p_i.v[2] - p_j.v[2] };
        double dot_product = v_ij[0]*r_ij[0] + v_ij[1]*r_ij[1] + v_ij[2]*r_ij[2];

        if (dot_product < 0) {
            return (-alpha_viscosity * c * h * dot_product) / (p_i.rho + p_j.rho + epsilon * h * h);
        }
        return 0.0;
    }

    // 弾性力の計算
    std::vector<double> computeElasticForce(const Particle& p_i, const Particle& p_j, const std::vector<double>& r_ij) {
        double r_len = sqrt(r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2]);
        double k_effective = k_solid * (r_len - h);
        if (r_len < h) {
            return { -k_effective * r_ij[0] / r_len, -k_effective * r_ij[1] / r_len, -k_effective * r_ij[2] / r_len };
        } else {
            return {0.0, 0.0, 0.0};
        }
    }

    // 粒子間の距離ベクトルを計算
    std::vector<double> distance(const std::vector<double>& r_i, const std::vector<double>& r_j) {
        return { r_j[0] - r_i[0], r_j[1] - r_i[1], r_j[2] - r_i[2] };
    }

    // 物理量をファイルに出力する関数
    void outputPhysicalQuantities(int step) {
        std::ofstream outfile;
        std::string filename = "output_step_" + std::to_string(step) + ".txt";
        outfile.open(filename);

        if (!outfile.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing.\n";
            return;
        }

        outfile << std::setprecision(5);  // 精度を指定（小数点以下5桁）
        outfile << "Step: " << step << "\n";
        outfile << "Index\tPosition\t\tVelocity\t\tDensity\tPressure\tTemperature\tEnthalpy\tState\n";
        for (int i = 0; i < particles.size(); ++i) {
            Particle& p = particles[i];
            outfile << i << "\t"
                    << "(" << p.r[0] << ", " << p.r[1] << ", " << p.r[2] << ")\t"
                    << "(" << p.v[0] << ", " << p.v[1] << ", " << p.v[2] << ")\t"
                    << p.rho << "\t"
                    << p.p << "\t"
                    << p.T << "\t"
                    << p.H << "\t"
                    << (p.state == ParticleState::Liquid ? "Liquid" : "Solid") << "\n";
        }

        outfile.close();
    }

};

// メイン関数
int main() {
    // 粒子数100で初期化
    SPHSimulation simulation(10,10,1,1);
    double dt = 0.01; // 時間刻み
    int steps = 1000; // シミュレーションステップ数

    // シミュレーション実行
    simulation.simulate(dt, steps);

    std::cout << "Simulation completed." << std::endl;

    return 0;
}
