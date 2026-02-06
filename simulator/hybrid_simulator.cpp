#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <functional>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <optional>
#include <queue>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "vedic_sutras_complete.hpp"

namespace simulator {

using vedic::BigInt;
using vedic::Rational;

constexpr long double kPi = 3.14159265358979323846264338327950288L;
constexpr long double kTwoPi = 2.0L * kPi;

struct SutraOutput {
    std::string name;
    Rational value;
    std::string detail;
    double angle;
};

struct SimulationModeResult {
    std::string mode;
    std::vector<SutraOutput> outputs;
    Rational classical_total;
    double quantum_expectation;
    std::vector<double> probabilities;
};

class ThreadPool {
public:
    explicit ThreadPool(size_t workers) : stop_(false) {
        if (workers == 0) {
            workers = 1;
        }
        for (size_t i = 0; i < workers; ++i) {
            threads_.emplace_back([this]() { this->worker_loop(); });
        }
    }

    ~ThreadPool() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stop_ = true;
        }
        cv_.notify_all();
        for (auto &t : threads_) {
            if (t.joinable()) {
                t.join();
            }
        }
    }

    template <typename Func>
    auto submit(Func f) -> std::future<decltype(f())> {
        using Result = decltype(f());
        auto task = std::make_shared<std::packaged_task<Result()>>(std::move(f));
        auto future = task->get_future();
        {
            std::lock_guard<std::mutex> lock(mutex_);
            queue_.push([task]() { (*task)(); });
        }
        cv_.notify_one();
        return future;
    }

private:
    void worker_loop() {
        while (true) {
            std::function<void()> job;
            {
                std::unique_lock<std::mutex> lock(mutex_);
                cv_.wait(lock, [this]() { return stop_ || !queue_.empty(); });
                if (stop_ && queue_.empty()) {
                    return;
                }
                job = std::move(queue_.front());
                queue_.pop();
            }
            job();
        }
    }

    std::vector<std::thread> threads_;
    std::queue<std::function<void()>> queue_;
    std::mutex mutex_;
    std::condition_variable cv_;
    bool stop_;
};

class QuantumState {
public:
    explicit QuantumState(size_t qubits) : qubits_(qubits) {
        size_t size = 1ULL << qubits_;
        amplitudes_.assign(size, {0.0, 0.0});
        amplitudes_[0] = {1.0, 0.0};
    }

    void apply_h(size_t q) {
        apply_single_qubit_gate(q, {{1.0 / std::sqrt(2.0), 0.0}, {1.0 / std::sqrt(2.0), 0.0}},
                                   {{1.0 / std::sqrt(2.0), 0.0}, {-1.0 / std::sqrt(2.0), 0.0}});
    }

    void apply_x(size_t q) {
        apply_single_qubit_gate(q, {{0.0, 0.0}, {1.0, 0.0}}, {{1.0, 0.0}, {0.0, 0.0}});
    }

    void apply_ry(size_t q, double theta) {
        double c = std::cos(theta / 2.0);
        double s = std::sin(theta / 2.0);
        apply_single_qubit_gate(q, {{c, 0.0}, {-s, 0.0}}, {{s, 0.0}, {c, 0.0}});
    }

    void apply_rz(size_t q, double theta) {
        std::complex<double> e0(std::cos(-theta / 2.0), std::sin(-theta / 2.0));
        std::complex<double> e1(std::cos(theta / 2.0), std::sin(theta / 2.0));
        apply_single_qubit_gate(q, {{e0}, {0.0, 0.0}}, {{0.0, 0.0}, {e1}});
    }

    void apply_cnot(size_t control, size_t target) {
        size_t size = amplitudes_.size();
        for (size_t i = 0; i < size; ++i) {
            if (((i >> control) & 1U) == 1U && ((i >> target) & 1U) == 0U) {
                size_t flipped = i ^ (1ULL << target);
                std::swap(amplitudes_[i], amplitudes_[flipped]);
            }
        }
    }

    std::vector<double> probabilities() const {
        std::vector<double> probs;
        probs.reserve(amplitudes_.size());
        for (const auto &amp : amplitudes_) {
            probs.push_back(std::norm(amp));
        }
        return probs;
    }

private:
    void apply_single_qubit_gate(
        size_t q,
        std::array<std::complex<double>, 2> row0,
        std::array<std::complex<double>, 2> row1
    ) {
        size_t size = amplitudes_.size();
        size_t step = 1ULL << q;
        for (size_t i = 0; i < size; i += 2 * step) {
            for (size_t j = 0; j < step; ++j) {
                size_t idx0 = i + j;
                size_t idx1 = idx0 + step;
                auto a0 = amplitudes_[idx0];
                auto a1 = amplitudes_[idx1];
                amplitudes_[idx0] = row0[0] * a0 + row0[1] * a1;
                amplitudes_[idx1] = row1[0] * a0 + row1[1] * a1;
            }
        }
    }

    size_t qubits_;
    std::vector<std::complex<double>> amplitudes_;
};

std::string to_string(const BigInt &value) {
    return value.convert_to<std::string>();
}

std::string to_string(const Rational &value) {
    std::string num = value.numerator().convert_to<std::string>();
    std::string den = value.denominator().convert_to<std::string>();
    if (den == "1") {
        return num;
    }
    return num + "/" + den;
}

double rational_to_angle(const Rational &value) {
    BigInt num = value.numerator();
    BigInt den = value.denominator();
    if (den == 0) {
        return 0.0;
    }
    BigInt bounded_mod = BigInt(360) * den;
    BigInt bounded = num % bounded_mod;
    if (bounded < 0) {
        bounded += bounded_mod;
    }
    long double numerator = bounded.convert_to<long double>();
    long double denominator = den.convert_to<long double>();
    long double degrees = numerator / denominator;
    long double radians = degrees * (kPi / 180.0L);
    long double angle = std::fmod(radians, kTwoPi);
    if (angle < 0) {
        angle += kTwoPi;
    }
    return static_cast<double>(angle);
}

struct SutraTask {
    std::string name;
    std::function<SutraOutput()> run;
};

std::vector<SutraTask> build_tasks() {
    using namespace vedic;
    std::vector<SutraTask> tasks;

    tasks.push_back({"S1_Ekadhikena", []() {
        auto result = S1_Ekadhikena::divide_by_nine_ender(BigInt(19), 200);
        int sum = 0;
        for (int digit : result.recurring) {
            sum += digit;
        }
        Rational value(sum);
        return SutraOutput{"S1_Ekadhikena", value, "Sum of recurring digits for 1/19", rational_to_angle(value)};
    }});

    tasks.push_back({"S2_Nikhilam", []() {
        auto result = S2_Nikhilam::multiply(BigInt(98), BigInt(97), BigInt(100));
        Rational value(result.product);
        return SutraOutput{"S2_Nikhilam", value, "Nikhilam product 98×97", rational_to_angle(value)};
    }});

    tasks.push_back({"S3_Urdhva", []() {
        auto result = S3_Urdhva::multiply(BigInt(123), BigInt(456));
        Rational value(result.product);
        return SutraOutput{"S3_Urdhva", value, "Urdhva multiplication 123×456", rational_to_angle(value)};
    }});

    tasks.push_back({"S4_Paravartya", []() {
        auto result = S4_Paravartya::divide(BigInt(9506), BigInt(98));
        Rational value(result.quotient);
        return SutraOutput{"S4_Paravartya", value, "Paravartya division 9506÷98", rational_to_angle(value)};
    }});

    tasks.push_back({"S5_Shunyam", []() {
        auto result = S5_Shunyam::solve_product_equality(Rational(1), Rational(2), Rational(3), Rational(1));
        Rational value = result.solutions.empty() ? Rational(0) : result.solutions.front();
        return SutraOutput{"S5_Shunyam", value, "Solve (x+1)(x+2)=(x+3)(x+1)", rational_to_angle(value)};
    }});

    tasks.push_back({"S6_Anurupye", []() {
        auto result = S6_Anurupye::solve_system_2x2(
            Rational(2), Rational(3), Rational(8),
            Rational(1), Rational(2), Rational(5)
        );
        Rational value = result.x.value_or(Rational(0)) + result.y.value_or(Rational(0));
        return SutraOutput{"S6_Anurupye", value, "Solve 2x+3y=8 and x+2y=5", rational_to_angle(value)};
    }});

    tasks.push_back({"S7_Sankalana", []() {
        auto result = S7_Sankalana::solve_by_elimination(
            Rational(1), Rational(1), Rational(5),
            Rational(1), Rational(-1), Rational(1)
        );
        Rational value = result.x.value_or(Rational(0)) * result.y.value_or(Rational(0));
        return SutraOutput{"S7_Sankalana", value, "Elimination result x*y for x+y=5, x-y=1", rational_to_angle(value)};
    }});

    tasks.push_back({"S8_Purana", []() {
        auto result = S8_Purana::complete_the_square(Rational(2), Rational(4), Rational(1));
        Rational value = result.h + result.k;
        return SutraOutput{"S8_Purana", value, "Complete square for 2x^2+4x+1", rational_to_angle(value)};
    }});

    tasks.push_back({"S9_Calana", []() {
        std::vector<Rational> poly = {Rational(1), Rational(-3), Rational(2)};
        auto deriv = S9_Calana::differentiate(poly);
        Rational value = S9_Calana::evaluate(deriv, Rational(2));
        return SutraOutput{"S9_Calana", value, "Derivative of 1-3x+2x^2 at x=2", rational_to_angle(value)};
    }});

    tasks.push_back({"S10_Yavadunam", []() {
        auto result = S10_Yavadunam::square(BigInt(997), BigInt(1000));
        Rational value(result.square);
        return SutraOutput{"S10_Yavadunam", value, "Square 997 using base 1000", rational_to_angle(value)};
    }});

    tasks.push_back({"S11_Vyashti", []() {
        std::vector<Rational> masses = {Rational(1), Rational(3), Rational(5)};
        Rational value = S11_Vyashti::total_energy(masses);
        return SutraOutput{"S11_Vyashti", value, "Total energy for masses 1,3,5", rational_to_angle(value)};
    }});

    tasks.push_back({"S12_Sesanyankena", []() {
        BigInt remainder = S12_Sesanyankena::mod_pow(BigInt(7), BigInt(128), BigInt(13));
        Rational value(remainder);
        return SutraOutput{"S12_Sesanyankena", value, "7^128 mod 13", rational_to_angle(value)};
    }});

    tasks.push_back({"S13_Sopantya", []() {
        auto convergent = S13_Sopantya::golden_ratio_convergent(12);
        Rational value(convergent.numerator, convergent.denominator);
        return SutraOutput{"S13_Sopantya", value, "Golden ratio convergent n=12", rational_to_angle(value)};
    }});

    tasks.push_back({"S14_Ekanyunena", []() {
        auto result = S14_Ekanyunena::multiply_by_nines(BigInt(12345), 4);
        Rational value(result.product);
        return SutraOutput{"S14_Ekanyunena", value, "12345×9999", rational_to_angle(value)};
    }});

    tasks.push_back({"S15_Gunitasamuccaya", []() {
        auto result = S15_Gunitasamuccaya::verify_distributive(
            Rational(2), Rational(3), Rational(4), Rational(5)
        );
        Rational value = result.left_side;
        return SutraOutput{"S15_Gunitasamuccaya", value, "(2+3)(4+5)", rational_to_angle(value)};
    }});

    tasks.push_back({"S16_Gunakasamuccaya", []() {
        BigInt a = 84;
        BigInt b = 126;
        BigInt gcd = vedic::util::gcd(a, b);
        Rational value(gcd);
        return SutraOutput{"S16_Gunakasamuccaya", value, "GCD(84,126)", rational_to_angle(value)};
    }});

    tasks.push_back({"US1_Anurupyena", []() {
        auto parts = US1_Anurupyena::divide_proportionally(Rational(250), Rational(3), Rational(7));
        Rational value = parts.first;
        return SutraOutput{"US1_Anurupyena", value, "250 split 3:7", rational_to_angle(value)};
    }});

    tasks.push_back({"US2_Shishyate", []() {
        auto cycle = US2_Shishyate::detect_linear_cycle(BigInt(5), BigInt(1), BigInt(97), BigInt(3));
        Rational value(static_cast<long long>(cycle.cycle_length));
        return SutraOutput{"US2_Shishyate", value, "Cycle length for (5x+1) mod 97", rational_to_angle(value)};
    }});

    tasks.push_back({"US3_Adyam", []() {
        std::vector<Rational> values = {Rational(1), Rational(3), Rational(6), Rational(10)};
        auto bounds = US3_Adyam::check_bounds(values);
        Rational value = bounds.max_bound;
        return SutraOutput{"US3_Adyam", value, "Max bound from endpoint scan", rational_to_angle(value)};
    }});

    tasks.push_back({"US4_Kevalaih", []() {
        BigInt product = US4_Kevalaih::multiply_by_7(BigInt(12345));
        Rational value(product);
        return SutraOutput{"US4_Kevalaih", value, "12345×7 via complements", rational_to_angle(value)};
    }});

    tasks.push_back({"US5_Vestanam", []() {
        auto osc = US5_Vestanam::find_negative_osculator(BigInt(7));
        Rational value = osc ? Rational(*osc) : Rational(0);
        return SutraOutput{"US5_Vestanam", value, "Negative osculator for 7", rational_to_angle(value)};
    }});

    tasks.push_back({"US6_Yavadunam_Squared", []() {
        auto result = US6_Yavadunam_Squared::square_by_deficiency(BigInt(97), BigInt(100));
        Rational value(result.result);
        return SutraOutput{"US6_Yavadunam_Squared", value, "97^2 by deficiency", rational_to_angle(value)};
    }});

    tasks.push_back({"US7_Yavadunam_Extended", []() {
        auto result = US7_Yavadunam_Extended::square_extended(BigInt(103), BigInt(100));
        Rational value(result.result);
        return SutraOutput{"US7_Yavadunam_Extended", value, "103^2 extended", rational_to_angle(value)};
    }});

    tasks.push_back({"US8_Antyayor", []() {
        auto result = US8_Antyayor::multiply_sum_to_ten(BigInt(43), BigInt(47));
        Rational value(result.product);
        return SutraOutput{"US8_Antyayor", value, "43×47 (sum to ten)", rational_to_angle(value)};
    }});

    tasks.push_back({"US9_Antyayoreva", []() {
        int last_digit = US9_Antyayoreva::last_digit_of_power(BigInt(7), BigInt(100));
        Rational value(last_digit);
        return SutraOutput{"US9_Antyayoreva", value, "Last digit of 7^100", rational_to_angle(value)};
    }});

    tasks.push_back({"US10_Samuccayagunitah", []() {
        std::vector<Rational> a = {Rational(1), Rational(2), Rational(3)};
        std::vector<Rational> b = {Rational(4), Rational(5), Rational(6)};
        Rational value = US10_Samuccayagunitah::dot_product(a, b);
        return SutraOutput{"US10_Samuccayagunitah", value, "Dot product [1,2,3]·[4,5,6]", rational_to_angle(value)};
    }});

    tasks.push_back({"US11_Lopanasthapana", []() {
        std::vector<std::vector<Rational>> A = {
            {Rational(1), Rational(1), Rational(1)},
            {Rational(2), Rational(-1), Rational(1)},
            {Rational(3), Rational(2), Rational(-1)}
        };
        std::vector<Rational> b = {Rational(6), Rational(3), Rational(4)};
        auto step = US11_Lopanasthapana::eliminate_variable(A, b, 0);
        Rational value = step.reduced_matrix.empty() ? Rational(0) : step.reduced_matrix[0][0];
        return SutraOutput{"US11_Lopanasthapana", value, "Eliminate x from 3x3 system", rational_to_angle(value)};
    }});

    tasks.push_back({"US12_Vilokanam", []() {
        auto result = US12_Vilokanam::check_difference_of_squares(BigInt(45));
        Rational value = result.parameters.empty() ? Rational(0) : result.parameters[0];
        return SutraOutput{"US12_Vilokanam", value, "Difference of squares check for 45", rational_to_angle(value)};
    }});

    tasks.push_back({"US13_Gunitasamuccaya_Samuccayagunitah", []() {
        std::vector<Rational> a = {Rational(1), Rational(2), Rational(3)};
        std::vector<Rational> b = {Rational(4), Rational(1), Rational(2)};
        auto result = US13_Gunitasamuccaya_Samuccayagunitah::verify_consistency(a, b);
        Rational value = result.product_of_sums;
        return SutraOutput{"US13_Gunitasamuccaya_Samuccayagunitah", value, "Consistency of sums/products", rational_to_angle(value)};
    }});

    return tasks;
}

Rational aggregate_classical(const std::vector<SutraOutput> &outputs) {
    Rational total(0);
    for (const auto &out : outputs) {
        total += vedic::abs_rational(out.value);
    }
    return total;
}

SimulationModeResult run_serial(const std::vector<SutraTask> &tasks, size_t qubits) {
    SimulationModeResult result;
    result.mode = "serial";
    for (const auto &task : tasks) {
        result.outputs.push_back(task.run());
    }
    result.classical_total = aggregate_classical(result.outputs);

    QuantumState state(qubits);
    for (size_t q = 0; q < qubits; ++q) {
        state.apply_h(q);
    }

    for (size_t i = 0; i < result.outputs.size(); ++i) {
        size_t q = i % qubits;
        double angle = result.outputs[i].angle;
        if (i % 2 == 0) {
            state.apply_ry(q, angle);
        } else {
            state.apply_rz(q, angle);
        }
        if (i % 5 == 0 && qubits > 1) {
            size_t target = (q + 1) % qubits;
            state.apply_cnot(q, target);
        }
        if (i % 7 == 0) {
            state.apply_x(q);
        }
    }

    result.probabilities = state.probabilities();
    result.quantum_expectation = 0.0;
    for (size_t i = 0; i < result.probabilities.size(); ++i) {
        result.quantum_expectation += static_cast<double>(i) * result.probabilities[i];
    }
    return result;
}

SimulationModeResult run_concurrent(const std::vector<SutraTask> &tasks, size_t qubits) {
    SimulationModeResult result;
    result.mode = "concurrent";
    std::vector<std::future<SutraOutput>> futures;
    futures.reserve(tasks.size());

    for (const auto &task : tasks) {
        futures.push_back(std::async(std::launch::async, task.run));
    }
    for (auto &future : futures) {
        result.outputs.push_back(future.get());
    }

    result.classical_total = aggregate_classical(result.outputs);

    QuantumState state(qubits);
    for (size_t q = 0; q < qubits; ++q) {
        state.apply_h(q);
    }

    for (size_t i = 0; i < result.outputs.size(); ++i) {
        size_t q = (result.outputs.size() - 1 - i) % qubits;
        double angle = result.outputs[i].angle;
        if (i % 3 == 0) {
            state.apply_rz(q, angle);
        } else {
            state.apply_ry(q, angle);
        }
        if (i % 4 == 0 && qubits > 1) {
            size_t target = (q + qubits - 1) % qubits;
            state.apply_cnot(q, target);
        }
    }

    result.probabilities = state.probabilities();
    result.quantum_expectation = 0.0;
    for (size_t i = 0; i < result.probabilities.size(); ++i) {
        result.quantum_expectation += static_cast<double>(i) * result.probabilities[i];
    }
    return result;
}

SimulationModeResult run_parallel(const std::vector<SutraTask> &tasks, size_t qubits) {
    SimulationModeResult result;
    result.mode = "parallel";

    size_t workers = std::min(tasks.size(), static_cast<size_t>(std::max(1u, std::thread::hardware_concurrency())));
    ThreadPool pool(workers);
    std::vector<std::future<SutraOutput>> futures;
    futures.reserve(tasks.size());

    for (const auto &task : tasks) {
        futures.push_back(pool.submit(task.run));
    }
    for (auto &future : futures) {
        result.outputs.push_back(future.get());
    }

    result.classical_total = aggregate_classical(result.outputs);

    QuantumState state(qubits);
    for (size_t q = 0; q < qubits; ++q) {
        state.apply_h(q);
    }

    for (size_t i = 0; i < result.outputs.size(); ++i) {
        size_t q = (i * 2) % qubits;
        double angle = result.outputs[i].angle;
        state.apply_ry(q, angle);
        if (i % 6 == 0 && qubits > 1) {
            size_t target = (q + 2) % qubits;
            if (target != q) {
                state.apply_cnot(q, target);
            }
        }
    }

    result.probabilities = state.probabilities();
    result.quantum_expectation = 0.0;
    for (size_t i = 0; i < result.probabilities.size(); ++i) {
        result.quantum_expectation += static_cast<double>(i) * result.probabilities[i];
    }
    return result;
}

void write_report(
    const std::string &path,
    const std::vector<SimulationModeResult> &modes,
    size_t qubits,
    size_t shots
) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Unable to open report file: " + path);
    }

    out << "{\n";
    out << "  \"qubits\": " << qubits << ",\n";
    out << "  \"shots\": " << shots << ",\n";
    out << "  \"modes\": [\n";
    for (size_t i = 0; i < modes.size(); ++i) {
        const auto &mode = modes[i];
        out << "    {\n";
        out << "      \"mode\": \"" << mode.mode << "\",\n";
        out << "      \"classical_total\": \"" << to_string(mode.classical_total) << "\",\n";
        out << "      \"quantum_expectation\": " << std::setprecision(12) << mode.quantum_expectation << ",\n";
        out << "      \"sutras\": [\n";
        for (size_t j = 0; j < mode.outputs.size(); ++j) {
            const auto &sutra = mode.outputs[j];
            out << "        {\n";
            out << "          \"name\": \"" << sutra.name << "\",\n";
            out << "          \"value\": \"" << to_string(sutra.value) << "\",\n";
            out << "          \"detail\": \"" << sutra.detail << "\",\n";
            out << "          \"angle\": " << std::setprecision(12) << sutra.angle << "\n";
            out << "        }";
            if (j + 1 < mode.outputs.size()) {
                out << ",";
            }
            out << "\n";
        }
        out << "      ],\n";
        out << "      \"probabilities\": [";
        for (size_t j = 0; j < mode.probabilities.size(); ++j) {
            out << std::setprecision(12) << mode.probabilities[j];
            if (j + 1 < mode.probabilities.size()) {
                out << ", ";
            }
        }
        out << "]\n";
        out << "    }";
        if (i + 1 < modes.size()) {
            out << ",";
        }
        out << "\n";
    }
    out << "  ]\n";
    out << "}\n";
}

void print_usage(const char *argv0) {
    std::cout << "Usage: " << argv0 << " --mode <serial|concurrent|parallel|all> --qubits <n> --shots <n> --report <path>\n";
}

} // namespace simulator

int main(int argc, char **argv) {
    using namespace simulator;

    std::string mode = "all";
    size_t qubits = 4;
    size_t shots = 2048;
    std::string report_path = "simulation_report.json";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--mode" && i + 1 < argc) {
            mode = argv[++i];
        } else if (arg == "--qubits" && i + 1 < argc) {
            qubits = static_cast<size_t>(std::stoul(argv[++i]));
        } else if (arg == "--shots" && i + 1 < argc) {
            shots = static_cast<size_t>(std::stoul(argv[++i]));
        } else if (arg == "--report" && i + 1 < argc) {
            report_path = argv[++i];
        } else if (arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (qubits == 0 || qubits > 8) {
        std::cerr << "Qubit count must be between 1 and 8 for this simulator.\n";
        return 1;
    }

    auto tasks = build_tasks();
    std::vector<SimulationModeResult> results;

    if (mode == "serial" || mode == "all") {
        results.push_back(run_serial(tasks, qubits));
    }
    if (mode == "concurrent" || mode == "all") {
        results.push_back(run_concurrent(tasks, qubits));
    }
    if (mode == "parallel" || mode == "all") {
        results.push_back(run_parallel(tasks, qubits));
    }
    if (results.empty()) {
        std::cerr << "Invalid mode: " << mode << "\n";
        return 1;
    }

    try {
        write_report(report_path, results, qubits, shots);
    } catch (const std::exception &ex) {
        std::cerr << "Failed to write report: " << ex.what() << "\n";
        return 1;
    }

    std::cout << "Hybrid simulation complete. Report written to " << report_path << "\n";
    return 0;
}
