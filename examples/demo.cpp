#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

#include "battery_soc267/estimator.hpp"

namespace {

using battery_soc267::MeasurementStatus;
using battery_soc267::ModelParams;
using battery_soc267::ProcessNoise;
using battery_soc267::SocEstimator;
using battery_soc267::StateVector;
using battery_soc267::StepResult;

const char* StatusName(MeasurementStatus s) {
  switch (s) {
    case MeasurementStatus::kNoMeasurement: return "no-meas ";
    case MeasurementStatus::kAccepted:      return "accepted";
    case MeasurementStatus::kRejected:      return "REJECTED";
  }
  return "?";
}

// Simulated plant: integrates the same 2-RC dynamics with a true sensor bias.
struct Plant {
  ModelParams p;
  double soc;
  double v1 = 0.0;
  double v2 = 0.0;
  double true_bias;  // sensor reads i_true + true_bias

  double Step(double dt, double i_true) {
    soc -= i_true * dt / (p.capacity_ah * 3600.0);
    const double a1 = std::exp(-dt / (p.r1 * p.c1));
    const double a2 = std::exp(-dt / (p.r2 * p.c2));
    v1 = v1 * a1 + i_true * p.r1 * (1.0 - a1);
    v2 = v2 * a2 + i_true * p.r2 * (1.0 - a2);
    return battery_soc267::Ocv(p, soc) - p.r0 * i_true - v1 - v2;
  }
};

void PrintStep(int k, double dt, double i_meas, bool has_v, double v_meas,
               const StepResult& r, const Plant& plant) {
  std::cout << std::fixed << std::setprecision(5)
            << "k=" << std::setw(3) << k
            << " dt=" << dt << "s i_meas=" << std::setw(7) << i_meas << "A";
  if (has_v)
    std::cout << " v_meas=" << v_meas << "V";
  else
    std::cout << " v_meas=  (none)";
  std::cout << " | " << StatusName(r.status)
            << " soc=" << r.state(0) << " (true " << plant.soc << ")"
            << " bias=" << r.state(3) << "A (true " << plant.true_bias << ")"
            << " v_pred=" << r.predicted_voltage;
  if (r.status != MeasurementStatus::kNoMeasurement)
    std::cout << " innov=" << r.innovation << " nis=" << r.nis;
  std::cout << "\n";
}

}  // namespace

int main() {
  ModelParams params;
  params.capacity_ah = 50.0;
  params.r0 = 0.01;
  params.r1 = 0.005;
  params.c1 = 2000.0;
  params.r2 = 0.004;
  params.c2 = 20000.0;
  params.ocv_soc = {0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0};
  params.ocv_voltage = {3.00, 3.20, 3.35, 3.50, 3.62, 3.75, 3.90, 4.15};

  // Estimator starts with a wrong bias guess (0 A) and slightly wrong SOC.
  StateVector x0;
  x0 << 0.75, 0.0, 0.0, 0.0;
  Eigen::Matrix<double, 4, 4> p0 = Eigen::Matrix<double, 4, 4>::Zero();
  p0.diagonal() << 1e-3, 1e-2, 1e-2, 1e-2;

  ProcessNoise q;
  q.soc = 1e-8;
  q.rc1 = 1e-7;
  q.rc2 = 1e-7;
  q.bias = 1e-5;  // lets the bias state converge to the sensor offset

  const double r_voltage = 1e-4;   // 10 mV std-dev voltage noise
  const double nis_gate = 9.0;     // reject measurements with NIS > 9

  SocEstimator est;
  std::string error;
  if (!est.Initialize(params, x0, p0, q, r_voltage, nis_gate, error)) {
    std::cerr << "init failed: " << error << "\n";
    return 1;
  }

  // Plant: true SOC 0.80, sensor offset +0.4 A (reads high), 5 mV noise.
  Plant plant{params, 0.80, 0.0, 0.0, 0.4};
  std::mt19937 rng(42);
  std::normal_distribution<double> voltage_noise(0.0, 0.005);

  // SOC and bias are only weakly jointly observable (a constant voltage
  // offset fits both), so joint convergence needs a longer horizon.
  std::cout << "=== Phase 1: +/-10 A cycles, bias +0.4 A, 30000 x 1 s ===\n";
  const double dt = 1.0;
  int k = 0;
  for (; k < 30000; ++k) {
    // Alternate discharge/charge every 60 s for observability excitation.
    const double i_true = (k / 60) % 2 == 0 ? 10.0 : -10.0;
    const double v = plant.Step(dt, i_true) + voltage_noise(rng);
    const double i_meas = i_true + plant.true_bias;
    StepResult r;
    if (!est.Step(dt, i_meas, true, v, r, error)) {
      std::cerr << "step failed: " << error << "\n";
      return 1;
    }
    if (k < 3 || (k + 1) % 5000 == 0)
      PrintStep(k + 1, dt, i_meas, true, v, r, plant);
  }

  std::cout << "\n=== Phase 2: missing voltage (prediction only) ===\n";
  for (int m = 0; m < 3; ++m, ++k) {
    const double i_true = 10.0;
    const double i_meas = i_true + plant.true_bias;
    plant.Step(dt, i_true);
    StepResult r;
    if (!est.Step(dt, i_meas, false, 0.0, r, error)) {
      std::cerr << "step failed: " << error << "\n";
      return 1;
    }
    PrintStep(k + 1, dt, i_meas, false, 0.0, r, plant);
  }

  std::cout << "\n=== Phase 3: outlier voltage rejected by NIS gate ===\n";
  {
    const double i_true = 10.0;
    const double v = plant.Step(dt, i_true) + voltage_noise(rng);
    const double bad_v = v + 0.5;  // 500 mV spike
    StepResult r;
    if (!est.Step(dt, i_true + plant.true_bias, true, bad_v, r, error)) {
      std::cerr << "step failed: " << error << "\n";
      return 1;
    }
    PrintStep(++k, dt, i_true + plant.true_bias, true, bad_v, r, plant);
  }

  std::cout << "\n=== Phase 4: 8 A charge (negative current) ===\n";
  const double i_chg = -8.0;
  for (int m = 0; m < 20; ++m, ++k) {
    const double v = plant.Step(dt, i_chg) + voltage_noise(rng);
    StepResult r;
    if (!est.Step(dt, i_chg + plant.true_bias, true, v, r, error)) {
      std::cerr << "step failed: " << error << "\n";
      return 1;
    }
    if ((m + 1) % 5 == 0)
      PrintStep(k + 1, dt, i_chg + plant.true_bias, true, v, r, plant);
  }

  std::cout << "\n=== Phase 5: invalid input is rejected, state kept ===\n";
  {
    const auto before = est.state();
    StepResult r;
    if (est.Step(-1.0, 0.0, false, 0.0, r, error))
      std::cout << "unexpected success\n";
    else
      std::cout << "negative dt rejected: " << error << "\n";
    std::cout << "state unchanged: "
              << (est.state().isApprox(before) ? "yes" : "NO") << "\n";
  }

  std::cout << "\nFinal estimate: soc=" << est.soc()
            << " (true " << plant.soc << ")"
            << " bias=" << est.bias() << "A (true " << plant.true_bias
            << "A)\n";
  return 0;
}
