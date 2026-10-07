// Self-tests for battery_soc267. Returns non-zero on the first failure.
#include <cmath>
#include <iostream>
#include <string>

#include "battery_soc267/estimator.hpp"

namespace {

using namespace battery_soc267;

int failures = 0;

void Check(bool cond, const std::string& name) {
  if (cond) {
    std::cout << "[PASS] " << name << "\n";
  } else {
    std::cout << "[FAIL] " << name << "\n";
    ++failures;
  }
}

ModelParams MakeParams() {
  ModelParams p;
  p.capacity_ah = 50.0;
  p.r0 = 0.01;
  p.r1 = 0.005;
  p.c1 = 2000.0;
  p.r2 = 0.004;
  p.c2 = 20000.0;
  p.ocv_soc = {0.0, 0.5, 1.0};
  p.ocv_voltage = {3.0, 3.6, 4.2};
  return p;
}

ProcessNoise MakeNoise() {
  ProcessNoise q;
  q.soc = 1e-8;
  q.rc1 = 1e-7;
  q.rc2 = 1e-7;
  q.bias = 1e-6;
  return q;
}

SocEstimator MakeEstimator(const StateVector& x0) {
  SocEstimator est;
  std::string err;
  bool ok = est.Initialize(MakeParams(), x0, 1e-2, MakeNoise(), 1e-4, 9.0,
                           err);
  Check(ok, "estimator initializes: " + err);
  return est;
}

}  // namespace

int main() {
  const ModelParams p = MakeParams();

  // --- Model validation ---
  {
    std::string err;
    ModelParams bad = p;
    bad.capacity_ah = -1.0;
    Check(!ValidateParams(bad, err), "negative capacity rejected");
    bad = p;
    bad.c1 = 0.0;
    Check(!ValidateParams(bad, err), "zero C1 rejected");
    bad = p;
    bad.ocv_soc = {0.0, 0.5, 0.5};
    Check(!ValidateParams(bad, err), "non-increasing SOC nodes rejected");
    bad = p;
    bad.ocv_voltage = {3.0, 3.6, 3.5};
    Check(!ValidateParams(bad, err), "decreasing OCV rejected");
    Check(ValidateParams(p, err), "valid params accepted");
  }

  // --- OCV interpolation and derivative convention ---
  {
    Check(std::fabs(Ocv(p, 0.25) - 3.3) < 1e-12, "OCV midpoint interp");
    Check(std::fabs(Ocv(p, 0.0) - 3.0) < 1e-12, "OCV first node");
    Check(std::fabs(Ocv(p, 1.0) - 4.2) < 1e-12, "OCV last node");
    // Right-segment slope at interior node 0.5: (4.2-3.6)/0.5 = 1.2
    Check(std::fabs(OcvDerivative(p, 0.5) - 1.2) < 1e-12,
          "node derivative uses right segment");
    // Last node uses left segment slope: 1.2 here as well; use 0.0 node:
    // right segment slope at node 0.0: (3.6-3.0)/0.5 = 1.2
    Check(std::fabs(OcvDerivative(p, 0.0) - 1.2) < 1e-12,
          "first node derivative uses right segment");
    ModelParams p2 = p;
    p2.ocv_voltage = {3.0, 3.6, 3.8};  // last segment slope 0.4
    Check(std::fabs(OcvDerivative(p2, 1.0) - 0.4) < 1e-12,
          "last node derivative uses left segment");
  }

  // --- Prediction: SOC coulomb counting, RC exponential, bias constant ---
  {
    StateVector x;
    x << 0.5, 0.1, -0.05, 0.2;  // bias 0.2 A
    const StateVector xn = PredictState(p, x, 1.2, 10.0);  // i_true = 1.0 A
    const double dsoc = 1.0 * 10.0 / (50.0 * 3600.0);
    Check(std::fabs(xn(0) - (0.5 - dsoc)) < 1e-15,
          "SOC uses true current (meas - bias)");
    const double a1 = std::exp(-10.0 / (p.r1 * p.c1));
    Check(std::fabs(xn(1) - (0.1 * a1 + 1.0 * p.r1 * (1.0 - a1))) < 1e-15,
          "RC1 exact exponential update");
    Check(xn(3) == 0.2, "bias prediction unchanged");
  }

  // --- Terminal voltage and Jacobian bias coupling ---
  {
    StateVector x;
    x << 0.5, 0.1, 0.05, 0.3;
    const double v = TerminalVoltage(p, x, 1.3);  // i_true = 1.0
    Check(std::fabs(v - (3.6 - 0.01 * 1.0 - 0.1 - 0.05)) < 1e-12,
          "terminal voltage equation");
    const auto h = TerminalVoltageJacobian(p, x);
    Check(h(0, 3) == p.r0, "Jacobian includes bias coupling +R0");
    Check(h(0, 1) == -1.0 && h(0, 2) == -1.0, "Jacobian RC terms");
  }

  // --- Filter: prediction-only step grows covariance, keeps mean dynamics ---
  {
    StateVector x0;
    x0 << 0.5, 0.0, 0.0, 0.0;
    SocEstimator est = MakeEstimator(x0);
    StepResult r;
    std::string err;
    Check(est.Step(1.0, 5.0, false, 0.0, r, err),
          "prediction-only step succeeds");
    Check(r.status == MeasurementStatus::kNoMeasurement,
          "status is no-measurement");
    Check(r.covariance(0, 0) > 1e-2, "process noise added to covariance");
    Check(std::fabs(r.state(0) - (0.5 - 5.0 / (50.0 * 3600.0))) < 1e-15,
          "predicted SOC matches coulomb counting");
  }

  // --- Filter: measurement update reduces covariance and reports NIS ---
  {
    StateVector x0;
    x0 << 0.5, 0.0, 0.0, 0.0;
    SocEstimator est = MakeEstimator(x0);
    StepResult r;
    std::string err;
    const double v_meas = TerminalVoltage(p, est.state(), 5.0) + 0.01;
    Check(est.Step(1.0, 5.0, true, v_meas, r, err),
          "update step succeeds");
    Check(r.status == MeasurementStatus::kAccepted, "measurement accepted");
    Check(std::fabs(r.innovation - (v_meas - r.predicted_voltage)) < 1e-12,
          "innovation value");
    Check(r.nis > 0.0, "NIS computed");
    Check(r.covariance(0, 0) < 1e-2, "covariance shrinks after update");
  }

  // --- Filter: NIS gate rejects outliers but keeps the prediction ---
  {
    StateVector x0;
    x0 << 0.5, 0.0, 0.0, 0.0;
    SocEstimator est = MakeEstimator(x0);
    StepResult r;
    std::string err;
    const double v_bad = TerminalVoltage(p, est.state(), 5.0) + 5.0;
    Check(est.Step(1.0, 5.0, true, v_bad, r, err), "outlier step succeeds");
    Check(r.status == MeasurementStatus::kRejected, "outlier rejected");
    Check(std::fabs(r.state(0) - est.state()(0)) < 1e-15 &&
              std::fabs(r.state(0) - (0.5 - 5.0 / (50.0 * 3600.0))) < 1e-15,
          "rejected step keeps pure prediction");
  }

  // --- Bias convergence: sensor offset is identified, not SOC distortion ---
  {
    StateVector x0;
    x0 << 0.6, 0.0, 0.0, 0.0;
    SocEstimator est;
    {
      ProcessNoise q = MakeNoise();
      q.bias = 1e-4;
      std::string init_err;
      Check(est.Initialize(p, x0, 1e-2, q, 1e-6, 9.0, init_err),
            "bias estimator initializes");
    }
    const double true_bias = 0.5;
    double soc_true = 0.6, v1 = 0.0, v2 = 0.0;
    std::string err;
    bool steps_ok = true;
    for (int k = 0; k < 20000; ++k) {
      const double i_true = (k % 100 < 50) ? 5.0 : -5.0;  // excitation
      soc_true -= i_true * 1.0 / (50.0 * 3600.0);
      const double a1 = std::exp(-1.0 / (p.r1 * p.c1));
      const double a2 = std::exp(-1.0 / (p.r2 * p.c2));
      v1 = v1 * a1 + i_true * p.r1 * (1.0 - a1);
      v2 = v2 * a2 + i_true * p.r2 * (1.0 - a2);
      const double v = Ocv(p, soc_true) - p.r0 * i_true - v1 - v2;
      StepResult r;
      steps_ok = est.Step(1.0, i_true + true_bias, true, v, r, err) &&
                 steps_ok;
    }
    Check(steps_ok, "bias run steps ok");
    Check(std::fabs(est.bias() - true_bias) < 0.05,
          "bias estimate converges to sensor offset");
    Check(std::fabs(est.soc() - soc_true) < 0.01,
          "SOC estimate tracks true SOC");
  }

  // --- Error handling: invalid input leaves state untouched ---
  {
    StateVector x0;
    x0 << 0.5, 0.0, 0.0, 0.0;
    SocEstimator est = MakeEstimator(x0);
    const auto before = est.state();
    StepResult r;
    std::string err;
    Check(!est.Step(0.0, 5.0, false, 0.0, r, err), "zero dt rejected");
    Check(!est.Step(-1.0, 5.0, false, 0.0, r, err), "negative dt rejected");
    Check(!est.Step(1.0, std::nan(""), false, 0.0, r, err),
          "NaN current rejected");
    Check(!est.Step(1.0, 5.0, true, std::nan(""), r, err),
          "NaN voltage rejected");
    Check(est.state().isApprox(before), "state unchanged after errors");
  }

  // --- Error handling: SOC leaving [0,1] fails without truncation ---
  {
    StateVector x0;
    x0 << 0.01, 0.0, 0.0, 0.0;
    SocEstimator est = MakeEstimator(x0);
    const auto before = est.state();
    StepResult r;
    std::string err;
    // 50 Ah cell, 500 A for 10 s -> way past empty.
    Check(!est.Step(10.0, 500.0, false, 0.0, r, err),
          "SOC out-of-range prediction fails");
    Check(est.state().isApprox(before), "state unchanged after SOC violation");
    Check(est.soc() == 0.01, "no truncation/clamping applied");
  }

  // --- Init validation ---
  {
    std::string err;
    SocEstimator est;
    StateVector x0;
    x0 << 1.5, 0.0, 0.0, 0.0;  // SOC outside table
    Check(!est.Initialize(p, x0, 1e-2, MakeNoise(), 1e-4, 9.0, err),
          "init SOC outside table rejected");
    x0 << 0.5, 0.0, 0.0, 0.0;
    Check(!est.Initialize(p, x0, 1e-2, MakeNoise(), -1.0, 9.0, err),
          "non-positive voltage noise rejected");
    Check(!est.Initialize(p, x0, 1e-2, MakeNoise(), 1e-4, 0.0, err),
          "non-positive NIS threshold rejected");
    ProcessNoise q = MakeNoise();
    q.bias = -1.0;
    Check(!est.Initialize(p, x0, 1e-2, q, 1e-4, 9.0, err),
          "negative process noise rejected");
    StateCovariance pbad = StateCovariance::Identity();
    pbad(0, 1) = 0.5;  // asymmetric
    Check(!est.Initialize(p, x0, pbad, MakeNoise(), 1e-4, 9.0, err),
          "non-SPD covariance rejected");
  }

  std::cout << (failures == 0 ? "ALL TESTS PASSED\n" : "FAILURES PRESENT\n");
  return failures == 0 ? 0 : 1;
}
