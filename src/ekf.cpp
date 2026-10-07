#include "battery_soc267/ekf.hpp"

#include <cmath>

namespace battery_soc267 {
namespace {

bool IsFinite(double v) { return std::isfinite(v); }

bool IsSpd(const StateCovariance& p) {
  if (!p.allFinite()) return false;
  if (!p.isApprox(p.transpose(), 1e-10)) return false;
  Eigen::LLT<StateCovariance> llt(p);
  return llt.info() == Eigen::Success;
}

bool SocInRange(const ModelParams& params, double soc) {
  return soc >= params.ocv_soc.front() && soc <= params.ocv_soc.back();
}

}  // namespace

bool AugmentedEkf::Initialize(const ModelParams& params,
                              const StateVector& initial_state,
                              const StateCovariance& initial_covariance,
                              const ProcessNoise& process_noise,
                              double voltage_noise_variance,
                              double nis_threshold, std::string& error) {
  if (!ValidateParams(params, error)) return false;
  if (!initial_state.allFinite()) {
    error = "initial state must be finite";
    return false;
  }
  if (!SocInRange(params, initial_state(kIndexSoc))) {
    error = "initial SOC outside OCV table range";
    return false;
  }
  if (!IsSpd(initial_covariance)) {
    error = "initial covariance must be symmetric positive definite";
    return false;
  }
  const auto q = process_noise.Diagonal();
  for (int i = 0; i < 4; ++i) {
    if (!IsFinite(q(i)) || q(i) < 0.0) {
      error = "process noise densities must be non-negative and finite";
      return false;
    }
  }
  if (!IsFinite(voltage_noise_variance) || voltage_noise_variance <= 0.0) {
    error = "voltage noise variance must be positive and finite";
    return false;
  }
  if (!IsFinite(nis_threshold) || nis_threshold <= 0.0) {
    error = "NIS threshold must be positive and finite";
    return false;
  }

  params_ = params;
  state_ = initial_state;
  covariance_ = initial_covariance;
  process_noise_ = process_noise;
  voltage_noise_variance_ = voltage_noise_variance;
  nis_threshold_ = nis_threshold;
  return true;
}

bool AugmentedEkf::Step(double dt, double current_meas,
                        const double* terminal_voltage_meas,
                        StepResult& result, std::string& error) {
  if (!IsFinite(dt) || dt <= 0.0) {
    error = "dt must be positive and finite";
    return false;
  }
  if (!IsFinite(current_meas)) {
    error = "current measurement must be finite";
    return false;
  }
  if (terminal_voltage_meas != nullptr &&
      !IsFinite(*terminal_voltage_meas)) {
    error = "voltage measurement must be finite";
    return false;
  }

  // Work on copies so any failure leaves the filter untouched.
  StateVector state = PredictState(params_, state_, current_meas, dt);
  if (!state.allFinite()) {
    error = "numerical failure during state prediction";
    return false;
  }
  if (!SocInRange(params_, state(kIndexSoc))) {
    error = "predicted SOC out of OCV table range";
    return false;
  }

  StateCovariance cov = covariance_;
  cov.diagonal() += process_noise_.Diagonal() * dt;
  if (!cov.allFinite()) {
    error = "numerical failure during covariance prediction";
    return false;
  }

  result.state = state;
  result.covariance = cov;
  result.predicted_voltage =
      TerminalVoltage(params_, state, current_meas);
  result.status = MeasurementStatus::kNoMeasurement;
  result.innovation = 0.0;
  result.nis = 0.0;

  if (terminal_voltage_meas != nullptr) {
    const MeasurementJacobian h = TerminalVoltageJacobian(params_, state);
    const double s =
        (h * cov * h.transpose())(0, 0) + voltage_noise_variance_;
    if (!IsFinite(s) || s <= 0.0) {
      error = "numerical failure: non-positive innovation variance";
      return false;
    }
    const double innovation =
        *terminal_voltage_meas - result.predicted_voltage;
    const double nis = innovation * innovation / s;
    result.innovation = innovation;
    result.nis = nis;

    if (nis > nis_threshold_) {
      result.status = MeasurementStatus::kRejected;
    } else {
      const Eigen::Matrix<double, 4, 1> k =
          cov * h.transpose() / s;
      StateVector updated = state + k * innovation;
      if (!updated.allFinite()) {
        error = "numerical failure during state update";
        return false;
      }
      if (!SocInRange(params_, updated(kIndexSoc))) {
        error = "corrected SOC out of OCV table range";
        return false;
      }
      // Joseph form: P = (I-KH)P(I-KH)' + K R K'
      StateCovariance ikh = StateCovariance::Identity() - k * h;
      StateCovariance updated_cov =
          ikh * cov * ikh.transpose() +
          voltage_noise_variance_ * (k * k.transpose());
      if (!updated_cov.allFinite()) {
        error = "numerical failure during covariance update";
        return false;
      }
      // Enforce exact symmetry lost through floating-point roundoff.
      updated_cov =
          0.5 * (updated_cov + updated_cov.transpose()).eval();
      state = updated;
      cov = updated_cov;
      result.state = state;
      result.covariance = cov;
      result.status = MeasurementStatus::kAccepted;
    }
  }

  state_ = state;
  covariance_ = cov;
  return true;
}

}  // namespace battery_soc267

