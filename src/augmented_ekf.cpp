#include "battery_soc267/augmented_ekf.hpp"

#include <cmath>

namespace battery_soc267 {
namespace {

bool IsFinite(double value) {
  return std::isfinite(value);
}

}  // namespace

AugmentedEkf::AugmentedEkf(BatteryModel model, FilterConfig config)
    : model_(std::move(model)), config_(std::move(config)) {
  for (int i = 0; i < config_.process_noise_density.size(); ++i) {
    const double density = config_.process_noise_density[i];
    if (!std::isfinite(density) || density < 0.0) {
      throw BatterySocError(
          "process noise density entries must be non-negative and finite");
    }
  }
  if (!IsFinite(config_.voltage_noise_variance) ||
      config_.voltage_noise_variance <= 0.0) {
    throw BatterySocError("voltage noise variance must be positive and finite");
  }
  if (!IsFinite(config_.nis_threshold) || config_.nis_threshold <= 0.0) {
    throw BatterySocError("NIS threshold must be positive and finite");
  }
}

Prediction AugmentedEkf::Predict(const StateVector& state,
                                 const StateMatrix& covariance,
                                 double measured_current_a,
                                 double duration_s) const {
  if (!covariance.allFinite()) {
    throw BatterySocError("covariance contains a non-finite value");
  }
  Prediction prediction;
  prediction.state =
      model_.PredictState(state, measured_current_a, duration_s);
  const StateMatrix transition =
      model_.StateJacobian(duration_s);
  const StateMatrix process_noise =
      config_.process_noise_density.asDiagonal() * duration_s;
  prediction.covariance =
      transition * covariance * transition.transpose() + process_noise;
  prediction.predicted_voltage_v =
      model_.TerminalVoltage(prediction.state, measured_current_a);
  if (!transition.allFinite() || !prediction.covariance.allFinite() ||
      !std::isfinite(prediction.predicted_voltage_v)) {
    throw BatterySocError("non-finite value occurred during prediction");
  }
  return prediction;
}

Correction AugmentedEkf::Correct(const Prediction& prediction,
                                 double measured_voltage_v) const {
  if (!std::isfinite(measured_voltage_v)) {
    throw BatterySocError("measured terminal voltage must be finite");
  }
  Correction correction;
  correction.measurement_jacobian =
      model_.MeasurementJacobian(prediction.state);
  correction.innovation_v =
      measured_voltage_v - prediction.predicted_voltage_v;
  const Eigen::Vector4d cross_covariance =
      prediction.covariance * correction.measurement_jacobian.transpose();
  correction.innovation_variance_v2 =
      correction.measurement_jacobian * cross_covariance +
      config_.voltage_noise_variance;
  if (!(correction.innovation_variance_v2 > 0.0) ||
      !std::isfinite(correction.innovation_variance_v2) ||
      !cross_covariance.allFinite()) {
    throw BatterySocError("innovation variance must be positive and finite");
  }

  const Eigen::Vector4d kalman_gain =
      cross_covariance / correction.innovation_variance_v2;
  correction.state = prediction.state +
                     kalman_gain * correction.innovation_v;
  if (!correction.state.allFinite() || correction.state[0] < 0.0 ||
      correction.state[0] > 1.0) {
    throw BatterySocError("corrected SOC is outside [0, 1]");
  }

  const StateMatrix identity_minus = StateMatrix::Identity() -
                                     kalman_gain *
                                         correction.measurement_jacobian;
  correction.covariance =
      identity_minus * prediction.covariance * identity_minus.transpose() +
      config_.voltage_noise_variance * kalman_gain * kalman_gain.transpose();
  correction.covariance =
      0.5 * (correction.covariance + correction.covariance.transpose());
  correction.nis = correction.innovation_v * correction.innovation_v /
                   correction.innovation_variance_v2;
  if (!correction.covariance.allFinite() || !std::isfinite(correction.nis)) {
    throw BatterySocError("non-finite value occurred during correction");
  }
  return correction;
}

}  // namespace battery_soc267
