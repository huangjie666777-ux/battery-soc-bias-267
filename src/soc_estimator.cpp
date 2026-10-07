#include "battery_soc267/soc_estimator.hpp"

#include <cmath>

namespace battery_soc267 {

SocEstimator::SocEstimator(StateVector initial_state,
                           StateMatrix initial_covariance,
                           ModelConfig model_config,
                           FilterConfig filter_config)
    : ekf_(BatteryModel(std::move(model_config)), std::move(filter_config)),
      state_(std::move(initial_state)),
      covariance_(std::move(initial_covariance)) {
  ValidateInitialState(state_);
  ValidateCovariance(covariance_);
}

void SocEstimator::ValidateInitialState(const StateVector& state) {
  if (!state.allFinite()) {
    throw BatterySocError("initial state contains a non-finite value");
  }
  if (state[0] < 0.0 || state[0] > 1.0) {
    throw BatterySocError("initial SOC must be in [0, 1]");
  }
}

void SocEstimator::ValidateCovariance(const StateMatrix& covariance) {
  if (!covariance.allFinite()) {
    throw BatterySocError("initial covariance contains a non-finite value");
  }
  const StateMatrix symmetric = 0.5 * (covariance + covariance.transpose());
  if (!covariance.isApprox(symmetric, 1.0e-10)) {
    throw BatterySocError("initial covariance must be symmetric");
  }
  Eigen::SelfAdjointEigenSolver<StateMatrix> solver(symmetric);
  if (solver.info() != Eigen::Success ||
      solver.eigenvalues().minCoeff() <= 0.0) {
    throw BatterySocError("initial covariance must be positive definite");
  }
}

void SocEstimator::ValidateStepInput(const StepInput& input) {
  if (!std::isfinite(input.duration_s) || input.duration_s <= 0.0) {
    throw BatterySocError("step duration must be positive and finite");
  }
  if (!std::isfinite(input.current_a)) {
    throw BatterySocError("measured current must be finite");
  }
  if (input.terminal_voltage_v.has_value() &&
      !std::isfinite(*input.terminal_voltage_v)) {
    throw BatterySocError("measured terminal voltage must be finite");
  }
}

void SocEstimator::CheckFiniteState(const StateVector& state,
                                    const StateMatrix& covariance,
                                    const char* phase) {
  if (!state.allFinite() || !covariance.allFinite() || state[0] < 0.0 ||
      state[0] > 1.0) {
    throw BatterySocError(std::string("invalid state after ") + phase);
  }
}

StepResult SocEstimator::Step(const StepInput& input) {
  ValidateStepInput(input);
  const StateVector previous_state = state_;
  const StateMatrix previous_covariance = covariance_;

  Prediction prediction;
  try {
    prediction = ekf_.Predict(state_, covariance_, input.current_a,
                              input.duration_s);
  } catch (...) {
    state_ = previous_state;
    covariance_ = previous_covariance;
    throw;
  }

  StepResult result;
  result.predicted_voltage_v = prediction.predicted_voltage_v;
  if (!input.terminal_voltage_v.has_value()) {
    result.state = prediction.state;
    result.covariance = prediction.covariance;
    result.measurement_status = MeasurementStatus::kPredictedOnly;
    CheckFiniteState(result.state, result.covariance, "prediction");
    state_ = result.state;
    covariance_ = result.covariance;
    return result;
  }

  Correction correction;
  try {
    correction = ekf_.Correct(prediction, *input.terminal_voltage_v);
  } catch (...) {
    state_ = previous_state;
    covariance_ = previous_covariance;
    throw;
  }

  if (correction.nis > ekf_.config().nis_threshold) {
    result.state = prediction.state;
    result.covariance = prediction.covariance;
    result.measurement_status = MeasurementStatus::kRejected;
  } else {
    result.state = correction.state;
    result.covariance = correction.covariance;
    result.measurement_status = MeasurementStatus::kUpdated;
  }
  result.innovation_v = correction.innovation_v;
  result.nis = correction.nis;
  CheckFiniteState(result.state, result.covariance, "measurement update");
  state_ = result.state;
  covariance_ = result.covariance;
  return result;
}

}  // namespace battery_soc267
