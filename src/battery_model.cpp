#include "battery_soc267/battery_model.hpp"

#include <algorithm>
#include <cmath>

namespace battery_soc267 {
namespace {

bool IsFinite(double value) {
  return std::isfinite(value);
}

bool IsPositiveFinite(double value) {
  return std::isfinite(value) && value > 0.0;
}

}  // namespace

BatteryModel::BatteryModel(ModelConfig config) : config_(std::move(config)) {
  ValidateConfig();
}

void BatteryModel::ValidateConfig() const {
  if (!IsPositiveFinite(config_.capacity_ah)) {
    throw BatterySocError("capacity_ah must be positive and finite");
  }
  const double parameters[] = {config_.r0, config_.r1, config_.c1,
                               config_.r2, config_.c2};
  for (double value : parameters) {
    if (!IsPositiveFinite(value)) {
      throw BatterySocError("R0/R1/C1/R2/C2 must be positive and finite");
    }
  }
  if (config_.ocv_table.size() < 2) {
    throw BatterySocError("OCV table must contain at least two nodes");
  }
  for (const OcvPoint& point : config_.ocv_table) {
    if (!IsFinite(point.soc) || !IsFinite(point.voltage)) {
      throw BatterySocError("OCV table contains a non-finite value");
    }
  }
  if (config_.ocv_table.front().soc != 0.0 ||
      config_.ocv_table.back().soc != 1.0) {
    throw BatterySocError("OCV SOC nodes must cover exactly 0 to 1");
  }
  for (std::size_t i = 1; i < config_.ocv_table.size(); ++i) {
    if (!(config_.ocv_table[i].soc > config_.ocv_table[i - 1].soc)) {
      throw BatterySocError("OCV SOC nodes must be strictly increasing");
    }
    if (config_.ocv_table[i].voltage < config_.ocv_table[i - 1].voltage) {
      throw BatterySocError("OCV voltages must be non-decreasing");
    }
  }
}

std::size_t BatteryModel::SegmentIndex(double soc) const {
  if (soc < 0.0 || soc > 1.0) {
    throw BatterySocError("SOC is outside [0, 1]");
  }
  auto it = std::upper_bound(
      config_.ocv_table.begin(), config_.ocv_table.end(), soc,
      [](double value, const OcvPoint& point) { return value < point.soc; });
  if (it == config_.ocv_table.begin()) {
    return 0;
  }
  if (it == config_.ocv_table.end()) {
    return config_.ocv_table.size() - 2;
  }
  return static_cast<std::size_t>(it - config_.ocv_table.begin() - 1);
}

double BatteryModel::OpenCircuitVoltage(double soc) const {
  const std::size_t i = SegmentIndex(soc);
  const double soc0 = config_.ocv_table[i].soc;
  const double soc1 = config_.ocv_table[i + 1].soc;
  const double v0 = config_.ocv_table[i].voltage;
  const double v1 = config_.ocv_table[i + 1].voltage;
  const double fraction = (soc - soc0) / (soc1 - soc0);
  return v0 + fraction * (v1 - v0);
}

double BatteryModel::OpenCircuitVoltageDerivative(double soc) const {
  const std::size_t i = SegmentIndex(soc);
  const std::size_t derivative_index =
      soc >= config_.ocv_table.back().soc ? config_.ocv_table.size() - 2 : i;
  const double soc0 = config_.ocv_table[derivative_index].soc;
  const double soc1 = config_.ocv_table[derivative_index + 1].soc;
  const double v0 = config_.ocv_table[derivative_index].voltage;
  const double v1 = config_.ocv_table[derivative_index + 1].voltage;
  return (v1 - v0) / (soc1 - soc0);
}

StateVector BatteryModel::PredictState(const StateVector& state,
                                       double measured_current_a,
                                       double duration_s) const {
  if (!IsPositiveFinite(duration_s)) {
    throw BatterySocError("step duration must be positive and finite");
  }
  if (!IsFinite(measured_current_a) || !state.allFinite()) {
    throw BatterySocError("prediction input contains a non-finite value");
  }
  if (state[0] < 0.0 || state[0] > 1.0) {
    throw BatterySocError("SOC is outside [0, 1]");
  }

  const double true_current = measured_current_a - state[3];
  const double capacity_as = config_.capacity_ah * 3600.0;
  const double tau1 = config_.r1 * config_.c1;
  const double tau2 = config_.r2 * config_.c2;
  const double e1 = std::exp(-duration_s / tau1);
  const double e2 = std::exp(-duration_s / tau2);

  StateVector predicted = state;
  predicted[0] -= true_current * duration_s / capacity_as;
  predicted[1] = state[1] * e1 +
                 config_.r1 * true_current * (1.0 - e1);
  predicted[2] = state[2] * e2 +
                 config_.r2 * true_current * (1.0 - e2);
  predicted[3] = state[3];
  if (!predicted.allFinite() || predicted[0] < 0.0 || predicted[0] > 1.0) {
    throw BatterySocError("predicted SOC is outside [0, 1]");
  }
  return predicted;
}

StateMatrix BatteryModel::StateJacobian(double duration_s) const {
  const double tau1 = config_.r1 * config_.c1;
  const double tau2 = config_.r2 * config_.c2;
  const double e1 = std::exp(-duration_s / tau1);
  const double e2 = std::exp(-duration_s / tau2);
  StateMatrix jacobian = StateMatrix::Identity();
  jacobian(0, 3) = duration_s / (config_.capacity_ah * 3600.0);
  jacobian(1, 3) = -config_.r1 * (1.0 - e1);
  jacobian(2, 3) = -config_.r2 * (1.0 - e2);
  return jacobian;
}

double BatteryModel::TerminalVoltage(const StateVector& state,
                                     double measured_current_a) const {
  if (!state.allFinite() || !IsFinite(measured_current_a)) {
    throw BatterySocError("measurement model input contains a non-finite value");
  }
  if (state[0] < 0.0 || state[0] > 1.0) {
    throw BatterySocError("SOC is outside [0, 1]");
  }
  const double true_current = measured_current_a - state[3];
  return OpenCircuitVoltage(state[0]) - config_.r0 * true_current -
         state[1] - state[2];
}

Eigen::RowVector4d BatteryModel::MeasurementJacobian(
    const StateVector& state) const {
  Eigen::RowVector4d jacobian;
  jacobian << OpenCircuitVoltageDerivative(state[0]), -1.0, -1.0,
      config_.r0;
  return jacobian;
}

}  // namespace battery_soc267
