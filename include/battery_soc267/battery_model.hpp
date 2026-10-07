#pragma once

#include "battery_soc267/types.hpp"

namespace battery_soc267 {

class BatteryModel {
 public:
  explicit BatteryModel(ModelConfig config);

  StateVector PredictState(const StateVector& state,
                           double measured_current_a,
                           double duration_s) const;
  StateMatrix StateJacobian(double duration_s) const;
  double TerminalVoltage(const StateVector& state,
                         double measured_current_a) const;
  Eigen::RowVector4d MeasurementJacobian(const StateVector& state) const;

  double OpenCircuitVoltage(double soc) const;
  double OpenCircuitVoltageDerivative(double soc) const;
  const ModelConfig& config() const { return config_; }

 private:
  void ValidateConfig() const;
  std::size_t SegmentIndex(double soc) const;

  ModelConfig config_;
};

}  // namespace battery_soc267
