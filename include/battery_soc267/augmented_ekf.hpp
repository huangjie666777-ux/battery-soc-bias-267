#pragma once

#include "battery_soc267/battery_model.hpp"

namespace battery_soc267 {

struct Prediction {
  StateVector state = StateVector::Zero();
  StateMatrix covariance = StateMatrix::Zero();
  double predicted_voltage_v = 0.0;
};

struct Correction {
  StateVector state = StateVector::Zero();
  StateMatrix covariance = StateMatrix::Zero();
  double innovation_v = 0.0;
  double innovation_variance_v2 = 0.0;
  double nis = 0.0;
  Eigen::RowVector4d measurement_jacobian = Eigen::RowVector4d::Zero();
};

class AugmentedEkf {
 public:
  AugmentedEkf(BatteryModel model, FilterConfig config);

  Prediction Predict(const StateVector& state,
                     const StateMatrix& covariance,
                     double measured_current_a,
                     double duration_s) const;
  Correction Correct(const Prediction& prediction,
                     double measured_voltage_v) const;

  const BatteryModel& model() const { return model_; }
  const FilterConfig& config() const { return config_; }

 private:
  BatteryModel model_;
  FilterConfig config_;
};

}  // namespace battery_soc267
