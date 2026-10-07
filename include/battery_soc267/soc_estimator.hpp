#pragma once

#include "battery_soc267/augmented_ekf.hpp"

namespace battery_soc267 {

class SocEstimator {
 public:
  SocEstimator(StateVector initial_state,
               StateMatrix initial_covariance,
               ModelConfig model_config,
               FilterConfig filter_config);

  StepResult Step(const StepInput& input);

  const StateVector& state() const { return state_; }
  const StateMatrix& covariance() const { return covariance_; }
  const BatteryModel& model() const { return ekf_.model(); }

 private:
  static void ValidateInitialState(const StateVector& state);
  static void ValidateCovariance(const StateMatrix& covariance);
  static void ValidateStepInput(const StepInput& input);
  static void CheckFiniteState(const StateVector& state,
                               const StateMatrix& covariance,
                               const char* phase);

  AugmentedEkf ekf_;
  StateVector state_;
  StateMatrix covariance_;
};

}  // namespace battery_soc267
