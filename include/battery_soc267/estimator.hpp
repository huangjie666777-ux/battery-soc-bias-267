#pragma once

#include <string>

#include "battery_soc267/ekf.hpp"
#include "battery_soc267/model.hpp"

namespace battery_soc267 {

// Public facade bundling model configuration, filter initialization and
// per-step execution for joint SOC / current-bias estimation.
class SocEstimator {
 public:
  SocEstimator() = default;

  // Configures and initializes the estimator. See AugmentedEkf::Initialize
  // for the validation rules. Returns false and fills "error" on invalid
  // input; the estimator then stays uninitialized.
  bool Initialize(const ModelParams& params, const StateVector& initial_state,
                  const StateCovariance& initial_covariance,
                  const ProcessNoise& process_noise,
                  double voltage_noise_variance, double nis_threshold,
                  std::string& error);

  // Convenience overload with an isotropic initial covariance.
  bool Initialize(const ModelParams& params, const StateVector& initial_state,
                  double initial_variance, const ProcessNoise& process_noise,
                  double voltage_noise_variance, double nis_threshold,
                  std::string& error);

  bool initialized() const { return initialized_; }

  // One estimation step. "terminal_voltage_meas" may be nullptr (prediction
  // only). On failure returns false, fills "error", and the estimate is
  // unchanged.
  bool Step(double dt, double current_meas,
            const double* terminal_voltage_meas, StepResult& result,
            std::string& error);

  // Overload with optional-style voltage flag.
  bool Step(double dt, double current_meas, bool has_voltage,
            double terminal_voltage_meas, StepResult& result,
            std::string& error);

  const StateVector& state() const { return ekf_.state(); }
  const StateCovariance& covariance() const { return ekf_.covariance(); }
  double soc() const { return ekf_.state()(kIndexSoc); }
  double bias() const { return ekf_.state()(kIndexBias); }

 private:
  AugmentedEkf ekf_;
  bool initialized_ = false;
};

}  // namespace battery_soc267

