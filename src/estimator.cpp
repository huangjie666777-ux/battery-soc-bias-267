#include "battery_soc267/estimator.hpp"

#include <cmath>

namespace battery_soc267 {

bool SocEstimator::Initialize(const ModelParams& params,
                              const StateVector& initial_state,
                              const StateCovariance& initial_covariance,
                              const ProcessNoise& process_noise,
                              double voltage_noise_variance,
                              double nis_threshold, std::string& error) {
  AugmentedEkf ekf;
  if (!ekf.Initialize(params, initial_state, initial_covariance,
                      process_noise, voltage_noise_variance, nis_threshold,
                      error)) {
    return false;
  }
  ekf_ = ekf;
  initialized_ = true;
  return true;
}

bool SocEstimator::Initialize(const ModelParams& params,
                              const StateVector& initial_state,
                              double initial_variance,
                              const ProcessNoise& process_noise,
                              double voltage_noise_variance,
                              double nis_threshold, std::string& error) {
  if (!(initial_variance > 0.0) || !std::isfinite(initial_variance)) {
    error = "initial variance must be positive and finite";
    return false;
  }
  return Initialize(params, initial_state,
                    StateCovariance::Identity() * initial_variance,
                    process_noise, voltage_noise_variance, nis_threshold,
                    error);
}

bool SocEstimator::Step(double dt, double current_meas,
                        const double* terminal_voltage_meas,
                        StepResult& result, std::string& error) {
  if (!initialized_) {
    error = "estimator not initialized";
    return false;
  }
  return ekf_.Step(dt, current_meas, terminal_voltage_meas, result, error);
}

bool SocEstimator::Step(double dt, double current_meas, bool has_voltage,
                        double terminal_voltage_meas, StepResult& result,
                        std::string& error) {
  return Step(dt, current_meas,
              has_voltage ? &terminal_voltage_meas : nullptr, result, error);
}

}  // namespace battery_soc267
