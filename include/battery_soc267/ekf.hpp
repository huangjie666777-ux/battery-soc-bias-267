#pragma once

#include "battery_soc267/model.hpp"

namespace battery_soc267 {

// Continuous-time process noise spectral densities, one per state component
// (SOC, v_rc1, v_rc2, bias). All entries must be non-negative and finite.
struct ProcessNoise {
  double soc = 0.0;
  double rc1 = 0.0;
  double rc2 = 0.0;
  double bias = 0.0;

  Eigen::Matrix<double, 4, 1> Diagonal() const {
    return {soc, rc1, rc2, bias};
  }
};

// Outcome of a single filter step.
enum class MeasurementStatus {
  kNoMeasurement,  // prediction only (no terminal voltage supplied)
  kAccepted,       // measurement passed the NIS gate and was fused
  kRejected,       // NIS exceeded the threshold; prediction kept, update skipped
};

struct StepResult {
  StateVector state = StateVector::Zero();
  StateCovariance covariance = StateCovariance::Zero();
  double predicted_voltage = 0.0;  // model terminal voltage for this step [V]
  MeasurementStatus status = MeasurementStatus::kNoMeasurement;
  double innovation = 0.0;  // measured - predicted voltage [V], if measured
  double nis = 0.0;         // normalized innovation squared, if measured
};

// Augmented extended Kalman filter over [soc, v_rc1, v_rc2, bias].
//
// Prediction uses the exact analytical discretization (see PredictState) and
// adds diag(process_noise) * dt to the covariance. The measurement update
// uses the terminal-voltage model with the full Jacobian (including the bias
// coupling), an NIS gate, and the Joseph covariance form.
//
// All operations are transactional: on any error the filter state and
// covariance are left untouched and Step() returns false.
class AugmentedEkf {
 public:
  // Initializes the filter. All noise densities must be non-negative,
  // voltage_noise_variance and nis_threshold positive, the initial covariance
  // symmetric positive definite, and the initial SOC within the OCV table.
  bool Initialize(const ModelParams& params, const StateVector& initial_state,
                  const StateCovariance& initial_covariance,
                  const ProcessNoise& process_noise,
                  double voltage_noise_variance, double nis_threshold,
                  std::string& error);

  // Runs one step: dt > 0 [s], constant measured current [A], and an optional
  // terminal voltage [V]. On success fills "result" and returns true; on
  // failure fills "error", returns false, and the filter is unchanged.
  bool Step(double dt, double current_meas,
            const double* terminal_voltage_meas, StepResult& result,
            std::string& error);

  const StateVector& state() const { return state_; }
  const StateCovariance& covariance() const { return covariance_; }

 private:
  ModelParams params_;
  StateVector state_ = StateVector::Zero();
  StateCovariance covariance_ = StateCovariance::Zero();
  ProcessNoise process_noise_;
  double voltage_noise_variance_ = 0.0;
  double nis_threshold_ = 0.0;
};

}  // namespace battery_soc267

