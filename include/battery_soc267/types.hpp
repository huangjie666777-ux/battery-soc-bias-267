#pragma once

#include <Eigen/Dense>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace battery_soc267 {

inline constexpr int kStateSize = 4;

using StateVector = Eigen::Matrix<double, kStateSize, 1>;
using StateMatrix = Eigen::Matrix<double, kStateSize, kStateSize>;
using ProcessNoiseVector = Eigen::Matrix<double, 4, 1>;

class BatterySocError : public std::runtime_error {
 public:
  explicit BatterySocError(const std::string& message)
      : std::runtime_error(message) {}
};

struct OcvPoint {
  double soc = 0.0;
  double voltage = 0.0;
};

struct ModelConfig {
  double capacity_ah = 0.0;
  double r0 = 0.0;
  double r1 = 0.0;
  double c1 = 0.0;
  double r2 = 0.0;
  double c2 = 0.0;
  std::vector<OcvPoint> ocv_table;
};

struct FilterConfig {
  ProcessNoiseVector process_noise_density = ProcessNoiseVector::Zero();
  double voltage_noise_variance = 0.0;
  double nis_threshold = 0.0;
};

struct StepInput {
  double duration_s = 0.0;
  double current_a = 0.0;
  std::optional<double> terminal_voltage_v;
};

enum class MeasurementStatus {
  kPredictedOnly,
  kUpdated,
  kRejected
};

struct StepResult {
  StateVector state = StateVector::Zero();
  StateMatrix covariance = StateMatrix::Zero();
  double predicted_voltage_v = 0.0;
  MeasurementStatus measurement_status = MeasurementStatus::kPredictedOnly;
  std::optional<double> innovation_v;
  std::optional<double> nis;
};

}  // namespace battery_soc267
