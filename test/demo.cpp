#include "battery_soc267/soc_estimator.hpp"

#include <cassert>
#include <iomanip>
#include <iostream>
#include <string>

using battery_soc267::BatteryModel;
using battery_soc267::BatterySocError;
using battery_soc267::FilterConfig;
using battery_soc267::MeasurementStatus;
using battery_soc267::ModelConfig;
using battery_soc267::SocEstimator;
using battery_soc267::StateMatrix;
using battery_soc267::StateVector;
using battery_soc267::StepInput;
using battery_soc267::StepResult;

namespace {

ModelConfig MakeModelConfig() {
  ModelConfig config;
  config.capacity_ah = 2.0;
  config.r0 = 0.025;
  config.r1 = 0.018;
  config.c1 = 1200.0;
  config.r2 = 0.012;
  config.c2 = 9000.0;
  const std::vector<double> soc = {0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0};
  const std::vector<double> voltage = {3.20, 3.42, 3.60, 3.72,
                                       3.86, 4.05, 4.20};
  for (std::size_t i = 0; i < soc.size(); ++i) {
    config.ocv_table.push_back({soc[i], voltage[i]});
  }
  return config;
}

FilterConfig MakeFilterConfig() {
  FilterConfig config;
  config.process_noise_density << 1.0e-9, 1.0e-7, 1.0e-7, 1.0e-9;
  config.voltage_noise_variance = 1.0e-5;
  config.nis_threshold = 9.0;
  return config;
}

StateMatrix MakeCovariance() {
  StateMatrix covariance;
  covariance.setZero();
  covariance.diagonal() << 1.0e-4, 1.0e-5, 1.0e-5, 1.0e-4;
  return covariance;
}

std::string StatusName(MeasurementStatus status) {
  switch (status) {
    case MeasurementStatus::kPredictedOnly:
      return "prediction-only";
    case MeasurementStatus::kUpdated:
      return "updated";
    case MeasurementStatus::kRejected:
      return "rejected";
  }
  return "unknown";
}

void PrintStep(const std::string& label,
               const StepInput& input,
               const StepResult& result) {
  std::cout << label << "\n"
            << std::fixed << std::setprecision(6)
            << "  input: dt=" << input.duration_s << " s, I="
            << input.current_a << " A";
  if (input.terminal_voltage_v.has_value()) {
    std::cout << ", V=" << *input.terminal_voltage_v << " V";
  } else {
    std::cout << ", V=<absent>";
  }
  std::cout << "\n  state: SOC=" << result.state[0]
            << ", V1=" << result.state[1] << ", V2=" << result.state[2]
            << ", b=" << result.state[3] << " A\n"
            << "  predicted V=" << result.predicted_voltage_v
            << " V, status=" << StatusName(result.measurement_status);
  if (result.innovation_v.has_value()) {
    std::cout << ", innovation=" << *result.innovation_v << " V";
  }
  if (result.nis.has_value()) {
    std::cout << ", NIS=" << *result.nis;
  }
  std::cout << "\n\n";
}

void SimulateObservedStep(BatteryModel& model,
                          StateVector& truth,
                          double measured_current,
                          double dt) {
  truth = model.PredictState(truth, measured_current, dt);
}

}  // namespace

int main() {
  const ModelConfig model_config = MakeModelConfig();
  BatteryModel standalone_model(model_config);
  assert(std::abs(standalone_model.OpenCircuitVoltage(0.05) - 3.31) < 1e-12);
  assert(std::abs(standalone_model.OpenCircuitVoltageDerivative(0.05) - 2.2) <
         1e-12);
  assert(std::abs(standalone_model.OpenCircuitVoltageDerivative(1.0) - 1.5) <
         1e-12);

  StateVector truth;
  truth << 0.60, 0.0, 0.0, 0.02;
  StateVector initial;
  initial << 0.60, 0.0, 0.0, 0.0;
  SocEstimator estimator(initial, MakeCovariance(), model_config,
                         MakeFilterConfig());
  BatteryModel simulation_model(model_config);

  std::cout << std::boolalpha;
  std::cout << "battery_soc267 self-test and examples\n\n";

  const double dt = 10.0;
  const double discharge_measured = 2.0 + truth[3];
  SimulateObservedStep(simulation_model, truth, discharge_measured, dt);
  double voltage = simulation_model.TerminalVoltage(truth,
                                                    discharge_measured);
  StepResult discharge = estimator.Step({dt, discharge_measured, voltage});
  PrintStep("discharge with biased current sensor",
            {dt, discharge_measured, voltage}, discharge);
  assert(discharge.measurement_status == MeasurementStatus::kUpdated);

  const double charge_measured = -1.5 + truth[3];
  SimulateObservedStep(simulation_model, truth, charge_measured, dt);
  voltage = simulation_model.TerminalVoltage(truth, charge_measured);
  StepResult charge = estimator.Step({dt, charge_measured, voltage});
  PrintStep("charge with biased current sensor",
            {dt, charge_measured, voltage}, charge);
  assert(charge.measurement_status == MeasurementStatus::kUpdated);

  const StateVector before_missing = estimator.state();
  StepResult missing = estimator.Step({dt, charge_measured, std::nullopt});
  PrintStep("missing voltage measurement: prediction only",
            {dt, charge_measured, std::nullopt}, missing);
  assert(missing.measurement_status == MeasurementStatus::kPredictedOnly);
  assert(!missing.innovation_v.has_value());
  assert(!missing.nis.has_value());
  assert((estimator.state().array() - missing.state.array()).abs().maxCoeff() <
         1e-15);
  (void)before_missing;

  const double abnormal_voltage = voltage - 0.5;
  StepResult rejected = estimator.Step({dt, charge_measured,
                                        abnormal_voltage});
  PrintStep("abnormal voltage: measurement rejected",
            {dt, charge_measured, abnormal_voltage}, rejected);
  assert(rejected.measurement_status == MeasurementStatus::kRejected);
  assert(rejected.nis.value() > MakeFilterConfig().nis_threshold);

  StateVector zero_bias_initial;
  zero_bias_initial << 0.80, 0.0, 0.0, 0.0;
  SocEstimator zero_bias_estimator(zero_bias_initial, MakeCovariance(),
                                  model_config, MakeFilterConfig());
  StateVector zero_bias_truth = zero_bias_initial;
  const double zero_bias_current = 1.0;
  zero_bias_truth = simulation_model.PredictState(
      zero_bias_truth, zero_bias_current, dt);
  voltage = simulation_model.TerminalVoltage(zero_bias_truth,
                                             zero_bias_current);
  StepResult zero_bias = zero_bias_estimator.Step(
      {dt, zero_bias_current, voltage});
  PrintStep("zero current bias", {dt, zero_bias_current, voltage}, zero_bias);
  assert(zero_bias.measurement_status == MeasurementStatus::kUpdated);
  assert(std::abs(zero_bias.state[3]) < 1.0e-3);

  const StateVector state_before_error = estimator.state();
  const StateMatrix covariance_before_error = estimator.covariance();
  bool caught_invalid_duration = false;
  try {
    estimator.Step({0.0, 1.0, std::nullopt});
  } catch (const BatterySocError& error) {
    caught_invalid_duration = true;
    std::cout << "invalid input error: " << error.what() << "\n";
  }
  assert(caught_invalid_duration);
  assert(estimator.state().isApprox(state_before_error));
  assert(estimator.covariance().isApprox(covariance_before_error));

  bool caught_boundary_error = false;
  try {
    estimator.Step({100000.0, 100.0, std::nullopt});
  } catch (const BatterySocError& error) {
    caught_boundary_error = true;
    std::cout << "SOC boundary error: " << error.what() << "\n\n";
  }
  assert(caught_boundary_error);
  assert(estimator.state().isApprox(state_before_error));
  assert(estimator.covariance().isApprox(covariance_before_error));

  std::cout << "all self-tests passed\n";
  return 0;
}
