#pragma once

#include <Eigen/Dense>

#include <string>
#include <vector>

namespace battery_soc267 {

// Second-order RC equivalent-circuit model with constant parameters.
// Units: current in A (discharge positive), voltage in V, time in s.
struct ModelParams {
  double capacity_ah = 0.0;  // usable capacity [Ah], positive finite
  double r0 = 0.0;           // ohmic resistance [Ohm]
  double r1 = 0.0;           // polarization resistance 1 [Ohm]
  double c1 = 0.0;           // polarization capacitance 1 [F]
  double r2 = 0.0;           // polarization resistance 2 [Ohm]
  double c2 = 0.0;           // polarization capacitance 2 [F]
  std::vector<double> ocv_soc;      // SOC breakpoints in [0,1], strictly increasing
  std::vector<double> ocv_voltage;  // OCV values [V], non-decreasing
};

// Validates the parameter set. Returns true on success, otherwise fills
// "error" with a human-readable message and returns false.
bool ValidateParams(const ModelParams& params, std::string& error);

// Piecewise-linear OCV interpolation. "soc" must lie within the table range.
double Ocv(const ModelParams& params, double soc);

// OCV derivative dV/dSOC. At interior nodes the slope of the segment to the
// right is used; at the last node the slope of the segment to the left.
double OcvDerivative(const ModelParams& params, double soc);

// State of the augmented filter: [soc, v_rc1, v_rc2, current_bias].
// The true current is (measured current - bias).
using StateVector = Eigen::Matrix<double, 4, 1>;
using StateCovariance = Eigen::Matrix<double, 4, 4>;
using MeasurementJacobian = Eigen::Matrix<double, 1, 4>;

constexpr int kIndexSoc = 0;
constexpr int kIndexRc1 = 1;
constexpr int kIndexRc2 = 2;
constexpr int kIndexBias = 3;

// Exact (analytical) discrete-time prediction over "dt" seconds with constant
// measured current "current_meas" [A]:
//   soc'  = soc - i_true * dt / (capacity_ah * 3600)
//   v_k'  = v_k * exp(-dt/tau_k) + i_true * R_k * (1 - exp(-dt/tau_k))
//   bias' = bias
// where i_true = current_meas - bias and tau_k = R_k * C_k.
StateVector PredictState(const ModelParams& params, const StateVector& state,
                         double current_meas, double dt);

// Terminal voltage [V]: Vt = OCV(soc) - R0 * i_true - v_rc1 - v_rc2.
double TerminalVoltage(const ModelParams& params, const StateVector& state,
                       double current_meas);

// Jacobian of TerminalVoltage with respect to the state, including the bias
// coupling: dVt/dbias = +R0 (since i_true = i_meas - bias).
MeasurementJacobian TerminalVoltageJacobian(const ModelParams& params,
                                            const StateVector& state);

}  // namespace battery_soc267

