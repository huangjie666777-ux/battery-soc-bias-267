#include "battery_soc267/model.hpp"

#include <cmath>

namespace battery_soc267 {
namespace {

bool IsPositiveFinite(double v) { return std::isfinite(v) && v > 0.0; }

// Index of the segment whose right-slope applies at "soc" (right segment at
// interior nodes, last segment at the final node).
std::size_t SegmentIndex(const ModelParams& params, double soc) {
  const auto& xs = params.ocv_soc;
  if (soc >= xs.back()) return xs.size() - 2;
  // First node uses the segment to its right as well.
  std::size_t lo = 0;
  while (lo + 1 < xs.size() - 1 && soc >= xs[lo + 1]) ++lo;
  return lo;
}

}  // namespace

bool ValidateParams(const ModelParams& params, std::string& error) {
  if (!IsPositiveFinite(params.capacity_ah)) {
    error = "capacity_ah must be positive and finite";
    return false;
  }
  if (!IsPositiveFinite(params.r0) || !IsPositiveFinite(params.r1) ||
      !IsPositiveFinite(params.c1) || !IsPositiveFinite(params.r2) ||
      !IsPositiveFinite(params.c2)) {
    error = "R0/R1/C1/R2/C2 must be positive and finite";
    return false;
  }
  const auto& xs = params.ocv_soc;
  const auto& ys = params.ocv_voltage;
  if (xs.size() < 2 || xs.size() != ys.size()) {
    error = "OCV table needs >= 2 nodes and matching SOC/voltage lengths";
    return false;
  }
  if (xs.front() < 0.0 || xs.back() > 1.0) {
    error = "OCV SOC nodes must cover a range within [0, 1]";
    return false;
  }
  for (std::size_t i = 0; i < xs.size(); ++i) {
    if (!std::isfinite(xs[i]) || !std::isfinite(ys[i])) {
      error = "OCV table entries must be finite";
      return false;
    }
    if (i > 0) {
      if (!(xs[i] > xs[i - 1])) {
        error = "OCV SOC nodes must be strictly increasing";
        return false;
      }
      if (ys[i] < ys[i - 1]) {
        error = "OCV voltages must be non-decreasing";
        return false;
      }
    }
  }
  return true;
}

double Ocv(const ModelParams& params, double soc) {
  const auto& xs = params.ocv_soc;
  const auto& ys = params.ocv_voltage;
  const std::size_t seg = SegmentIndex(params, soc);
  const double t = (soc - xs[seg]) / (xs[seg + 1] - xs[seg]);
  return ys[seg] + t * (ys[seg + 1] - ys[seg]);
}

double OcvDerivative(const ModelParams& params, double soc) {
  const auto& xs = params.ocv_soc;
  const auto& ys = params.ocv_voltage;
  const std::size_t seg = SegmentIndex(params, soc);
  return (ys[seg + 1] - ys[seg]) / (xs[seg + 1] - xs[seg]);
}

StateVector PredictState(const ModelParams& params, const StateVector& state,
                         double current_meas, double dt) {
  const double i_true = current_meas - state(kIndexBias);
  StateVector next;
  next(kIndexSoc) =
      state(kIndexSoc) - i_true * dt / (params.capacity_ah * 3600.0);
  const double tau1 = params.r1 * params.c1;
  const double tau2 = params.r2 * params.c2;
  const double a1 = std::exp(-dt / tau1);
  const double a2 = std::exp(-dt / tau2);
  next(kIndexRc1) = state(kIndexRc1) * a1 + i_true * params.r1 * (1.0 - a1);
  next(kIndexRc2) = state(kIndexRc2) * a2 + i_true * params.r2 * (1.0 - a2);
  next(kIndexBias) = state(kIndexBias);
  return next;
}

double TerminalVoltage(const ModelParams& params, const StateVector& state,
                       double current_meas) {
  const double i_true = current_meas - state(kIndexBias);
  return Ocv(params, state(kIndexSoc)) - params.r0 * i_true -
         state(kIndexRc1) - state(kIndexRc2);
}

MeasurementJacobian TerminalVoltageJacobian(const ModelParams& params,
                                            const StateVector& state) {
  MeasurementJacobian h;
  h(0, kIndexSoc) = OcvDerivative(params, state(kIndexSoc));
  h(0, kIndexRc1) = -1.0;
  h(0, kIndexRc2) = -1.0;
  h(0, kIndexBias) = params.r0;  // d(-R0*(i_meas - b))/db = +R0
  return h;
}

}  // namespace battery_soc267

