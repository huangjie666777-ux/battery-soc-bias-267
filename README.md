# battery_soc267

`battery_soc267` is a small C++17 library for jointly estimating one cell's
SOC and constant current-sensor bias. It uses an augmented extended Kalman
filter with a second-order RC equivalent-circuit model and Eigen 3.4.0 from
`third_party/eigen3`. There is no frontend or HTTP component.

Sign convention: current is amperes, voltage is volts, time is seconds, and
discharge current is positive. The true current is `reading - b`.

## Model

State order:

```text
[SOC, V1, V2, b]
```

Terminal voltage:

```text
V = OCV(SOC) - R0 * (I_reading - b) - V1 - V2
```

The OCV table must include SOC nodes `0` and `1`, use strictly increasing SOC,
and contain non-decreasing voltages. Lookup is piecewise linear. At an interior
node the right segment derivative is used; the final node uses the left segment.
Capacity, `R0`, `R1`, `C1`, `R2`, and `C2` must all be positive and finite.

Prediction uses exact ampere-hour integration and exponential RC updates for a
constant current over a positive step duration. The bias is constant in
prediction. Process covariance adds a diagonal covariance equal to the four
non-negative process-noise densities multiplied by duration.

When voltage is absent, the step is prediction-only. When present, innovation
and NIS are calculated. A NIS above the configured positive threshold rejects
the voltage but retains prediction; otherwise an EKF correction is applied.
The covariance correction uses Joseph form. Predicted or corrected SOC outside
`[0, 1]`, invalid inputs, or non-finite numerical results raise
`battery_soc267::BatterySocError`; in `SocEstimator::Step` such failures do not
mutate the estimator state. SOC is never truncated.

## Build and self-test

```sh
make
make test
```

The self-test executable prints charge, discharge, zero-bias, missing-
measurement, abnormal-voltage, and invalid-input examples with inputs,
estimated states, predicted voltage, innovation, NIS, and measurement status.

## Minimal interface

```cpp
#include "battery_soc267/soc_estimator.hpp"

battery_soc267::ModelConfig model;
model.capacity_ah = 2.0;
model.r0 = 0.025;
model.r1 = 0.018;
model.c1 = 1200.0;
model.r2 = 0.012;
model.c2 = 9000.0;
model.ocv_table = {{0.0, 3.2}, {1.0, 4.2}};

battery_soc267::FilterConfig filter;
filter.process_noise_density << 1e-9, 1e-7, 1e-7, 1e-9;
filter.voltage_noise_variance = 1e-5;
filter.nis_threshold = 9.0;

battery_soc267::StateVector x;
x << 0.6, 0.0, 0.0, 0.0;
battery_soc267::StateMatrix P =
    battery_soc267::StateMatrix::Zero();
P.diagonal() << 1e-4, 1e-5, 1e-5, 1e-4;

battery_soc267::SocEstimator estimator(x, P, model, filter);
auto result = estimator.Step({1.0, 2.0, 3.72});
```
