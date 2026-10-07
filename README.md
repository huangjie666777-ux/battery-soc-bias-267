# battery_soc267

Single-cell SOC and current-bias joint estimation library for energy-storage
controllers. C++17, Eigen 3.4.0 (vendored in `third_party/eigen3`), built
with g++ and Make. No front end, no HTTP.

## Model

Second-order RC equivalent circuit with constant parameters, augmented with a
constant current-sensor bias. Units: A, V, s; discharge current is positive.

- State: `x = [soc, v_rc1, v_rc2, b]`, true current `i = i_meas - b`.
- Prediction (exact analytical discretization over `dt`):
  - `soc' = soc - i * dt / (capacity_ah * 3600)`
  - `v_k' = v_k * exp(-dt/tau_k) + i * R_k * (1 - exp(-dt/tau_k))`, `tau_k = R_k C_k`
  - `b' = b`
  - `P += diag(q) * dt` with the four non-negative process-noise densities.
- Measurement: `Vt = OCV(soc) - R0 * i - v_rc1 - v_rc2`.
  - OCV: piecewise-linear table over SOC in [0, 1], non-decreasing, strictly
    increasing SOC nodes. Node derivatives use the right segment; the last
    node uses the left segment.
  - Jacobian `H = [dOCV/dsoc, -1, -1, +R0]` includes the bias coupling, so
    voltage information corrects the bias state, not only the SOC.

## Filter

Augmented EKF (`src/ekf.cpp`):

- Missing terminal voltage -> prediction only.
- With a voltage: NIS = innovation^2 / innovation variance. If NIS exceeds
  the configured threshold the measurement is rejected and the prediction is
  kept; otherwise the state is corrected and the covariance updated in Joseph
  form.
- Any illegal input, numerical failure, or predicted/corrected SOC outside
  the OCV table range makes the step fail with an error message; the filter
  state is left unchanged (no truncation/clamping).

## Layout

- `include/battery_soc267/model.hpp`, `src/model.cpp` — parameters, OCV
  table, prediction and measurement equations with Jacobians.
- `include/battery_soc267/ekf.hpp`, `src/ekf.cpp` — augmented EKF
  (prediction, NIS gating, Joseph update).
- `include/battery_soc267/estimator.hpp`, `src/estimator.cpp` — public
  facade (`SocEstimator`).
- `examples/demo.cpp` — charge/discharge run with sensor bias, missing
  measurements, outlier rejection and invalid-input handling, printing
  inputs, estimates and diagnostics.
- `tests/selftest.cpp` — self-tests for the model, filter and interface.

## Build and run

```sh
make            # builds build/libbattery_soc267.a, selftest and demo
make test       # runs the self-tests
make demo       # runs the example
```

## Usage sketch

```cpp
battery_soc267::ModelParams params;   // fill capacity, R0/R1/C1/R2/C2, OCV table
battery_soc267::StateVector x0;       // [soc, v1, v2, bias]
battery_soc267::ProcessNoise q;       // four non-negative densities
battery_soc267::SocEstimator est;
std::string error;
est.Initialize(params, x0, /*initial_variance=*/1e-2, q,
               /*voltage_noise_variance=*/1e-4, /*nis_threshold=*/9.0, error);

battery_soc267::StepResult r;
// With a voltage measurement:
est.Step(/*dt=*/1.0, /*current_meas=*/10.0, /*has_voltage=*/true,
         /*voltage=*/3.71, r, error);
// Prediction only:
est.Step(1.0, 10.0, /*has_voltage=*/false, 0.0, r, error);
```

`StepResult` carries the state, covariance, predicted terminal voltage,
measurement status (none/accepted/rejected) and, when measured, the
innovation and NIS.

