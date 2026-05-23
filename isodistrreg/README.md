# isodistrreg: Rust crate

Core Rust library for Isotonic Distributional Regression (IDR) and Survival-IDR
(S-IDR). See the [main README](https://github.com/AlexanderHenzi/isodistrreg)
for background and references.

## Features

| Feature | Default | Description |
|---------|---------|-------------|
| `partial-order` | off | Partial-order covariates (requires OSQP) |
| `subagging` | off | Subagging for partial-order IDR |

## Usage

```toml
[dependencies]
isodistrreg = "0.4"
```

Enable optional features as needed:

```toml
[dependencies]
isodistrreg = { version = "0.4", features = ["partial-order"] }
```

## Precision

`IsotonicDistributionalRegressionFit` is generic over its covariate and response
scalar types via the associated types `X` and `Y` (`f32` or `f64`), which are
deduced from the slices passed to `fit`. Sample weights have their own
independent `Float` parameter `W`. CDF outputs are always `f32` — that is the
precision the algorithm body computes and stores.
