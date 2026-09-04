# isodistrreg: Rust crate

<p align="center"><strong>Available on <a href="https://crates.io/crates/isodistrreg">crates.io</a></strong></p>

Core Rust library for Isotonic Distributional Regression (IDR) and Survival-IDR
(S-IDR). See the [main README](https://github.com/AlexanderHenzi/isodistrreg)
for background and references.

## Features

| Feature | Default | Description |
|---------|---------|-------------|
| `partial-order` | off | Partial-order covariates (multivariate) |
| `subagging` | off | Multiple fits on subsets of the data |

## Usage

```toml
[dependencies]
isodistrreg = "0.6"
```

Enable optional features as needed:

```toml
[dependencies]
isodistrreg = { version = "0.6", features = ["partial-order"] }
```

## Precision

`IsotonicDistributionalRegressionFit` is generic over its covariate and response
scalar types via the associated types `X` and `Y` (`f32` or `f64`), which are
deduced from the slices passed to `fit`. Sample weights have their own
independent `Float` parameter `W`. CDF outputs are always `f32` — that is the
precision the algorithm body computes and stores.
