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
