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
isodistrreg = { git = "https://github.com/AlexanderHenzi/isodistrreg" }
```

Enable optional features as needed:

```toml
[dependencies]
isodistrreg = { git = "https://github.com/AlexanderHenzi/isodistrreg", features = ["partial-order"] }
```
