# isodistrreg: Isotonic Distributional Regression

This repository provides Rust, Python and R implementations of Isotonic
Distributional Regression (IDR), including the Survival-IDR (S-IDR) extension
for right-censored outcomes as encountered in survival analysis.

## Language bindings

| Language | Path | Install |
|----------|------|---------|
| **Python** | [`bindings/python`](bindings/python/README.md) | `pip install isodistrreg` |
| **R** | [`bindings/R`](bindings/R/README.md) | `devtools::install_github("AlexanderHenzi/isodistrreg")` or CRAN |
| **Rust** | [`isodistrreg`](isodistrreg/README.md) | Add `isodistrreg` to your `Cargo.toml` |

## References

Henzi, A., Ziegel, J.F. and Gneiting, T. (2021). Isotonic distributional
regression. *J R Stat Soc Series B*, 83: 963–993.
<https://doi.org/10.1111/rssb.12450>

Henzi, A., Moesching, A. and Duembgen, L. (2022). Accelerating the
Pool-Adjacent-Violators Algorithm for Isotonic Distributional Regression.
*Methodol Comput Appl Probab*.
<https://doi.org/10.1007/s11009-022-09937-2>
