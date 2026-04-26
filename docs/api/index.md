# API reference

TWISTED is laid out as a small flat package with one file per architectural layer. This page is the entry point to the function-by-function reference.

## Module map

| Module | Role |
|--------|------|
| [`model_setup`](model-setup.md) | Data classes the user constructs: `GenomicSetup`, `ModelSetup`, `BindingProtein`, `Model`, `SimulationSetupAndState` |
| [`biol_methods`](biol-methods.md) | Pure scalar physics/biology functions: torque laws, single-RNAP velocity, Lk dynamics, recruitment rate |
| [`model_dynamics`](model-dynamics.md) | Vectorised dynamics over the whole DNA: state vector ↔ dict conversions, segment attributes, ODE right-hand side, integration window |
| [`simulate_dynamics`](simulate-dynamics.md) | Outer Gillespie loop with full event dispatch |
| [`utilities`](utilities.md) | Pure helpers: gene-file IO, segment/spot lookup, steric checks, sampling helpers |

## Dependency graph

```
utilities  ──►  model_setup ──►  biol_methods ──►  model_dynamics ──►  simulate_dynamics
                                                            ▲
                                                            └─── utilities
```

Cycles are avoided: `biol_methods` depends only on `model_setup`; `model_dynamics` consumes both; `simulate_dynamics` orchestrates the whole.

## Unit conventions

| Quantity | Unit | Notes |
|----------|------|------|
| Length | nanometres (nm) | `1 bp = 0.34 nm` |
| Time | seconds (s) | |
| Rates | s⁻¹ | |
| Force | piconewtons (pN) | |
| Energy | pN·nm | `kBT ≈ 4.114 pN·nm` at 300 K |
| Torque | pN·nm | |
| Linking number | dimensionless | `h_dna = 2π / w₀ ≈ 10.5 bp/turn → 3.57 nm/turn` |

## Conventions for callable arguments

- All public dynamics functions take a `Model` plus arrays it does not own. The arrays live in either the user's calling code or the integrator state vector.
- Segments are always indexed **right to left** along the DNA. RNAPs are stored in the same order.
- The "state vector" handed to the SciPy ODE solver concatenates `[x_RNAPs, theta_RNAPs, Lk_segments, A]`, where `A` is the cumulative event propensity.
