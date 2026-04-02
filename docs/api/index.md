# API Reference Overview

The TWISTED codebase is organized into five modules. Each module is documented on its own page.

| Module | File | Description |
|--------|------|-------------|
| [`model_setup`](model-setup.md) | `model_setup.py` | Core data classes |
| [`biol_methods`](biol-methods.md) | `biol_methods.py` | Biophysical equations |
| [`model_dynamics`](model-dynamics.md) | `model_dynamics.py` | ODE integration and state management |
| [`simulate_dynamics`](simulate-dynamics.md) | `simulate_dynamics.py` | Main simulation loop |
| [`utilities`](utilities.md) | `utilities.py` | I/O, stochastic helpers, geometry |

---

## Dependency Graph

```
utilities.py ←──────────────────────────┐
    └── model_setup.py                  │
            ├── biol_methods.py          │
            └── model_dynamics.py ───────┤
                    └── simulate_dynamics.py
```

All modules use `from <module> import *`. The primary chain is `utilities → model_setup → biol_methods → model_dynamics → simulate_dynamics`. In addition, `model_dynamics` and `simulate_dynamics` import directly from `utilities`, and `utilities` imports from `model_setup` (circular).

---

## Type Conventions

All physical quantities use the following units unless otherwise noted:

| Quantity | Unit |
|----------|------|
| Position | nm |
| Length | nm |
| Time | s |
| Rate | s⁻¹ |
| Torque | pN·nm |
| Force | pN |
| Velocity | nm/s |
| Angular velocity | rad/s |
| Linking number | dimensionless (turns) |
| Supercoiling density σ | dimensionless |
