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
utilities.py
    └── model_setup.py
            └── biol_methods.py
                    └── model_dynamics.py
                                └── simulate_dynamics.py
```

Each module imports all symbols from the module above it via `from <module> import *`.

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
