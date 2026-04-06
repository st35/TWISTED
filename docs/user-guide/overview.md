# Overview

TWISTED simulates the coupled dynamics of RNA polymerase (RNAP) elongation and DNA supercoiling on a one-dimensional DNA molecule. The simulation couples:

1. **Continuous mechanics** — RNAP positions, angular positions, and DNA linking numbers evolve according to ODEs that are integrated with an adaptive RK45 solver.
2. **Discrete stochastic events** — RNAP recruitment, supercoiling relaxation, and topoisomerase binding/unbinding are selected stochastically using a Gillespie-style algorithm.

---

## Architecture

```
┌─────────────────────────────────────────────────┐
│                  simulate_dynamics               │
│  (main simulation loop — simulate_dynamics.py)  │
└───────────────┬─────────────────────────────────┘
                │ calls
     ┌──────────┴──────────┐
     │                     │
  integrate()         event handling
  (model_dynamics.py)  (simulate_dynamics.py)
     │
     ▼
  model_dynamics()   ←── biol_methods.py
  (ODE right-hand side)
```

The central object is `Model`, which holds all dynamic state. `GenomicSetup` and `ModelSetup` provide static configuration (gene geometry and physical parameters, respectively).

---

## Simulation Loop

Each iteration of the main loop performs:

1. **Extract state vectors** from the model's internal dictionaries (`get_state_vectors_from_dicts`).
2. **Draw a random number** `p0 ∈ (0, 1)` for event-time sampling.
3. **Integrate** the ODEs forward until the cumulative propensity satisfies `A(t) > ln(1/p0)` (`integrate`). Completed RNAPs that have traversed their gene are removed during integration (`update_state_vector_to_remove_dead_RNAPs`).
4. **Update model dictionaries** from the new state vector.
5. **Recalculate segment attributes** (lengths, supercoiling densities, torques, DNA state).
6. **Draw a second random number** `p1 ∈ (0, 1)` for event-type selection.
7. **Select and execute an event** from the propensity-weighted distribution.
8. **Check termination condition** (time-based or event-count-based).

---

## State Representation

The DNA is divided into segments separated by RNAPs. RNAPs are ordered right-to-left internally. The state vector has the form:

```
state_vector = [x₁, x₂, …, xₙ,          # RNAP positions (nm)
                θ₁, θ₂, …, θₙ,          # RNAP angular positions (rad)
                Lk₀, Lk₁, …, Lkₙ,       # Linking numbers of n+1 segments
                A]                        # Cumulative propensity
```

where `n` is the current number of active RNAPs and there are `n+1` DNA segments (the segments between consecutive RNAPs and the clamps).

---

## Module Summary

| Module | Responsibility |
|--------|---------------|
| `model_setup.py` | Data classes: `GenomicSetup`, `ModelSetup`, `BindingProtein`, `Model`, `SimulationSetupAndState` |
| `biol_methods.py` | Physical equations: RNAP velocity/angular velocity, torque, topoisomerase effects |
| `model_dynamics.py` | State vector management, ODE right-hand side, event rate calculation, integration |
| `simulate_dynamics.py` | Main simulation loop and event dispatch |
| `utilities.py` | I/O, stochastic helpers, steric-hindrance checks |

---

## Data Flow

```
GenomicSetup + ModelSetup
          │
          ▼
        Model  ──────────────────────────────────────────┐
          │                                              │
          ▼                                              │
  get_state_vectors_from_dicts()                         │
          │                                              │
          ▼                                              │
  integrate() ──► model_dynamics() ──► biol_methods      │
          │                                              │
          ▼                                              │
  event selected & executed ──────────────────────────►  │
          │                                              │
          ▼                                              │
  update_dicts_from_state_vector()  ◄────────────────────┘
```
