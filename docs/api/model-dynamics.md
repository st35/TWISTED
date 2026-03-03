# `model_dynamics` — Integration and State Management

**File:** `model_dynamics.py`
**Imports:** `model_setup`, `biol_methods`, `utilities`, `scipy`, `numpy`

This module manages the state vector, computes the ODE right-hand side, calculates event rates, and integrates the system forward in time.

---

## State Vector Conventions

RNAPs are ordered **right-to-left** in the state vector (descending x position). DNA segments are indexed from right to left: segment 0 is between the right clamp and the rightmost RNAP; segment `n` (last) is between the leftmost RNAP and the left clamp.

```
state_vector = [x₁, …, xₙ,          # positions (nm), right-to-left
                θ₁, …, θₙ,          # angular positions (rad)
                Lk₀, …, Lkₙ,       # linking numbers (n+1 segments)
                A]                   # cumulative propensity (dimensionless)
```

---

## State Vector Access

### `get_state_vectors_from_dicts`

```python
get_state_vectors_from_dicts(model: Model) -> tuple[list[int], list[float]]
```

Converts the per-gene `x_dict`/`theta_dict` dictionaries and `Lk` list into a flat state vector suitable for ODE integration. Returns `(RNAP_gene_index, state_vector)`.

`RNAP_gene_index[i]` is the gene index of the i-th RNAP in the state vector.

### `update_dicts_from_state_vector`

```python
update_dicts_from_state_vector(
    model: Model,
    RNAP_gene_index: list[int],
    state_vector: list[float]
) -> None
```

Writes the state vector back into the model's `x_dict`, `theta_dict`, and `Lk`. Reverses the ordering for negative-strand genes.

---

## Segment Attributes

### `calculate_segments_attributes`

```python
calculate_segments_attributes(
    model: Model,
    RNAP_gene_index: list[int],
    state_vector: list[float]
) -> tuple[
    list[float],  # lengths (nm)
    list[float],  # sigmas
    list[float],  # torques (pN·nm)
    list[int],    # DNA states
    list[float],  # writhe fractions
    list[float],  # plectoneme thresholds
]
```

Computes the geometric and mechanical properties of all DNA segments given the current state vector. Calls `get_prokaryotic_torque` for each segment.

---

## Velocity Calculations

### `get_RNAP_velocities`

```python
get_RNAP_velocities(
    model: Model,
    state_vector: list[float],
    RNAP_gene_index: list[int],
    segments_lengths: list[float],
    segments_torques: list[float]
) -> list[float]
```

Returns linear velocities (nm/s) for all RNAPs. Also zeroes velocities for RNAPs blocked by bound topoisomerases in `topoisomerase_based` mode.

### `get_RNAP_angular_velocities`

```python
get_RNAP_angular_velocities(
    model: Model,
    RNAP_gene_index: list[int],
    state_vector: list[float],
    segments_lengths: list[float],
    segments_torques: list[float],
    RNAP_velocities: list[float]
) -> list[float]
```

Returns angular velocities (rad/s) for all RNAPs using the distances from their respective TSSes.

### `get_segments_Lk_dynamics`

```python
get_segments_Lk_dynamics(
    model: Model,
    dtheta_dt: list[float],
    segments_lengths: list[float],
    segments_sigmas: list[float],
    segments_torques: list[float],
    segments_writhe_fractions: list[float]
) -> list[float]
```

Returns $\dot{Lk}$ (turns/s) for all segments. In `topoisomerase_based` mode, adds topoisomerase-driven contributions from `get_TOP1_effect_on_Lk_dynamics` and `get_TOP2_effect_on_Lk_dynamics`.

---

## Event Rate Calculations

### `get_RNAP_recruitment_rates`

```python
get_RNAP_recruitment_rates(
    model: Model,
    RNAP_gene_index: list[int],
    state_vector: list[float],
    segments_lengths: list[float],
    segments_sigmas: list[float]
) -> list[float]
```

Returns recruitment rates (s⁻¹) for each gene. A rate of 0 is returned if the TSS is sterically blocked or the promoter is OFF.

### `get_TOPO_binding_rates`

```python
get_TOPO_binding_rates(
    model: Model,
    segments_lengths: list[float],
    segments_sigmas: list[float]
) -> list[float]
```

Returns binding rates for all unbound topoisomerases. Returns an empty list for modes without explicit topoisomerase tracking.

### `get_TOPO_unbinding_rates`

```python
get_TOPO_unbinding_rates(
    model: Model,
    segments_lengths: list[float],
    segments_sigmas: list[float]
) -> list[float]
```

Returns unbinding rates for all bound topoisomerases.

### `get_events_rates`

```python
get_events_rates(
    model: Model,
    RNAP_gene_index: list[int],
    state_vector: list[float]
) -> tuple[list[float], list[int]]
```

Assembles the complete rates vector and returns `(rates_vector, events_indices)`. The rates vector layout:

```
[RNAP_recruitment_rates...]   # indices 0 to events_indices[0]-1
[model_observation_rate]      # index events_indices[0]
[global_relaxation_rate]      # index events_indices[1]
[local_relaxation_rates...]   # indices events_indices[2] to events_indices[3]-1
[TOPO_activity_rates...]      # indices events_indices[3] to events_indices[4]-1
[TOPO_binding_rates...]       # indices events_indices[4] to events_indices[5]-1
[TOPO_unbinding_rates...]     # indices events_indices[5] to events_indices[6]-1
```

`events_indices` is a 7-element list of cumulative boundary indices.

---

## ODE Right-Hand Side

### `model_dynamics`

```python
model_dynamics(
    t: float,
    state_vector: list[float],
    RNAP_gene_index: list[int],
    model: Model,
    simulation_setup_and_state: SimulationSetupAndState
) -> list[float]
```

The scipy-compatible ODE function. Returns the time derivative of the full state vector:

```
d/dt(state_vector) = [dx/dt..., dθ/dt..., dLk/dt..., da0/dt]
```

where `da0/dt = Σ rates_i(t)` is the instantaneous total propensity (used to track cumulative propensity for Gillespie timing).

---

## Integration

### `integrate`

```python
integrate(
    model: Model,
    simulation_setup_and_state: SimulationSetupAndState,
    t_start: float,
    state_vector: list[float],
    RNAP_gene_index: list[int],
    p0: float,
    print_at_each_integration_step: Callable | None
) -> tuple[float, float]
```

Integrates the ODE from `t_start` until the cumulative propensity increment satisfies:

$$\int_{t_{\mathrm{start}}}^{t_{\mathrm{event}}} a_0(t')\,dt' \geq \ln(1/p_0)$$

Uses `scipy.integrate.solve_ivp` with `method='RK45'`. Calls `update_state_vector_to_remove_dead_RNAPs` at each ODE output time point.

Returns `(dt, cumulative_propensity_at_event)` where `dt = t_event - t_start`.

---

## State Update Utilities

### `update_Lk_vector_after_RNAP_recruitment`

```python
update_Lk_vector_after_RNAP_recruitment(
    model: Model,
    TSS_index: int,
    RNAP_gene_index: list[int],
    state_vector: list[float],
    segments_lengths: list[float],
    segments_sigmas: list[float]
) -> None
```

Splits the DNA segment containing the TSS into two segments and assigns linking numbers such that the supercoiling density is conserved. Updates `model.Lk` in-place.

### `update_state_vector_to_remove_dead_RNAPs`

```python
update_state_vector_to_remove_dead_RNAPs(
    model: Model,
    RNAP_gene_index: list[int],
    t: float,
    state_vector: list[float],
    simulation_setup_and_state: SimulationSetupAndState
) -> None
```

Removes RNAPs that have passed the end of their gene. The linking numbers of the segments flanking a departing RNAP are summed into a single merged segment. Simultaneously records exit times, positions, and completion counts in `simulation_setup_and_state`. Modifies `RNAP_gene_index` and `state_vector` in-place.

### `are_RNAPs_alive`

```python
are_RNAPs_alive(
    model: Model,
    RNAP_gene_index: list[int],
    state_vector: list[float]
) -> list[int]
```

Returns a list of binary flags (`1` = active, `0` = finished) for all RNAPs. An RNAP is considered finished when it has moved beyond `TSS + gene_length` (positive strand) or below `TSS - gene_length` (negative strand).
