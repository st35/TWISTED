# `model_dynamics`

Vectorised dynamics over the whole DNA. These functions assemble the segment-level inputs, call the per-segment `biol_methods`, and produce the time derivatives that drive `solve_ivp`.

---

## State-vector helpers

### `get_state_vectors_from_dicts`

```python
get_state_vectors_from_dicts(model: Model) -> tuple[list[int], list[float]]
```

Flatten `model.x_dict`, `model.theta_dict`, `model.Lk` into the integrator state vector. Returns `(RNAP_gene_index, state_vector)`. RNAPs are ordered **right to left** along the DNA.

### `update_dicts_from_state_vector`

```python
update_dicts_from_state_vector(model, RNAP_gene_index, state_vector) -> None
```

Inverse of the above. Writes the integrator output back into `model.x_dict`, `model.theta_dict`, `model.Lk`.

---

## Segment attributes

### `calculate_segments_attributes`

```python
calculate_segments_attributes(
    model, RNAP_gene_index, state_vector,
) -> tuple[list[float], list[float], list[float], list[int], list[float], list[float]]
```

Returns `(segments_lengths, segments_sigmas, segments_torques, segments_dna_states, segments_writhe_fractions, segments_plectoneme_thresholds)`. For eukaryotic models, the per-segment nucleosome occupancy fraction is computed via [`get_nucleosome_occupied_fraction_per_segment`](utilities.md#get_nucleosome_occupied_fraction_per_segment) and passed into the eukaryotic torque law.

---

## Velocities and rates

### `get_RNAP_velocities`

```python
get_RNAP_velocities(model, state_vector, RNAP_gene_index, segments_lengths, segments_torques) -> list[float]
```

Per-RNAP linear velocity (nm/s) including the steric-hindrance soft ramp against the nearest obstacle in the direction of travel.

### `get_RNAP_angular_velocities`

```python
get_RNAP_angular_velocities(model, RNAP_gene_index, state_vector, segments_lengths, segments_torques, RNAP_velocities) -> list[float]
```

Per-RNAP angular velocity (rad/s).

### `get_segments_Lk_dynamics`

```python
get_segments_Lk_dynamics(model, RNAP_gene_index, state_vector, dx_dt, dtheta_dt, segments_lengths, segments_sigmas, segments_torques, segments_writhe_fractions) -> list[float]
```

Per-segment `dLk/dt`. Iterates over segments right-to-left, dispatches by left/right barrier type and clamp status.

### `get_RNAP_recruitment_rates`

```python
get_RNAP_recruitment_rates(model, RNAP_gene_index, state_vector, segments_lengths, segments_sigmas) -> list[float]
```

Per-gene recruitment rate. Currently ignores `TSS_sigma`.

### `get_mRNA_degradation_rates`

```python
get_mRNA_degradation_rates(model: Model) -> list[float]
```

Per-gene mRNA degradation rate.

### `get_binding_proteins_on_rates`

```python
get_binding_proteins_on_rates(model, segments_lengths, segments_sigmas) -> list[list[float]]
```

Per-species per-segment on-rates. Each species' total on-rate is `sum(per_segment[k])`.

### `get_binding_proteins_off_rates`

```python
get_binding_proteins_off_rates(model, segments_lengths, segments_sigmas) -> list[list[float]]
```

Per-species per-molecule off-rates. Each species' total off-rate is `sum(per_molecule[k])`.

---

## RNAP lifecycle

### `are_RNAPs_alive`

```python
are_RNAPs_alive(model, RNAP_gene_index, state_vector) -> list[int]
```

Per-RNAP boolean (1 alive, 0 finished). An RNAP is dead when its position has moved past `TSS + gene_length` (`+1` strand) or before `TSS − gene_length` (`−1` strand), modulo direction.

### `update_state_vector_to_remove_dead_RNAPs`

```python
update_state_vector_to_remove_dead_RNAPs(model, RNAP_gene_index, t, state_vector, simulation_setup_and_state) -> None
```

Removes finished RNAPs from `state_vector` and `RNAP_gene_index`, merges the two adjacent segments' linking numbers, increments the finished-transcription counters, and (if `mRNA_dynamics_mode == 1`) increments the mRNA count.

---

## Lk bookkeeping for discrete events

### `update_Lk_vector_after_RNAP_or_protein_recruitment`

```python
update_Lk_vector_after_RNAP_or_protein_recruitment(model, recruitment_position, RNAP_gene_index, state_vector, segments_lengths, segments_sigmas) -> None
```

Splits the segment that contains `recruitment_position` into two new sub-segments and partitions the parent `Lk` so that both sub-segments share the parent's `σ`.

### `update_Lk_vector_after_protein_unbinding`

```python
update_Lk_vector_after_protein_unbinding(model, unbinding_position, RNAP_gene_index, state_vector, segments_lengths, segments_sigmas) -> None
```

For a topological-barrier protein leaving: merges the segments on either side and adds their linking numbers.

---

## Event rate vector

### `get_events_rates`

```python
get_events_rates(model, RNAP_gene_index, state_vector) -> tuple[list[float], list[int]]
```

Returns `(rates_vector, events_indices)`. `events_indices` is a list of cumulative split points marking the boundaries between event blocks. See [Events and propensities](../user-guide/events-and-propensities.md) for the block layout.

---

## ODE right-hand side

### `model_dynamics`

```python
model_dynamics(t, state_vector, RNAP_gene_index, model, simulation_setup_and_state) -> list[float]
```

The function handed to `scipy.integrate.solve_ivp`. Returns the time derivatives of `state_vector`, where the last component is `dA/dt = sum(rates)`.

---

## Integration window

### `integrate`

```python
integrate(model, simulation_setup_and_state, t_start, state_vector, RNAP_gene_index, p0, print_at_each_integration_step) -> tuple[float, float]
```

Advances the ODE integrator until the cumulative propensity `A(t)` exceeds the event threshold `ln(1/p0)`, removing dead RNAPs in between chunks. Returns `(dt, a0)` — the elapsed time to the event and the propensity total at the event time.

See [Theory → Coupling continuous and discrete dynamics](../theory/dna-mechanics.md#10-coupling-continuous-and-discrete-dynamics).
