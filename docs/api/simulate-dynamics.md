# `simulate_dynamics` — Main Simulation Loop

**File:** `simulate_dynamics.py`
**Imports:** `model_setup`, `model_dynamics`, `utilities`

---

## `simulate_dynamics`

```python
simulate_dynamics(
    model: Model,
    simulation_setup_and_state: SimulationSetupAndState,
    print_at_each_integration_step: Callable | None = None,
    print_at_each_simulation_step: Callable | None = None,
    print_at_end_of_simulation: Callable | None = None
) -> None
```

Runs the main simulation loop. Modifies `model` and `simulation_setup_and_state` in-place. Returns `None`; all results are accessed from `simulation_setup_and_state` after the call.

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `model` | `Model` | The model to simulate |
| `simulation_setup_and_state` | `SimulationSetupAndState` | Controls termination and stores results |
| `print_at_each_integration_step` | `Callable or None` | Called at the start of each integration interval (approximately every `RNAP_alive_status_check_interval` seconds); receives `(model, sim, t, state_vector)` |
| `print_at_each_simulation_step` | `Callable or None` | Called at the start of each Gillespie step; receives `(model, sim)` |
| `print_at_end_of_simulation` | `Callable or None` | Called once when the termination criterion is met; receives `(model, sim)` |

### Algorithm

Each iteration of the loop:

1. Extract the current state vectors.
2. Draw `p0 ~ U(0,1)` for event-time sampling.
3. Integrate ODEs until the cumulative propensity satisfies the Gillespie criterion.
4. Update model dictionaries from the new state vector.
5. Recalculate segment attributes.
6. Draw `p1 ~ U(0,1)` and select an event by its propensity weight.
7. Execute the selected event (see table below).
8. Check termination condition.

### Event Dispatch Table

| Event range | Event type | Action |
|-------------|-----------|--------|
| `< events_indices[0]` | **RNAP recruitment** | Adds RNAP at TSS position; splits `Lk` vector; records recruitment time. Fails silently if TSS is sterically blocked or recruitment cap is reached |
| `events_indices[0]` | **Model observation** | No-op (dummy event to advance clock) |
| `events_indices[1]` | **Global supercoiling relaxation** | Resets `Lk` of all or one segment(s) to relaxed value; mode-dependent |
| `events_indices[2]–events_indices[3]` | **Type-specific relaxation** | Selectively resets positively or negatively supercoiled segments |
| `events_indices[3]–events_indices[4]` | **Topoisomerase approximated** | TOP1: reset `Lk` if no writhe; TOP2: reset to plectoneme threshold if writhe present |
| `events_indices[4]–events_indices[5]` | **mRNA degradation** | Decrements `model.mRNA_counts` for the selected gene by 1. Raises `ValueError` if the gene has zero mRNA (should not be selected when count is 0, as the rate is 0) |
| `events_indices[5]–events_indices[6]` | **Binding protein binding** | Selects a segment proportional to per-segment on-rates, places the protein at a random position within that segment, and appends to `model.binding_proteins_positions` |
| `events_indices[6]–events_indices[7]` | **Binding protein unbinding** | Selects which bound protein to release (proportional to per-molecule off-rates) and removes it from `model.binding_proteins_positions` |

### RNAP Recruitment Details

An RNAP recruitment event succeeds only when **all** of the following conditions hold:

1. The TSS is not sterically blocked by another RNAP within `RNAP_diameter` nm.
2. No bound nucleosome (with `is_steric_barrier_to_RNAPs=True`) is within `(RNAP_diameter + per_nucleosome_DNA_length + nucleosome_linker_length) / 2` nm of the TSS — **unless** the nucleosome has `can_be_displaced_at_TSS_by_RNAP=True` (see below).
3. No other bound protein (with `is_steric_barrier_to_RNAPs=True`) is within `(RNAP_diameter + generic_binding_protein_diameter) / 2` nm of the TSS — **unless** the protein has `can_be_displaced_at_TSS_by_RNAP=True` (see below).
4. The `max_RNAPs_to_recruit` cap for that gene has not been reached.

If recruitment fails, the event iteration is still consumed (time advances).

#### Protein Displacement at TSS

When a binding protein blocks the TSS but has `can_be_displaced_at_TSS_by_RNAP=True`, the RNAP recruitment **succeeds** and the blocking protein is removed from `model.binding_proteins_positions`. This models scenarios where an incoming RNAP can evict an obstacle — for example, nucleosome displacement during transcription initiation in eukaryotic chromatin. The displacement is controlled per protein type via the `can_be_displaced_at_TSS_by_RNAP` attribute on `BindingProtein`, and for auto-created nucleosomes via the `nucleosomes_can_be_displaced_at_TSS_by_RNAP` keyword on `GenomicSetup`.

### Supercoiling Relaxation Details

**`global_overall`:** All segment linking numbers are immediately reset to their relaxed values.

**`global_per_segment`:** One segment is chosen with probability proportional to its fractional length; its `Lk` is reset.

**`global_by_type`:** The first rate in `local_supercoiling_relaxation_rates` governs relaxation of all positively supercoiled segments; the second governs all negatively supercoiled segments. Both happen globally across all segments simultaneously.

**`per_segment_by_type`:** Same as `global_by_type` but one segment is chosen per event.

**`topoisomerase_approximated`:** A segment is chosen proportional to length. TOP1 resets `Lk = Lk₀` only if writhe is absent. TOP2 resets `Lk = Lk₀(1 + σ_s)` only if writhe is present (relaxes to the plectoneme-formation threshold, emulating strand passage).

**`topoisomerase_based`:** *Not yet implemented.*

### Termination

**Time-based (`simulation_end_mode=0`):** Loop exits when `curr_simulation_time >= simulation_end_time`.

**Event-based (`simulation_end_mode=1`):** Loop exits when `RNAPs_finished_transcription[i] >= simulation_end_event_counts[i]` for **all** genes `i`.

### Example Callbacks

```python
# Log every Gillespie step
def step_log(model, sim):
    n_rnaps = sum(len(x) for x in model.x_dict)
    print(f't={sim.curr_simulation_time:.2f}  active_RNAPs={n_rnaps}')

# Snapshot at end
def final_log(model, sim):
    print('Simulation complete.')
    print('Finished:', sim.RNAPs_finished_transcription)

simulate_dynamics(model, sim,
                  print_at_each_simulation_step=step_log,
                  print_at_end_of_simulation=final_log)
```
