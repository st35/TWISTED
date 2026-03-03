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
| `print_at_each_integration_step` | `Callable or None` | Called at every ODE evaluation point; receives `(model, sim, t, state_vector)` |
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
| `events_indices[4]–events_indices[5]` | **Topoisomerase binding** | Selects segment and position; binds topoisomerase if not sterically blocked |
| `events_indices[5]–events_indices[6]` | **Topoisomerase unbinding** | Releases bound topoisomerase; resets position and segment index to `-1` |
| `events_indices[6]–events_indices[7]` | **mRNA degradation** | Decrements `model.mRNA_counts` for the selected gene by 1. Raises `ValueError` if the gene has zero mRNA (should not be selected when count is 0, as the rate is 0) |

### RNAP Recruitment Details

An RNAP recruitment event succeeds only when **all** of the following conditions hold:

1. The TSS is not sterically blocked by another RNAP within `between_RNAPs_steric_effect_cutoff` nm.
2. No bound topoisomerase is within `RNAP_TOPO_steric_effect_cutoff` nm of the TSS (in `topoisomerase_based` mode).
3. The `max_RNAPs_to_recruit` cap for that gene has not been reached.

If recruitment fails, the event iteration is still consumed (time advances).

### Supercoiling Relaxation Details

**`global_overall`:** All segment linking numbers are immediately reset to their relaxed values.

**`global_per_segment`:** One segment is chosen with probability proportional to its fractional length; its `Lk` is reset.

**`global_by_type`:** The first rate in `local_supercoiling_relaxation_rates` governs relaxation of all positively supercoiled segments; the second governs all negatively supercoiled segments. Both happen globally across all segments simultaneously.

**`per_segment_by_type`:** Same as `global_by_type` but one segment is chosen per event.

**`topoisomerase_approximated`:** A segment is chosen proportional to length. TOP1 resets `Lk = Lk₀` only if writhe is absent. TOP2 resets `Lk = Lk₀(1 + σ_s)` only if writhe is present (relaxes to the plectoneme-formation threshold, emulating strand passage).

**`topoisomerase_based`:** Topoisomerase binding positions are drawn uniformly within the chosen segment; binding is refused if an RNAP is within the steric cutoff.

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
