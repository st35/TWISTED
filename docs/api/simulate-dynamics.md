# `simulate_dynamics`

## `simulate_dynamics`

```python
simulate_dynamics(
    model: Model,
    simulation_setup_and_state: SimulationSetupAndState,
    print_at_each_integration_step: Callable | None = None,
    print_at_each_simulation_step: Callable | None = None,
    print_at_end_of_simulation: Callable | None = None,
) -> None
```

The outer Gillespie loop. Repeatedly:

1. Builds the state vector from the current `model` dicts.
2. Draws a uniform `p0 ~ U(0, 1)` and calls [`integrate`](model-dynamics.md#integrate) to advance the ODE until the cumulative propensity passes `ln(1/p0)`.
3. Asks [`get_events_rates`](model-dynamics.md#get_events_rates) for the propensity vector and event-block split points.
4. Selects an event index proportional to its current rate via [`select_event_based_on_propensities`](utilities.md#select_event_based_on_propensities).
5. Dispatches to the appropriate handler.
6. Updates `model` from the state vector.
7. Calls `print_at_each_simulation_step(model, sim)` if provided.
8. Repeats until the termination criterion is met.

When the loop exits, `simulation_completed = True` and `print_at_end_of_simulation(model, sim)` is called.

### Termination

- `simulation_end_mode == 0` (time-based): exit when `curr_simulation_time >= simulation_end_time`.
- `simulation_end_mode == 1` (event-based): exit when `RNAPs_finished_transcription[i] >= simulation_end_event_counts[i]` for every gene `i`.

### Event dispatch table

The selected `event_index` falls into one of the cumulative blocks defined by `events_indices`. The dispatch is, in order:

| Block | Handler action |
|------|---------------|
| 0 — RNAP recruitment (per gene) | TSS steric check; if clear, append RNAP, split host segment, log recruitment time; if blocked by displaceable protein, remove that protein first |
| 1 — Observation event | No state change |
| 2 — Global SC relaxation | Reset every segment's `Lk` (`'global_overall'`) or one length-weighted segment's `Lk` (`'global_per_segment'`) |
| 3 — By-type SC relaxation | Reset every (or one length-weighted) segment of the appropriate sign |
| 4 — Topoisomerase-approximated | TOP1 event resets `Lk → Lk₀` on a writhe-free length-weighted segment; TOP2 event sets `Lk → Lk₀(1+σ_s)` on a plectonemic length-weighted segment |
| 5 — mRNA degradation (per gene) | Decrement `mRNA_counts[gene]` |
| 6 — Protein binding (per species) | Pick segment by per-segment on-rate weights; sample uniform position inside; if not blocked, append to `binding_proteins_positions`; split segment if topological barrier |
| 7 — Protein unbinding (per species) | Pick a bound molecule by per-molecule off-rate weights; remove from `binding_proteins_positions`; merge segments if topological barrier |

Each dispatch path may early-exit silently (e.g. recruitment blocked, TOP1 picked a plectonemic segment) — these still consume a Gillespie event.

See also: [Events and propensities](../user-guide/events-and-propensities.md), [Running simulations](../user-guide/simulation.md).
