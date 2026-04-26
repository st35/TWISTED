# Overview

This page describes the **architecture** of TWISTED: the contents of the central objects, the structure of the simulation loop, and the layout of the state vector. For scripting alone, the [Getting Started](../getting-started.md) guide is sufficient; this page is intended as a reference for extending the simulator, debugging, or writing non-trivial callbacks.

---

## The four user-facing objects

```
GenomicSetup    →   gene geometry, recruitment rates, chromatin type
ModelSetup      →   physical/biological constants and the relaxation mode
Model           →   bundles both, plus all dynamic state mutated during the run
SimulationSetupAndState   →   integrator settings, termination, results
```

`Model` and `SimulationSetupAndState` are *mutated in place* by `simulate_dynamics`. After the run, all results reside on one of them.

| Object | Defined in | Contents |
|--------|-----------|----------|
| `GenomicSetup` | [`model_setup.py`](../api/model-setup.md#class-genomicsetup) | gene names, TSSes, lengths, directions, recruitment rates, promoter mode, chromatin type, nucleosome parameters, `clamp_left`, `clamp_right` |
| `ModelSetup` | [`model_setup.py`](../api/model-setup.md#class-modelsetup) | DNA/torque parameters, RNAP velocity parameters, topoisomerase parameters, steric parameters, boundary conditions, relaxation mode, mRNA mode |
| `Model` | [`model_setup.py`](../api/model-setup.md#class-model) | `x_dict`, `theta_dict`, `Lk`, `promoter_status`, `mRNA_counts`, `binding_proteins`, `binding_proteins_positions` |
| `SimulationSetupAndState` | [`model_setup.py`](../api/model-setup.md#class-simulationsetupandstate) | termination criterion, integrator settings, `RNAPs_finished_transcription`, recruitment/exit times, exit positions, `curr_simulation_time` |

A separate `BindingProtein` class (also in `model_setup.py`) describes any DNA-binding protein species other than RNAP and is passed to `Model` at construction.

---

## Module map

```
utilities.py
   │
   ▼
model_setup.py        (data classes; imports utilities)
   │
   ▼
biol_methods.py       (physics: torque, velocity, dLk/dt; imports model_setup)
   │
   ▼
model_dynamics.py     (state ↔ dict, ODE RHS, event rates; imports the above)
   │
   ▼
simulate_dynamics.py  (main loop and event dispatch)
```

There is one circular import at the top: `utilities.py` imports `model_setup`, and `model_setup.py` imports `utilities`. This resolves at runtime because the symbols actually *used* by `utilities` are looked up lazily inside function bodies.

---

## The simulation loop

`simulate_dynamics` is a single `while True` loop. Each iteration performs:

1. **Optional callback.** If `print_at_each_simulation_step` was provided, it is called with `(model, sim)`.
2. **Build the state vector** from the model's dictionaries via [`get_state_vectors_from_dicts`](../api/model-dynamics.md#get_state_vectors_from_dicts). This also returns `RNAP_gene_index`, an integer per RNAP recording which gene it belongs to.
3. **Sample an event waiting threshold.** Draw `p₀ ~ U(0, 1)`; the next discrete event happens when the cumulative propensity since the start of integration exceeds `ln(1/p₀)`.
4. **Integrate** via [`integrate`](../api/model-dynamics.md#integrate). The ODE [`model_dynamics`](../api/model-dynamics.md#model_dynamics) advances `x`, `θ`, `Lk`, and a final cumulative-propensity variable simultaneously. Integration runs in chunks of `RNAP_alive_status_check_interval` seconds; after each chunk, RNAPs that finished are removed and their adjacent segments merged.
5. **Push the new state back** into the model dictionaries via [`update_dicts_from_state_vector`](../api/model-dynamics.md#update_dicts_from_state_vector). Increment `sim.curr_simulation_time` by the integration `dt`.
6. **Recompute segment attributes** via [`calculate_segments_attributes`](../api/model-dynamics.md#calculate_segments_attributes): segment lengths, `σ`, torques, DNA states, writhe fractions, and plectoneme thresholds, at the new state.
7. **Recompute all event rates** via [`get_events_rates`](../api/model-dynamics.md#get_events_rates) and select an event index by drawing `p₁ ~ U(0, 1)` and walking the propensities.
8. **Execute the chosen event.** This is one of: RNAP recruitment, observation (no-op), global SC relaxation, type-specific SC relaxation, TOP1/TOP2 event, mRNA degradation, binding-protein on event, binding-protein off event. See [Events and propensities](events-and-propensities.md) for the per-event semantics.
9. **Check termination** (time-based or event-count based). If satisfied, mark `sim.simulation_completed = True` and break.
10. **Optional end callback.** If `print_at_end_of_simulation` was provided, call it.

The integration callback `print_at_each_integration_step` is fired *inside* step 4, at the start of every chunk.

---

## State vector layout

The state vector passed to the ODE solver and exposed to the integration callback is a flat list:

```
state_vector = [ x₁, x₂, …, xₙ,            # n RNAP positions, in nm
                 θ₁, θ₂, …, θₙ,            # n RNAP angular positions, in rad
                 Lk₀, Lk₁, …, Lkₘ,         # m+1 segment linking numbers, in turns
                 A ]                        # cumulative propensity since t_start
```

where:

- `n` is the current number of *alive* RNAPs;
- the segments are indexed **right to left**, so `Lk₀` is the rightmost segment (between `clamp_right` and the rightmost RNAP/topological barrier);
- the number of segments `m + 1` equals (RNAP count) + (topological-barrier protein count) + 1; for a stretch of DNA with no RNAPs and no topological barriers there is a single segment;
- the trailing element `A` is the time-integrated total propensity used for Gillespie waiting-time sampling (see step 4 above).

Inside the model, RNAPs are stored per gene in `model.x_dict[gene_index]` and `model.theta_dict[gene_index]`. Conversion between the two representations is handled by [`get_state_vectors_from_dicts`](../api/model-dynamics.md#get_state_vectors_from_dicts) and [`update_dicts_from_state_vector`](../api/model-dynamics.md#update_dicts_from_state_vector). The list ordering inside each `model.x_dict[i]` follows the gene direction: for `+1`-strand genes, positions are stored in left-to-right order; for `-1`-strand genes, the lists are reversed when round-tripping through the state vector.

---

## What lives on `Model` after the run

| Attribute | Meaning |
|----------|--------|
| `model.x_dict[i]` | list of RNAP positions on gene `i` (nm) |
| `model.theta_dict[i]` | list of RNAP angular positions on gene `i` (rad) |
| `model.Lk` | list of segment linking numbers (turns), right-to-left |
| `model.promoter_status[i]` | 1 if promoter ON, 0 if OFF (always 1 for `'constitutive'`) |
| `model.mRNA_counts[i]` | mRNA count on gene `i` |
| `model.binding_proteins[k]` | `k`-th `BindingProtein` species (nucleosomes are at index 0 in eukaryotic mode) |
| `model.binding_proteins_positions[k]` | bound positions of species `k` (nm) |

**These attributes must not be modified directly during a simulation.** Use the helpers in `model_dynamics.py` (e.g. `update_Lk_vector_after_*`), which keep the segment count consistent with the protein/RNAP layout.

---

## What lives on `SimulationSetupAndState` after the run

| Attribute | Meaning |
|----------|--------|
| `sim.RNAPs_finished_transcription[i]` | count of completed transcripts on gene `i` |
| `sim.RNAP_recruitment_times[i]` | list of times (s) at which each RNAP on gene `i` was recruited |
| `sim.RNAP_exit_times[i]` | list of times (s) at which each RNAP on gene `i` exited |
| `sim.RNAPs_exit_positions[i]` | list of positions (nm) at which each RNAP exited |
| `sim.curr_simulation_time` | total simulated seconds |
| `sim.simulation_completed` | `True` if the loop exited because the criterion was met |

The convenience method `sim.calculate_RNAP_transcription_rates(model)` returns, for each gene, the list of mean transcription rates (bp/s) of all RNAPs that completed transcription on that gene.

---

## Further reading

- [Genomic setup](genomic-setup.md): `GenomicSetup` reference.
- [Model parameters](model-setup.md): `ModelSetup` parameter set.
- [Binding proteins](binding-proteins.md): `BindingProtein` semantics.
- [Relaxation modes](relaxation-modes.md): behaviour of each supercoiling-relaxation mode.
- [Running simulations](simulation.md): `SimulationSetupAndState` and callbacks.
- [Events and propensities](events-and-propensities.md): exhaustive list of event types in the Gillespie loop.
- [Theory → DNA mechanics](../theory/dna-mechanics.md): equations underlying the ODE.
