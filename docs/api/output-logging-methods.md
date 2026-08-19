# `output_logging_methods`

Ready-made callback functions for the three hooks accepted by [`simulate_dynamics`](simulate-dynamics.md): `print_at_each_integration_step`, `print_at_each_simulation_step` (Gillespie-step logging), and `print_at_end_of_simulation`. Each function writes tab-delimited rows to caller-supplied file handles, so nothing is written unless you pass a file for it.

---

## `print_at_each_integration_step`

```python
print_at_each_integration_step(
    model: Model,
    simulation_setup_and_state: SimulationSetupAndState,
    t: float,
    state_vector: list[float],
    RNAPs_pos_file: TextIO = None,
    RNAPs_dx_dt_file: TextIO = None,
    RNAPs_dtheta_dt_file: TextIO = None,
    sigma_file: TextIO = None,
    torque_file: TextIO = None,
    dLk_dt_file: TextIO = None,
    nucleosome_position_file: TextIO = None,
    promoter_status_file: TextIO = None,
    mRNA_count_file: TextIO = None,
    barrier_position_file: TextIO = None,
    binding_proteins_file: TextIO = None,
    screen_logging: bool = False,
) -> None
```

Matches the `(model, sim, t, state_vector)` signature of the `print_at_each_integration_step` hook. Pass it directly, or wrap it in a `lambda` to bind the file handles:

```python
simulate_dynamics(
    model, sim,
    print_at_each_integration_step=lambda m, s, t, sv: print_at_each_integration_step(
        m, s, t, sv, RNAPs_pos_file=pos_f, sigma_file=sigma_f,
    ),
)
```

Each row is `t` followed by one tab-separated value per RNAP or per segment (segments are ordered right to left, as everywhere else in the package). Only the arguments you supply a file for are written:

| Argument | Row content |
|----------|-------------|
| `RNAPs_pos_file` | RNAP positions (nm), one column per RNAP |
| `RNAPs_dx_dt_file` | RNAP linear velocities (nm/s) |
| `RNAPs_dtheta_dt_file` | RNAP angular velocities (rad/s) |
| `sigma_file` | Per-segment supercoiling density σ |
| `torque_file` | Per-segment torque (pN·nm) |
| `dLk_dt_file` | Per-segment linking number `Lk` (see note below) |
| `nucleosome_position_file` | Nucleosome positions (nm); only written for `chromatin_type == 'eukaryotic'` |
| `barrier_position_file` | Chromosome-barrier positions (nm); only written when `genomic_setup.are_multiple_chromosomes_present` |
| `binding_proteins_file` | Positions of all bound generic binding proteins, excluding the nucleosome and chromosome-barrier entries (which have their own files) |
| `promoter_status_file` | `model.promoter_status[i]` for every gene |
| `mRNA_count_file` | Per gene: number of RNAPs currently transcribing, `RNAPs_finished_transcription[i]`, then `model.mRNA_counts[i]` |

If `screen_logging` is `True`, a summary line (`t`, RNAP count, finished-transcription counts, mRNA counts, promoter statuses) is also printed to stdout.

!!! warning "`dLk_dt_file` naming"
    Despite the name, this argument currently writes the segment's linking number `Lk` (taken from `state_vector`), not its time derivative `dLk/dt` (which is available in `model_dynamics`'s return value). Confirm this is the intended behavior before relying on it for rate calculations.

---

## `print_at_each_Gillespie_step`

```python
print_at_each_Gillespie_step(
    model: Model,
    simulation_setup_and_state: SimulationSetupAndState,
    events_log_file: TextIO = None,
    log_nucleosomal_events: bool = False,
) -> None
```

Matches the `(model, sim)` signature of the `print_at_each_simulation_step` hook. Writes one row per Gillespie iteration to `events_log_file`: `curr_simulation_time`, `last_event_index` (or `'NA'` before the first event), and `last_event_type` (or `'NA'`). Rows whose `last_event_type` contains `'nucleosome'` are skipped unless `log_nucleosomal_events=True`, since nucleosome binding/unbinding events are typically far more frequent than other event types and can dominate the log.

---

## `print_at_end_of_simulation`

```python
print_at_end_of_simulation(
    model: Model,
    simulation_setup_and_state: SimulationSetupAndState,
    transcription_rates_file: TextIO = None,
    final_mRNA_counts_file: TextIO = None,
) -> None
```

Matches the `(model, sim)` signature of the `print_at_end_of_simulation` hook. Calls `simulation_setup_and_state.calculate_RNAP_transcription_rates(model)` and, if `transcription_rates_file` is given, writes one row per gene (`gene_name` followed by the per-RNAP transcription rates in bp/s). If `final_mRNA_counts_file` is given, writes a single row with `model.mRNA_counts` for all genes.
