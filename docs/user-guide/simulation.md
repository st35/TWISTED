# Running simulations

This page documents `SimulationSetupAndState`, the object that controls **how** the simulator runs and **collects all results**, together with the three callbacks accepted by `simulate_dynamics`.

```python
from model_setup import SimulationSetupAndState
from simulate_dynamics import simulate_dynamics

sim = SimulationSetupAndState(simulation_end_mode=0, simulation_end_criterion=300.0)
simulate_dynamics(model, sim)
```

---

## Constructor

```python
SimulationSetupAndState(
    simulation_end_mode: int,
    simulation_end_criterion: float | list[int],
    integration_method: str = 'RK23',
    integration_time_resolution: float = 0.1,
    integration_rtol: float = 1.0e-8,
    integration_atol: float = 1.0e-10,
    RNAP_alive_status_check_interval: float = 1.0,
    max_RNAPs_to_recruit: list[int] | None = None,
    Gillespie_random_seed: int = 42,
    everything_else_random_seed: int = 42,
)
```

| Argument | Meaning |
|----------|--------|
| `simulation_end_mode` | `0` for time-based, `1` for event-count based |
| `simulation_end_criterion` | `float` (seconds) if mode = 0; `list[int]` per gene if mode = 1 |
| `integration_method` | One of `'RK23'`, `'RK45'`, `'DOP853'`, `'Radau'`, `'BDF'`, `'LSODA'` |
| `integration_time_resolution` | spacing of `t_eval` points within an integration window (s) |
| `integration_rtol` | relative tolerance passed to `solve_ivp` (default `1e-8`) |
| `integration_atol` | absolute tolerance passed to `solve_ivp` (default `1e-10`) |
| `RNAP_alive_status_check_interval` | how often the integrator pauses to remove finished RNAPs and check for events (s) |
| `max_RNAPs_to_recruit` | optional list capping recruitment per gene |
| `Gillespie_random_seed` | integer seed for the Gillespie event-time and event-selection RNG (default `42`) |
| `everything_else_random_seed` | integer seed for all other stochastic choices: segment selection, binding positions, etc. (default `42`) |

The constructor does **not** take a `GenomicSetup`. Per-gene result lists are sized, and length checks against the gene list are performed, by `setup_simulation_state(genomic_setup)`, which `simulate_dynamics` calls automatically before the main loop. You normally do not need to call it yourself.

---

## Termination modes

### Time-based (`simulation_end_mode = 0`)

```python
sim = SimulationSetupAndState(
    simulation_end_mode=0,
    simulation_end_criterion=300.0,   # seconds
)
```

The loop exits as soon as `sim.curr_simulation_time >= simulation_end_criterion`.

### Event-count based (`simulation_end_mode = 1`)

```python
sim = SimulationSetupAndState(
    simulation_end_mode=1,
    simulation_end_criterion=[100, 50],   # one entry per gene
)
```

The loop exits as soon as **every** gene has produced at least its target number of completed transcripts. `sim.curr_simulation_time` ends at whatever value satisfies that condition.

The list length must equal the number of genes; this is checked at the start of `simulate_dynamics` (by `setup_simulation_state`).

### Limiting recruitment

```python
sim = SimulationSetupAndState(
    simulation_end_mode=1,
    simulation_end_criterion=[50],
    max_RNAPs_to_recruit=[50],
)
```

Once `len(sim.RNAP_recruitment_times[i]) >= max_RNAPs_to_recruit[i]`, further recruitment events on gene `i` silently fail (the Gillespie event still fires, but no RNAP is added). When event-count mode is selected, `max_RNAPs_to_recruit[i]` must be `>=` `simulation_end_event_counts[i]`; otherwise the run could never complete, and `setup_simulation_state` raises before the main loop starts.

---

## Integrator settings

| Setting | Default | Effect |
|---------|--------|-------|
| `integration_method` | `'RK23'` | The SciPy `solve_ivp` method to use. Other methods emit a `UserWarning` because higher-order solvers can produce non-physical intermediate states (e.g. RNAPs overlapping) and crash the run. |
| `integration_time_resolution` | `0.1` s | Spacing of `t_eval` points within each integration window. Affects only the temporal resolution at which `print_at_each_integration_step` observes the state; the dynamics are unchanged. |
| `RNAP_alive_status_check_interval` | `1.0` s | The integrator advances in chunks of this size, then removes RNAPs that finished and re-evaluates whether the cumulative propensity has crossed the event threshold. Smaller values are safer (fewer mis-detections of fast events) but slower. |

---

## Reading results

### Counts and timing

| Attribute | Type | Meaning |
|----------|------|--------|
| `sim.RNAPs_finished_transcription` | `list[int]` | Number of RNAPs that completed transcription, per gene |
| `sim.RNAP_recruitment_times[i]` | `list[float]` | Times (s) at which each RNAP on gene `i` was recruited |
| `sim.RNAP_exit_times[i]` | `list[float]` | Times (s) at which each completed RNAP on gene `i` exited |
| `sim.RNAPs_exit_positions[i]` | `list[float]` | Positions (nm) at which each completed RNAP on gene `i` exited |
| `sim.curr_simulation_time` | `float` | Total simulated seconds |
| `sim.simulation_completed` | `bool` | `True` iff the loop exited because the criterion was met |

For an alive RNAP at the end of a time-based run, only `RNAP_recruitment_times[i]` will contain its entry; no exit time or exit position is recorded.

### Transcription rates

```python
rates = sim.calculate_RNAP_transcription_rates(model)
```

Returns `list[list[float]]`: for each gene, the per-RNAP mean transcription rate in **bp/s**. Computed only for RNAPs that completed transcription:

```
rate = (exit_position − TSS) × sign / 0.34 / (exit_time − recruitment_time)        # bp/s
```

This method asserts that the distance covered is at least `gene_length`, so a non-empty list always corresponds to genuinely completed transcripts.

### Inter-recruitment intervals

```python
import numpy as np
intervals = np.diff(sim.RNAP_recruitment_times[0])
mean_interval = float(intervals.mean())
```

---

## Callbacks

`simulate_dynamics` accepts three optional callables:

```python
simulate_dynamics(
    model,
    sim,
    print_at_each_integration_step=None,   # (model, sim, t, state_vector)
    print_at_each_simulation_step=None,    # (model, sim)
    print_at_end_of_simulation=None,       # (model, sim)
)
```

| Callback | Signature | Fires |
|----------|----------|------|
| `print_at_each_integration_step` | `(model, sim, t, state_vector)` | Once at the start of every integration window of length `RNAP_alive_status_check_interval`. The `state_vector` is the live integrator state (RNAP positions, angles, segment `Lk`, cumulative propensity). |
| `print_at_each_simulation_step` | `(model, sim)` | Once at the top of each iteration of the outer Gillespie loop, **before** event selection. |
| `print_at_end_of_simulation` | `(model, sim)` | Once when the loop exits. |

The integration callback is suitable for logging time series of RNAP positions and supercoiling densities. The Gillespie callback is suitable for logging per-event quantities (mRNA counts, recruitment-success rates). The end callback is suitable for summary output.

### Worked logging snippet

```python
from typing import TextIO

def log_step(x_file: TextIO, sigma_file: TextIO,
             model, sim, t: float, state_vector: list[float]):
    RNAP_gene_index, _ = get_state_vectors_from_dicts(model)
    n = len(RNAP_gene_index)
    seg_lengths, seg_sigmas, *_ = calculate_segments_attributes(model, RNAP_gene_index, state_vector)

    x_file.write(str(t))
    for i in range(n):
        x_file.write('\t' + str(state_vector[i]))
    x_file.write('\n')

    sigma_file.write(str(t))
    for s in seg_sigmas:
        sigma_file.write('\t' + str(s))
    sigma_file.write('\n')

with open('x.log', 'w') as xf, open('sigma.log', 'w') as sf:
    simulate_dynamics(
        model, sim,
        print_at_each_integration_step=lambda m, s, t, sv: log_step(xf, sf, m, s, t, sv),
    )
```

The number of segments changes whenever an RNAP is recruited or exits, so the per-row width of `sigma.log` varies. Either pad on read, or post-process by joining `t` to the segment boundaries reconstructed from `x.log`.

---

## Complete example

```python
from model_setup import GenomicSetup, ModelSetup, Model, SimulationSetupAndState
from simulate_dynamics import simulate_dynamics
import numpy as np

genomic_setup = GenomicSetup(
    chromatin_type='prokaryotic',
    gene_names=['geneA', 'geneB'],
    TSSes=[340.0, 4420.0],
    gene_lengths=[3400.0, 3400.0],
    gene_directions=[1, 1],
    RNAP_on_rates=[0.02, 0.02],
    promoter_mode='constitutive',
    buffer_length=4420.0,
)
model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.01,
    TOP2_effective_relaxation_rate=0.005,
)
model = Model(genomic_setup, model_setup)
sim = SimulationSetupAndState(
    simulation_end_mode=1,
    simulation_end_criterion=[100, 100],
    max_RNAPs_to_recruit=[100, 100],
)

simulate_dynamics(model, sim)

for i, name in enumerate(genomic_setup.gene_names):
    rates = sim.calculate_RNAP_transcription_rates(model)[i]
    print(f'{name}: mean transcription rate = {np.mean(rates):.1f} bp/s '
          f'over {len(rates)} completed transcripts')
print(f'Total simulated time: {sim.curr_simulation_time:.1f} s')
```
