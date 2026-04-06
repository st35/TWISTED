# Running Simulations

This page explains how to assemble the simulation objects, configure termination conditions, extract results, and add optional progress callbacks.

---

## Object Assembly

A complete simulation requires four objects:

```
GenomicSetup  +  ModelSetup
        │               │
        └───────┬───────┘
                ▼
             Model
                │
                ▼
    SimulationSetupAndState
                │
                ▼
       simulate_dynamics()
```

---

## `Model`

`Model` bundles a `GenomicSetup` and a `ModelSetup` and initialises all dynamic state variables.

```python
model = Model(genomic_setup, model_setup)
```

After construction, `model` holds:

| Attribute | Description |
|-----------|-------------|
| `model.x_dict` | `list[list[float]]` — RNAP positions per gene (nm) |
| `model.theta_dict` | `list[list[float]]` — RNAP angular positions per gene (rad) |
| `model.Lk` | `list[float]` — linking numbers of all DNA segments |
| `model.promoter_status` | `list[int]` — `1` (ON) or `0` (OFF) per gene |
| `model.mRNA_counts` | `list[int]` — mRNA copy numbers per gene |
| `model.binding_proteins` | `list[BindingProtein]` — binding protein types (auto-includes nucleosomes for eukaryotic) |
| `model.binding_proteins_positions` | `list[list[float]]` — positions (nm) of bound proteins per type |
| `model.topoisomerase_*` | Topoisomerase state (only when `mode == 'topoisomerase_based'`) |

**Do not modify these attributes directly** during a simulation; use the provided functions in `model_dynamics.py`.

---

## `SimulationSetupAndState`

```python
SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode,
    simulation_end_criterion,
    integration_method='RK23',
    integration_time_resolution=0.1,
    RNAP_alive_status_check_interval=1.0,
    max_RNAPs_to_recruit=None,
)
```

### Termination Modes

#### Time-based (`simulation_end_mode=0`)

```python
sim = SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode=0,
    simulation_end_criterion=1000.0,   # stop after 1000 s
)
```

#### Event-based (`simulation_end_mode=1`)

```python
sim = SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode=1,
    simulation_end_criterion=[100, 50],  # per-gene RNAP completion counts
)
```

The simulation stops when **all** genes have reached their respective target.

### Limiting RNAP Recruitment

```python
sim = SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode=1,
    simulation_end_criterion=[50],
    max_RNAPs_to_recruit=[50],   # recruit at most 50 RNAPs
)
```

`max_RNAPs_to_recruit[i]` must be ≥ `simulation_end_criterion[i]` in event-based mode.

### Integration Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `integration_time_resolution` | 0.1 s | Time step used for `t_eval` points within each ODE integration interval |
| `RNAP_alive_status_check_interval` | 1.0 s | Interval at which RNAP completion status is checked during integration |

Decrease `integration_time_resolution` for higher temporal resolution in output or to detect fast events more accurately. Smaller values increase computation time.

---

## Running the Simulation

```python
from simulate_dynamics import simulate_dynamics

simulate_dynamics(
    model,
    simulation_setup_and_state,
    print_at_each_integration_step=None,
    print_at_each_simulation_step=None,
    print_at_end_of_simulation=None,
)
```

All three `print_*` arguments accept an optional callable with signature:

```python
def my_callback(model: Model, sim: SimulationSetupAndState) -> None:
    ...
```

(The `print_at_each_integration_step` callback also receives `t: float` and `state_vector: list[float]` arguments.)

### Example: Progress Logging

```python
def log_progress(model, sim):
    print(f't={sim.curr_simulation_time:.1f}s  '
          f'finished={sim.RNAPs_finished_transcription}')

simulate_dynamics(model, sim, print_at_each_simulation_step=log_progress)
```

---

## Extracting Results

After the simulation completes, results are stored in `SimulationSetupAndState`:

| Attribute | Type | Description |
|-----------|------|-------------|
| `sim.RNAPs_finished_transcription` | `list[int]` | Number of RNAPs that completed transcription per gene |
| `sim.RNAP_recruitment_times` | `list[list[float]]` | Recruitment time of each RNAP per gene (s) |
| `sim.RNAP_exit_times` | `list[list[float]]` | Transcription completion time per RNAP per gene (s) |
| `sim.RNAPs_exit_positions` | `list[list[float]]` | Position where each RNAP exited the gene (nm) |
| `sim.curr_simulation_time` | `float` | Final simulation time (s) |
| `sim.simulation_completed` | `bool` | `True` if termination criterion was met |

### Transcription Rates

```python
# Returns list[list[float]] — transcription rate in bp/s for each RNAP per gene
rates = sim.calculate_RNAP_transcription_rates(model)
```

### Inter-RNAP Intervals

```python
import numpy as np

# Recruitment intervals (s) for gene 0
intervals = np.diff(sim.RNAP_recruitment_times[0])
mean_interval = np.mean(intervals)
```

---

## Complete Example with Two Genes

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
    supercoiling_relaxation_dynamics_mode='global_overall',
    global_supercoiling_relaxation_rate=0.1,
)

model = Model(genomic_setup, model_setup)
sim = SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode=1,
    simulation_end_criterion=[100, 100],
    max_RNAPs_to_recruit=[100, 100],
)

simulate_dynamics(model, sim)

for i, name in enumerate(genomic_setup.gene_names):
    rates = sim.calculate_RNAP_transcription_rates(model)[i]
    print(f'{name}: mean transcription rate = {np.mean(rates):.1f} bp/s')
```
