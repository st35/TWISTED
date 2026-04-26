# Getting Started

This page covers the path from a fresh checkout to a first running simulation. The reader will, by the end, understand:

1. how to install TWISTED;
2. the **four objects** every simulation requires and their relationships;
3. the meaning of each parameter in the minimal example;
4. how to retrieve the results.

Readers already familiar with the construction pattern may skip directly to [Tutorials](tutorials.md).

---

## 1. Install

TWISTED is a small pure-Python project with two runtime dependencies. Clone and install them:

```bash
git clone https://github.com/st35/TWISTED.git
cd TWISTED/code
pip install scipy numpy
```

There is no `setup.py`; the modules reside directly in `code/` and are imported by name. Either run Python from inside `code/`, or add it to `sys.path`:

```python
import sys
sys.path.append('/absolute/path/to/TWISTED/code')
```

To build and serve the documentation locally:

```bash
pip install -r docs/requirements.txt
mkdocs serve
```

| Package | Purpose |
|---------|---------|
| `scipy` | `solve_ivp` (RK23 by default) for the ODE component of the simulation |
| `numpy` | array bookkeeping inside the integrator wrapper |

---

## 2. Conceptual model

Consider a one-dimensional DNA molecule running from `clamp_left = 0` on the left to `clamp_right` (in nm) on the right. RNA polymerases reside on it. The molecule is partitioned by the RNAPs (and optionally by *topological-barrier* proteins) into a sequence of **segments**. Each segment carries its own linking number `Lk`, from which a supercoiling density `σ` and torque `τ` are derived.

The simulation alternates between two regimes:

- **Continuous** — an ODE solver advances RNAP positions `xᵢ`, RNAP angular positions `θᵢ`, and segment linking numbers `Lkⱼ` smoothly through time. The torque imbalance across an RNAP slows it and rotates it; rotation injects or removes twist.
- **Discrete** — at random times sampled from a Gillespie clock, a single event fires: a new RNAP is recruited at a TSS, an mRNA is degraded, a topoisomerase relaxes a segment, a binding protein binds or unbinds, and so on.

The integrator and Gillespie machinery are not invoked directly. Instead, four objects are assembled and passed to `simulate_dynamics`:

```
GenomicSetup    →   gene layout (geometry + recruitment rates)
ModelSetup      →   physical/biological constants and selected relaxation mode
Model           →   bundles the two and holds dynamic state
SimulationSetupAndState   →   termination criterion and integrator settings
                ↓
        simulate_dynamics(model, sim)
```

The `Model` and `SimulationSetupAndState` objects are *mutated in place* during the run. After `simulate_dynamics` returns, results are read from them.

---

## 3. Minimal example, line by line

```python
from model_setup import GenomicSetup, ModelSetup, Model, SimulationSetupAndState
from simulate_dynamics import simulate_dynamics
```

Three modules cover everything required at the top level.

### Step 1 — describe the DNA and the genes

```python
genomic_setup = GenomicSetup(
    chromatin_type='prokaryotic',
    gene_names=['lacZ'],
    TSSes=[340.0],
    gene_lengths=[3400.0],
    gene_directions=[1],
    RNAP_on_rates=[0.02],
    promoter_mode='constitutive',
    buffer_length=3400.0,
)
```

| Argument | Meaning |
|----------|---------|
| `chromatin_type` | `'prokaryotic'` selects the five-state Marko torque law; `'eukaryotic'` selects the chromatin torque law and adds nucleosomes (see [Tutorial 15](tutorials.md#15-eukaryotic-chromatin-with-nucleosomes)) |
| `gene_names` | identifier strings, one per gene |
| `TSSes` | transcription start sites in **nm** (multiply bp by 0.34) |
| `gene_lengths` | distance an RNAP must travel from its TSS to exit |
| `gene_directions` | `+1` (left → right) or `-1` (right → left) |
| `RNAP_on_rates` | recruitment rate when the promoter is ON, in s⁻¹ |
| `promoter_mode` | `'constitutive'` (always ON) or `'non-constitutive'` (see [User Guide → Not yet implemented](user-guide/not-yet-implemented.md)) |
| `buffer_length` | extra DNA past the gene, in nm; sets `clamp_right` |

The right end of the DNA is computed automatically as

```
clamp_right = TSSes[0] + gene_lengths[0] + buffer_length    # if gene 0 is on the + strand
clamp_right = TSSes[0] + buffer_length                      # if gene 0 is on the − strand
```

so when simulating multiple genes, **list the rightmost one first** or increase `buffer_length` accordingly. See [Tutorial 4](tutorials.md#4-multi-gene-tandem) for details.

### Step 2 — choose physical parameters

```python
model_setup = ModelSetup(
    v0=20.0,
    tau_c=12.0,
    supercoiling_relaxation_dynamics_mode='global_overall',
    global_supercoiling_relaxation_rate=0.1,
)
```

The defaults of `ModelSetup` are tuned for *E. coli*-like prokaryotic transcription. Two are overridden here:

- `v0` — maximum RNAP elongation velocity (nm/s). 20 nm/s ≈ 59 bp/s.
- `tau_c` — torque sensitivity in the velocity law `v = (v₀/2)(1 − tanh((τ_f − τ_b)/τ_c))`.

The third pair selects the simplest supercoiling-relaxation mode (`global_overall`): a single Poisson event at rate `0.1 / s` resets every segment's `Lk` to its relaxed value `Lk₀`. The five implemented modes are catalogued in [User Guide → Relaxation modes](user-guide/relaxation-modes.md) and compared in [Tutorial 7](tutorials.md#7-comparing-relaxation-modes).

### Step 3 — bundle into a model

```python
model = Model(genomic_setup, model_setup)
```

`Model` initialises all dynamic state: empty RNAP lists per gene, a single DNA segment of relaxed `Lk`, all promoters set to ON (for `'constitutive'`), zero mRNA counts, and any nucleosomes (for `'eukaryotic'`).

### Step 4 — say how to integrate and when to stop

```python
sim = SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode=0,             # 0 = time-based; 1 = event-count based
    simulation_end_criterion=500.0,    # 500 seconds
)
```

`simulation_end_mode=0` runs for a fixed number of seconds. `simulation_end_mode=1` stops when each gene has produced at least N completed transcripts (specified by a list `[N1, N2, ...]`). See [Tutorial 2](tutorials.md#2-event-based-termination).

Three optional settings control the integrator:

- `integration_method='RK23'` (default). Other SciPy solvers are accepted but trigger a warning, since they may produce non-physical intermediate states that violate steric constraints.
- `integration_time_resolution=0.1` — spacing of `t_eval` points within an integration interval (seconds).
- `RNAP_alive_status_check_interval=1.0` — every interval of this length, the integrator pauses to check which RNAPs have finished and to allow a discrete event to fire.

### Step 5 — run

```python
simulate_dynamics(model, sim)
```

Control returns when `sim.simulation_completed` becomes `True`.

### Step 6 — read the results

```python
print('Per-gene completed transcripts:', sim.RNAPs_finished_transcription)
print('Recruitment times (gene 0):',     sim.RNAP_recruitment_times[0])
print('Exit times (gene 0):',            sim.RNAP_exit_times[0])
print('Final mRNA counts:',              model.mRNA_counts)

rates = sim.calculate_RNAP_transcription_rates(model)
mean_rate = sum(rates[0]) / len(rates[0]) if rates[0] else None
print('Mean transcription rate (bp/s) on gene 0:', mean_rate)
```

The full list of attributes available on `sim` and `model` after a run is documented in [User Guide → Running simulations](user-guide/simulation.md).

---

## 4. Common variants

### Stop after N completed transcripts instead of after a time

```python
sim = SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode=1,
    simulation_end_criterion=[50],   # one entry per gene
)
```

### Cap the total number of recruitment events

```python
sim = SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode=1,
    simulation_end_criterion=[50],
    max_RNAPs_to_recruit=[50],       # never recruit more than 50 on this gene
)
```

In event-count mode `max_RNAPs_to_recruit[i]` must be ≥ `simulation_end_criterion[i]`.

### Two genes

```python
genomic_setup = GenomicSetup(
    chromatin_type='prokaryotic',
    gene_names=['geneA', 'geneB'],
    TSSes=[340.0, 4080.0],
    gene_lengths=[3400.0, 3400.0],
    gene_directions=[1, 1],
    RNAP_on_rates=[0.02, 0.02],
    promoter_mode='constitutive',
    buffer_length=4080.0,            # large enough to cover the rightmost gene
)
```

When multiple genes share DNA, RNAP–RNAP steric exclusion is enforced automatically through a soft tanh ramp; see [Theory](theory/dna-mechanics.md#4-rnap-velocity).

---

## 5. Further reading

- [Tutorials](tutorials.md) — sixteen worked examples, one per feature.
- [User Guide → Overview](user-guide/overview.md) — architecture and the simulation loop in detail.
- [User Guide → Model parameters](user-guide/model-setup.md) — the full `ModelSetup` parameter set.
- [User Guide → Relaxation modes](user-guide/relaxation-modes.md) — selection among the five implemented modes.
- [Theory → DNA mechanics](theory/dna-mechanics.md) — equations for torque, velocity, and `dLk/dt`.
