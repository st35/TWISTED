# Tutorials

This page provides step-by-step tutorials for running simulations with TWISTED.
Each example builds on the previous one, progressing from a minimal single-gene
simulation to more advanced multi-gene setups with output logging and different
relaxation modes.

---

## Prerequisites

Before starting, make sure the TWISTED `code/` directory is on your Python path:

```python
import sys
sys.path.append('/path/to/TWISTED/code')

from utilities import *
from model_setup import *
from simulate_dynamics import *
```

All tutorials also assume you have **NumPy** and **SciPy** installed (see
[Getting Started](getting-started.md)).

---

## Example 1 — Single Gene, Minimal Setup

The simplest possible simulation: one constitutive gene with global supercoiling
relaxation, run until a fixed number of RNAPs finish transcription.

### 1.1 Create a gene configuration file

TWISTED reads gene geometry from a **tab-delimited config file** where each line
specifies one gene with five columns:

| Column | Description | Units |
|--------|-------------|-------|
| 1 | Gene name | — |
| 2 | TSS position | bp |
| 3 | Gene length | bp |
| 4 | Direction (`+1` or `-1`) | — |
| 5 | RNAP on-rate | 1/s |

Create a file called `single_gene.config`:

```text
geneA	10000	5300	1	1.0
```

This places a 5300 bp gene on the positive strand with its TSS at 10 000 bp.
The RNAP on-rate is 1.0 s⁻¹ (this can be scaled via the `explicit_RNAP_on_rates` argument to `construct_genomic_setup`).

### 1.2 Build the genomic setup

```python
chromatin_type = 'prokaryotic'
promoter_mode = 'constitutive'

genomic_setup = construct_genomic_setup(
    'single_gene.config',
    chromatin_type,
    promoter_mode,
)
genomic_setup.print_genomic_setup()
```

`construct_genomic_setup` reads the config file, converts positions from bp to
nm internally (1 bp = 0.34 nm), and returns a `GenomicSetup` object. Calling
`print_genomic_setup()` prints a summary table of the loaded genes.

### 1.3 Configure model parameters

```python
model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='global_overall',
    global_supercoiling_relaxation_rate=0.008,   # 1/s
)
```

`'global_overall'` is the simplest relaxation mode: supercoiling density relaxes
uniformly across the entire DNA at the given rate. The remaining physical
parameters (`w0`, `chi`, `eta`, `v0`, `tau_c`, etc.) use their defaults,
which are appropriate for a standard prokaryotic system.

### 1.4 Create the model and simulation state

```python
model = Model(genomic_setup, model_setup)

# Event-based termination: stop after 50 RNAPs finish transcription on geneA
simulation_setup_and_state = SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode=1,        # 1 = event-based
    simulation_end_criterion=[50], # 50 completed transcriptions for geneA
)
```

### 1.5 Run the simulation

```python
simulate_dynamics(model, simulation_setup_and_state)
```

After the simulation completes, you can inspect the results:

```python
# Transcription rates (bp/s) for each RNAP that completed transcription
transcription_rates = simulation_setup_and_state.calculate_RNAP_transcription_rates(model)

avg_rate = (
    sum(transcription_rates[0]) / len(transcription_rates[0])
    if transcription_rates[0] else 'N/A'
)
print(f'Average transcription rate for geneA: {avg_rate} bp/s')
print(f'Total RNAPs finished: {simulation_setup_and_state.RNAPs_finished_transcription[0]}')
```

---

## Example 2 — Convergent Gene Pair with Supercoiling Logging

Two convergent genes on opposite strands, with callback functions that log
supercoiling density to files at each integration step.

### 2.1 Config file

Create `convergent_pair.config`:

```text
gene_R	25300	5300	-1	1.0
gene_F	10000	5300	1	1.0
```

`gene_F` is on the positive strand (direction = `+1`), `gene_R` is on the
reverse strand (direction = `-1`). Their 3′ ends face each other, making this a
**convergent** pair.

!!! note
    The first gene listed in the config file determines the right boundary of
    the DNA (`clamp_right`). List the gene that defines the rightmost extent
    first, or increase `buffer_length` to cover all genes.

### 2.2 Set up the model

```python
import random
random.seed(42)  # For reproducibility

from typing import TextIO

chromatin_type = 'prokaryotic'
promoter_mode = 'constitutive'

genomic_setup = construct_genomic_setup(
    'convergent_pair.config',
    chromatin_type,
    promoter_mode,
    explicit_RNAP_on_rates=[0.083, 0.083],  # ~5 RNAPs/min each
)
genomic_setup.print_genomic_setup()

model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.004,  # 1/s
    TOP2_effective_relaxation_rate=0.004,  # 1/s
)

model = Model(genomic_setup, model_setup)

# Event-based: stop after 20 RNAPs finish on each gene
simulation_setup_and_state = SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode=1,
    simulation_end_criterion=[20, 20],
)
```

### 2.3 Define logging callbacks

TWISTED provides three callback hooks in `simulate_dynamics`:

| Callback | Signature | Called when |
|----------|-----------|------------|
| `print_at_each_integration_step` | `(model, simulation_setup_and_state, t, state_vector)` | Start of each integration interval (approximately every `RNAP_alive_status_check_interval` seconds) |
| `print_at_each_simulation_step` | `(model, simulation_setup_and_state)` | Every Gillespie step (before event selection) |
| `print_at_end_of_simulation` | `(model, simulation_setup_and_state)` | Simulation ends |

Here is a callback that logs RNAP positions and supercoiling densities:

!!! note
    The `state_vector` passed to the callback reflects the current ODE solver
    state. We extract `RNAP_gene_index` from the model to identify which gene
    each RNAP belongs to, but use the passed `state_vector` for positions and
    segment calculations.

```python
def log_integration_step(
    x_file: TextIO,
    sigma_file: TextIO,
    model: Model,
    simulation_setup_and_state: SimulationSetupAndState,
    t: float,
    state_vector: list[float],
) -> None:
    RNAP_gene_index, _ = get_state_vectors_from_dicts(model)
    RNAP_count = len(RNAP_gene_index)
    segments_lengths, segments_sigmas, segments_torques, _, _, _ = (
        calculate_segments_attributes(model, RNAP_gene_index, state_vector)
    )

    # Write RNAP positions
    x_file.write(str(t))
    for i in range(RNAP_count):
        x_file.write('\t' + str(state_vector[i]))
    x_file.write('\n')

    # Write supercoiling densities per segment
    sigma_file.write(str(t))
    for sigma in segments_sigmas:
        sigma_file.write('\t' + str(sigma))
    sigma_file.write('\n')
```

### 2.4 Run with logging

```python
with (
    open('rnap_positions.log', 'w') as x_file,
    open('sigma_values.log', 'w') as sigma_file,
):
    simulate_dynamics(
        model,
        simulation_setup_and_state,
        print_at_each_integration_step=lambda model, ss, t, sv: (
            log_integration_step(x_file, sigma_file, model, ss, t, sv)
        ),
    )
```

After the simulation, `rnap_positions.log` contains tab-separated RNAP
positions at each time point and `sigma_values.log` contains supercoiling density
per segment. These can be loaded for plotting:

```python
import numpy as np

sigma_data = np.loadtxt('sigma_values.log')
time = sigma_data[:, 0]
sigmas = sigma_data[:, 1:]
```

---

## Example 3 — Three Genes with mRNA Degradation

A setup with three genes, mRNA degradation dynamics, and the
`topoisomerase_approximated` relaxation mode.

### 3.1 Config file

Create `three_genes.config`:

```text
geneC	35300	5300	-1	1.0
geneB	20000	5300	1	1.0
geneA	10000	5300	1	1.0
```

### 3.2 Set up and run

```python
import random
random.seed(2)

chromatin_type = 'prokaryotic'
promoter_mode = 'constitutive'

genomic_setup = construct_genomic_setup(
    'three_genes.config',
    chromatin_type,
    promoter_mode,
    explicit_RNAP_on_rates=[0.083, 0.083, 0.083],
)

model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.008,
    TOP2_effective_relaxation_rate=0.008,
    mRNA_dynamics_mode=1,             # Enable mRNA degradation
    mRNA_degradation_rate=0.05,       # 1/s
)

model = Model(genomic_setup, model_setup)

# Event-based: stop after 50 transcriptions completed per gene
simulation_setup_and_state = SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode=1,
    simulation_end_criterion=[50, 50, 50],
)
```

To track mRNA counts during the simulation, add a Gillespie-step callback:

```python
def log_gillespie_step(model, simulation_setup_and_state):
    t = simulation_setup_and_state.curr_simulation_time
    counts = model.mRNA_counts
    print(f't={t:.2f}s  mRNA counts: {counts}')

simulate_dynamics(
    model,
    simulation_setup_and_state,
    print_at_each_simulation_step=log_gillespie_step,
)
```

After completion, inspect per-gene results:

```python
transcription_rates = simulation_setup_and_state.calculate_RNAP_transcription_rates(model)
for i, name in enumerate(model.genomic_setup.gene_names):
    rates = transcription_rates[i]
    avg = sum(rates) / len(rates) if rates else 0
    print(f'{name}: {len(rates)} RNAPs completed, avg rate = {avg:.1f} bp/s')
print(f'Final mRNA counts: {model.mRNA_counts}')
```

---

## Example 4 — Time-Based Termination with Free Boundary Conditions

Simulate for a fixed duration (180 seconds) with free (unclamped) DNA ends
instead of the default clamped boundaries.

### 4.1 Set up

```python
import random
random.seed(1)

genomic_setup = construct_genomic_setup(
    'single_gene.config',   # Reuse the config from Example 1
    'prokaryotic',
    'constitutive',
)

model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.008,
    TOP2_effective_relaxation_rate=0.008,
    clamps_status=('free', 'free'),   # Both DNA ends are free
)

model = Model(genomic_setup, model_setup)

# Time-based termination: run for 180 seconds
simulation_setup_and_state = SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode=0,       # 0 = time-based
    simulation_end_criterion=180.0,  # seconds
)
```

### 4.2 Run and inspect

```python
simulate_dynamics(model, simulation_setup_and_state)

print(f'Simulation ended at t = {simulation_setup_and_state.curr_simulation_time:.2f} s')
print(f'RNAPs completed: {simulation_setup_and_state.RNAPs_finished_transcription}')
```

With free boundary conditions, supercoiling can dissipate at the DNA ends, so
you will typically observe higher transcription rates compared to the clamped
case.

---

## Example 5 — Explicit Topoisomerase Dynamics

The `topoisomerase_based` mode explicitly models individual topoisomerase
molecules binding and unbinding to DNA segments, with steric interactions between
topoisomerases and RNAPs.

### 5.1 Set up

```python
import random
random.seed(10)

genomic_setup = construct_genomic_setup(
    'single_gene.config',
    'prokaryotic',
    'constitutive',
    explicit_RNAP_on_rates=[0.083],
)

# Topoisomerase parameters
TOP1_count = 5
TOP2_count = 5
TOP1_on_rate = 0.083    # 1/s
TOP1_off_rate = 0.0008  # 1/s
TOP2_on_rate = 0.083    # 1/s
TOP2_off_rate = 0.0008  # 1/s

model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_based',
    topoisomerase_copy_numbers=[TOP1_count, TOP2_count],
    topoisomerase_on_off_rates=[
        (TOP1_on_rate, TOP1_off_rate),
        (TOP2_on_rate, TOP2_off_rate),
    ],
)

model = Model(genomic_setup, model_setup)

simulation_setup_and_state = SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode=1,
    simulation_end_criterion=[30],
)
```

### 5.2 Run

```python
simulate_dynamics(model, simulation_setup_and_state)

print(f'RNAPs completed: {simulation_setup_and_state.RNAPs_finished_transcription[0]}')
```

In `topoisomerase_based` mode, each topoisomerase molecule is tracked
individually. TOP1 relaxes torsional stress on non-plectonemic DNA (writhe
fraction = 0), while TOP2 relaxes plectonemic supercoiling (writhe fraction
> 0). Both enzymes act on positive and negative supercoiling, driving $\sigma$
toward zero. Their binding positions are subject to steric hindrance from
nearby RNAPs (exclusion distance `(RNAP_diameter + TOPO_diameter) / 2`, default
15 nm).

---

## Example 6 — Comparing Relaxation Modes

This example demonstrates how to sweep across different supercoiling relaxation
modes and compare transcription rates.

```python
import random

results = {}

# --- Mode 1: global_overall ---
random.seed(7)
genomic_setup = construct_genomic_setup('single_gene.config', 'prokaryotic', 'constitutive')
model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='global_overall',
    global_supercoiling_relaxation_rate=0.008,
)
model = Model(genomic_setup, model_setup)
ss = SimulationSetupAndState(genomic_setup, 1, [30])
simulate_dynamics(model, ss)
rates = ss.calculate_RNAP_transcription_rates(model)[0]
results['global_overall'] = sum(rates) / len(rates) if rates else 0

# --- Mode 2: global_by_type ---
random.seed(7)
genomic_setup = construct_genomic_setup('single_gene.config', 'prokaryotic', 'constitutive')
model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='global_by_type',
    local_supercoiling_relaxation_rates=[0.004, 0.004],  # [positive, negative]
)
model = Model(genomic_setup, model_setup)
ss = SimulationSetupAndState(genomic_setup, 1, [30])
simulate_dynamics(model, ss)
rates = ss.calculate_RNAP_transcription_rates(model)[0]
results['global_by_type'] = sum(rates) / len(rates) if rates else 0

# --- Mode 3: topoisomerase_approximated ---
random.seed(7)
genomic_setup = construct_genomic_setup('single_gene.config', 'prokaryotic', 'constitutive')
model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.004,
    TOP2_effective_relaxation_rate=0.004,
)
model = Model(genomic_setup, model_setup)
ss = SimulationSetupAndState(genomic_setup, 1, [30])
simulate_dynamics(model, ss)
rates = ss.calculate_RNAP_transcription_rates(model)[0]
results['topoisomerase_approximated'] = sum(rates) / len(rates) if rates else 0

# --- Print comparison ---
print(f"{'Mode':<30} {'Avg transcription rate (bp/s)':>30}")
print('-' * 62)
for mode, avg_rate in results.items():
    print(f'{mode:<30} {avg_rate:>30.2f}')
```

---

## Example 7 — Full Logging to Files

A complete example with logging of RNAP positions, supercoiling densities, and
mRNA counts to separate files, combined with time-based termination.

### 7.1 Set up

```python
import sys
sys.path.append('/path/to/TWISTED/code')

import random
random.seed(1)

from typing import TextIO

from utilities import *
from model_setup import *
from simulate_dynamics import *

chromatin_type = 'prokaryotic'
promoter_mode = 'constitutive'

genomic_setup = construct_genomic_setup(
    'convergent_pair.config',  # Two genes
    chromatin_type,
    promoter_mode,
)
genomic_setup.print_genomic_setup()

TOP1_rate = 8.33  # High topoisomerase activity (1/s)
TOP2_rate = 8.33

model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=TOP1_rate,
    TOP2_effective_relaxation_rate=TOP2_rate,
    mRNA_dynamics_mode=1,
    mRNA_degradation_rate=0.05,
)

model = Model(genomic_setup, model_setup)

simulation_setup_and_state = SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode=0,          # Time-based
    simulation_end_criterion=180.0, # 180 seconds
)
```

### 7.2 Define callbacks

```python
def log_integration(
    x_file: TextIO,
    sigma_file: TextIO,
    mRNA_file: TextIO,
    model: Model,
    simulation_setup_and_state: SimulationSetupAndState,
    t: float,
    state_vector: list[float],
) -> None:
    RNAP_gene_index, _ = get_state_vectors_from_dicts(model)
    RNAP_count = len(RNAP_gene_index)
    segments_lengths, segments_sigmas, _, _, _, _ = calculate_segments_attributes(
        model, RNAP_gene_index, state_vector
    )

    # RNAP positions
    x_file.write(str(t))
    for i in range(RNAP_count):
        x_file.write('\t' + str(state_vector[i]))
    x_file.write('\n')

    # Supercoiling densities
    sigma_file.write(str(t))
    for sigma in segments_sigmas:
        sigma_file.write('\t' + str(sigma))
    sigma_file.write('\n')

    # mRNA counts
    mRNA_file.write(str(t) + '\t')
    mRNA_file.write('\t'.join(str(c) for c in model.mRNA_counts))
    mRNA_file.write('\n')


def log_gillespie(model, simulation_setup_and_state):
    t = simulation_setup_and_state.curr_simulation_time
    rnap_counts = [len(model.x_dict[g]) for g in range(len(model.genomic_setup.gene_names))]
    finished = simulation_setup_and_state.RNAPs_finished_transcription
    print(f't={t:.2f}s  active RNAPs={rnap_counts}  finished={finished}')


def log_end(model, simulation_setup_and_state):
    rates = simulation_setup_and_state.calculate_RNAP_transcription_rates(model)
    for i, name in enumerate(model.genomic_setup.gene_names):
        avg = sum(rates[i]) / len(rates[i]) if rates[i] else 'N/A'
        print(f'{name}: avg transcription rate = {avg} bp/s')
```

### 7.3 Run

```python
with (
    open('x_values.log', 'w') as x_file,
    open('sigma_values.log', 'w') as sigma_file,
    open('mRNA_counts.log', 'w') as mRNA_file,
):
    simulate_dynamics(
        model,
        simulation_setup_and_state,
        print_at_each_integration_step=lambda model, ss, t, sv: (
            log_integration(x_file, sigma_file, mRNA_file, model, ss, t, sv)
        ),
        print_at_each_simulation_step=log_gillespie,
        print_at_end_of_simulation=log_end,
    )
```

The output files can then be loaded and plotted with any standard tool
(matplotlib, MATLAB, etc.).

---

## Example 8 — DNA-Binding Proteins

TWISTED can model additional DNA-binding proteins that stochastically bind and
unbind from DNA segments. Each protein type is defined as a `BindingProtein`
object and passed to the `Model` constructor.

### 8.1 Define a binding protein

```python
from model_setup import BindingProtein

# A hypothetical protein with 10 copies, binding at 0.001 / (s * nm)
# and unbinding at 0.01 / s
my_protein = BindingProtein(
    protein_name='ProteinX',
    total_copy_number=10,
    is_steric_barrier_to_RNAPs=False,
    is_topological_barrier=False,
    basal_on_rate=0.001,                 # s⁻¹ nm⁻¹
    basal_off_rate=0.01,                 # s⁻¹
)
```

The effective per-segment on-rate is `basal_on_rate × segment_length × n_unbound`,
so longer segments attract more binding. The off-rate is constant per bound
molecule.

### 8.2 Custom rate functions

You can supply optional `on_rate_func` and `off_rate_func` callables that
modulate the basal rates based on segment length and supercoiling density:

```python
# Protein that binds faster on negatively supercoiled DNA
def enhanced_binding(segment_length, segment_sigma):
    return 5.0 if segment_sigma < -0.02 else 1.0

my_protein = BindingProtein(
    protein_name='SC_sensor',
    total_copy_number=20,
    is_steric_barrier_to_RNAPs=False,
    is_topological_barrier=False,
    basal_on_rate=0.0005,
    basal_off_rate=0.02,
    on_rate_func=enhanced_binding,  # Multiplies basal_on_rate × segment_length
)
```

The effective on-rate becomes
`n_unbound × basal_on_rate × segment_length × on_rate_func(L, σ)`.

### 8.3 Run a simulation with binding proteins

```python
import random
random.seed(3)

genomic_setup = construct_genomic_setup(
    'single_gene.config',
    'prokaryotic',
    'constitutive',
)

model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.004,
    TOP2_effective_relaxation_rate=0.004,
)

# Pass binding proteins to the Model constructor
model = Model(genomic_setup, model_setup, binding_proteins=[my_protein])

simulation_setup_and_state = SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode=1,
    simulation_end_criterion=[30],
)

simulate_dynamics(model, simulation_setup_and_state)

# Inspect bound protein positions after the simulation
print(f'Bound {my_protein.protein_name} count: {len(model.binding_proteins_positions[0])}')
print(f'Positions: {model.binding_proteins_positions[0]}')
```

### 8.4 Multiple binding protein types

You can define several protein types with different kinetics and pass them all
to `Model`. Each type is tracked independently — positions for protein type `i`
are stored in `model.binding_proteins_positions[i]`.

```python
import random
random.seed(5)

# Fast-binding, fast-unbinding protein (high turnover)
transient_binder = BindingProtein(
    protein_name='Transient',
    total_copy_number=50,
    is_steric_barrier_to_RNAPs=False,
    is_topological_barrier=False,
    basal_on_rate=0.005,    # s⁻¹ nm⁻¹
    basal_off_rate=0.5,     # s⁻¹ (short residence time)
)

# Slow-binding, slow-unbinding protein (long-lived occupancy)
stable_binder = BindingProtein(
    protein_name='Stable',
    total_copy_number=5,
    is_steric_barrier_to_RNAPs=False,
    is_topological_barrier=False,
    basal_on_rate=0.0002,   # s⁻¹ nm⁻¹
    basal_off_rate=0.005,   # s⁻¹ (long residence time)
)

# Supercoiling-sensitive protein that unbinds faster under positive supercoiling
def sc_dependent_off(segment_length, segment_sigma):
    return 10.0 if segment_sigma > 0.02 else 1.0

sc_sensor = BindingProtein(
    protein_name='SC_sensor',
    total_copy_number=15,
    is_steric_barrier_to_RNAPs=False,
    is_topological_barrier=False,
    basal_on_rate=0.001,
    basal_off_rate=0.01,
    off_rate_func=sc_dependent_off,  # 10× faster unbinding on positively supercoiled DNA
)

# Set up and run
genomic_setup = construct_genomic_setup(
    'single_gene.config',
    'prokaryotic',
    'constitutive',
)

model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.004,
    TOP2_effective_relaxation_rate=0.004,
)

model = Model(
    genomic_setup,
    model_setup,
    binding_proteins=[transient_binder, stable_binder, sc_sensor],
)

simulation_setup_and_state = SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode=0,
    simulation_end_criterion=120.0,  # 120 seconds
)

simulate_dynamics(model, simulation_setup_and_state)

# Report bound counts for each protein type
for i, protein in enumerate(model.binding_proteins):
    n_bound = len(model.binding_proteins_positions[i])
    print(f'{protein.protein_name}: {n_bound} / {protein.total_copy_number} bound')
```

---

## Example 9 — Eukaryotic Simulation with Nucleosomes

TWISTED supports eukaryotic chromatin simulations where nucleosomes dynamically
bind and unbind from DNA, modulating the torque–supercoiling relationship through
a density-dependent buffering mechanism.

### 9.1 Config file

Create `eukaryotic_gene.config`:

```text
geneA	10000	5300	1	1.0
```

### 9.2 Set up a eukaryotic genomic setup

To use eukaryotic mode, set `chromatin_type='eukaryotic'` and optionally
provide nucleosome-related keyword arguments:

```python
import random
random.seed(42)

from utilities import *
from model_setup import *
from simulate_dynamics import *

genomic_setup = construct_genomic_setup(
    'eukaryotic_gene.config',
    chromatin_type='eukaryotic',
    per_nucleosome_DNA_length=147,       # bp (converted to nm internally)
    nucleosome_linker_length=30,         # bp (converted to nm internally)
    nucleosomes_are_steric_barriers_to_RNAPs=True,
)
genomic_setup.print_genomic_setup()
print(f'Total nucleosomes: {genomic_setup.get_total_nucleosome_count()}')
```

The nucleosome count is automatically computed by tiling the domain with
nucleosomes spaced by `per_nucleosome_DNA_length + nucleosome_linker_length`.
You can override this with an explicit count:

```python
# Or specify an explicit nucleosome count
genomic_setup_explicit = construct_genomic_setup(
    'eukaryotic_gene.config',
    chromatin_type='eukaryotic',
    nucleosome_count=50,  # Override automatic tiling
)
```

You can also provide custom rate functions that modulate nucleosome
binding and unbinding as a function of segment length and supercoiling density:

```python
# Custom nucleosome rate functions
def nucleosome_on_rate_modifier(segment_length, segment_sigma):
    # Suppress binding on highly negatively supercoiled segments
    return max(0.0, 1.0 + 10.0 * segment_sigma)

def nucleosome_off_rate_modifier(segment_length, segment_sigma):
    # Enhance unbinding on highly negatively supercoiled segments
    return max(1.0, 1.0 - 10.0 * segment_sigma)

genomic_setup_custom = construct_genomic_setup(
    'eukaryotic_gene.config',
    chromatin_type='eukaryotic',
    nucleosome_on_rate_func=nucleosome_on_rate_modifier,
    nucleosome_off_rate_func=nucleosome_off_rate_modifier,
)
```

These functions multiply the basal on/off rates (see
[BindingProtein](api/model-setup.md#bindingprotein) for details).

### 9.3 Configure model parameters

```python
model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.004,
    TOP2_effective_relaxation_rate=0.004,
)
```

### 9.4 Create the model

When you create a `Model` with a eukaryotic `GenomicSetup`, nucleosomes are
automatically added as the first entry in `model.binding_proteins`:

```python
model = Model(genomic_setup, model_setup)

# The nucleosome BindingProtein is at index 0
nucleosome_protein = model.binding_proteins[0]
print(f'Nucleosome protein: {nucleosome_protein.protein_name}')
print(f'Total copies: {nucleosome_protein.total_copy_number}')
print(f'Is nucleosome: {nucleosome_protein.is_a_nucleosome}')
print(f'Steric barrier to RNAPs: {nucleosome_protein.is_steric_barrier_to_RNAPs}')
```

!!! note
    Any additional binding proteins you pass to `Model(... binding_proteins=[...])` will
    appear at indices 1, 2, … after the auto-created nucleosome entry.

### 9.5 Run the simulation

```python
simulation_setup_and_state = SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode=1,
    simulation_end_criterion=[20],
)

simulate_dynamics(model, simulation_setup_and_state)

# Inspect results
transcription_rates = simulation_setup_and_state.calculate_RNAP_transcription_rates(model)
avg_rate = (
    sum(transcription_rates[0]) / len(transcription_rates[0])
    if transcription_rates[0] else 'N/A'
)
print(f'Average transcription rate: {avg_rate} bp/s')
print(f'Bound nucleosomes at end: {len(model.binding_proteins_positions[0])}')
```

The eukaryotic torque model features a **buffering regime** where nucleosomes
absorb positive supercoiling without increasing torque. The width of this
buffering regime scales with the local nucleosome density $\psi$ (fraction of
segment length occupied by nucleosomes). As RNAPs displace nucleosomes, $\psi$
decreases and the buffering capacity shrinks, allowing supercoiling to build up.

---

## Quick Reference

### Config file format

Tab-delimited, one gene per line:

```
name    TSS_bp    length_bp    direction    RNAP_on_rate
```

### Relaxation modes

| Mode | Required keyword arguments |
|------|---------------------------|
| `global_overall` | `global_supercoiling_relaxation_rate` |
| `global_per_segment` | `global_supercoiling_relaxation_rate` |
| `global_by_type` | `local_supercoiling_relaxation_rates` (list of 2 floats) |
| `per_segment_by_type` | `local_supercoiling_relaxation_rates` (list of 2 floats) |
| `topoisomerase_approximated` | `TOP1_effective_relaxation_rate`, `TOP2_effective_relaxation_rate` |
| `topoisomerase_based` | `topoisomerase_copy_numbers`, `topoisomerase_on_off_rates` |

### Simulation termination

| `simulation_end_mode` | `simulation_end_criterion` | Description |
|-----------------------|---------------------------|-------------|
| `0` | `float` (seconds) | Stop after a fixed time |
| `1` | `list[int]` (one per gene) | Stop after N RNAPs finish per gene |

### Callback signatures

```python
# Called at the start of each integration interval (~every RNAP_alive_status_check_interval seconds)
def on_integration_step(model, simulation_setup_and_state, t, state_vector): ...

# Called at every Gillespie step (before event selection)
def on_gillespie_step(model, simulation_setup_and_state): ...

# Called once when the simulation ends
def on_end(model, simulation_setup_and_state): ...
```
