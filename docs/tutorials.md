# Tutorials

Twenty worked examples, ordered roughly by complexity. Each section is self-contained and can be copy-pasted into a Python script. Earlier examples introduce idioms (config files, callbacks) that later ones reuse.

All examples assume the `code/` directory is importable:

```python
import sys
sys.path.append('/path/to/TWISTED/code')

from utilities import *
from model_setup import *
from simulate_dynamics import *
```

For reproducibility, pass explicit seeds to `SimulationSetupAndState`:

```python
sim = SimulationSetupAndState(
    ...,
    Gillespie_random_seed=42,
    everything_else_random_seed=42,
)
```

---

## Table of contents

1. [Single gene, minimal setup](#1-single-gene-minimal-setup)
2. [Event-based termination](#2-event-based-termination)
3. [Loading genes from a config file](#3-loading-genes-from-a-config-file)
4. [Multi-gene tandem](#4-multi-gene-tandem)
5. [Convergent gene pair](#5-convergent-gene-pair)
6. [Divergent gene pair](#6-divergent-gene-pair)
7. [Comparing relaxation modes](#7-comparing-relaxation-modes)
8. [Topoisomerase-approximated mode in depth](#8-topoisomerase-approximated-mode-in-depth)
9. [mRNA dynamics](#9-mrna-dynamics)
10. [Boundary conditions: clamped vs free](#10-boundary-conditions-clamped-vs-free)
11. [Callbacks and logging](#11-callbacks-and-logging)
12. [Generic DNA-binding proteins](#12-generic-dna-binding-proteins)
13. [Displaceable proteins at the TSS](#13-displaceable-proteins-at-the-tss)
14. [Topological-barrier proteins](#14-topological-barrier-proteins)
15. [Eukaryotic chromatin with nucleosomes](#15-eukaryotic-chromatin-with-nucleosomes)
16. [Reproducibility and integrator choice](#16-reproducibility-and-integrator-choice)
17. [Non-constitutive promoters](#17-non-constitutive-promoters)
18. [Multi-chromosome setup](#18-multi-chromosome-setup)
19. [Gene regulatory network: toggle switch](#19-gene-regulatory-network-toggle-switch)
20. [Saving and resuming a simulation](#20-saving-and-resuming-a-simulation)

---

## 1. Single gene, minimal setup

A single 10 kb gene transcribed for 300 simulated seconds with the simplest relaxation model.

```python
genomic_setup = GenomicSetup(
    chromatin_type='prokaryotic',
    gene_names=['geneA'],
    TSSes=[340.0],
    gene_lengths=[3400.0],
    gene_directions=[1],
    RNAP_on_rates=[0.02],
    promoter_mode='constitutive',
    buffer_length=3400.0,
    are_multiple_chromosomes_present=False,
    regulatory_network_information=read_regulatory_network(['geneA'], None),
)

model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='global_overall',
    global_supercoiling_relaxation_rate=0.1,
)

model = Model(genomic_setup, model_setup)
sim = SimulationSetupAndState(
    simulation_end_mode=0,
    simulation_end_criterion=300.0,
    Gillespie_random_seed=0,
    everything_else_random_seed=0,
)
simulate_dynamics(model, sim)

print('Completed transcripts:', sim.RNAPs_finished_transcription)
```

`'global_overall'` resets *every* DNA segment's `Lk` to its relaxed value `Lk₀` whenever a single Poisson event fires (rate `global_supercoiling_relaxation_rate`). It is the simplest mechanism for preventing supercoiling from accumulating without bound.

---

## 2. Event-based termination

Stop the simulation when each gene has produced a target number of completed transcripts:

```python
sim = SimulationSetupAndState(
    simulation_end_mode=1,
    simulation_end_criterion=[50],   # one entry per gene
)
```

In this mode `sim.curr_simulation_time` records the time at which the last gene crosses the threshold. To additionally limit the *total* number of recruitment attempts, supply `max_RNAPs_to_recruit`:

```python
sim = SimulationSetupAndState(
    simulation_end_mode=1,
    simulation_end_criterion=[50],
    max_RNAPs_to_recruit=[50],       # never recruit more than 50 on gene 0
)
```

The constructor enforces `max_RNAPs_to_recruit[i] >= simulation_end_criterion[i]` so that the run can in principle complete.

---

## 3. Loading genes from a config file

Specifying `TSSes`, `gene_lengths`, and so on by hand becomes cumbersome for more than a few genes. The helper `construct_genomic_setup` reads a tab-delimited file:

```text
# columns:  name   TSS_bp   length_bp   direction   RNAP_on_rate(1/s)
geneA       10000  5300     1           1.0
```

Note that **positions and lengths in the file are in bp** — the loader multiplies by `0.34` to convert to nm.

```python
genomic_setup = construct_genomic_setup(
    'single_gene.config',
    chromatin_type='prokaryotic',
    promoter_mode='constitutive',
)
genomic_setup.print_genomic_setup()
```

The optional `explicit_RNAP_on_rates` keyword multiplies the per-gene rate read from the file:

```python
genomic_setup = construct_genomic_setup(
    'single_gene.config',
    chromatin_type='prokaryotic',
    promoter_mode='constitutive',
    explicit_RNAP_on_rates=[0.05],   # final rate = file_rate × 0.05
)
```

This is useful for parameter sweeps in which the relative ranking between genes is fixed in the file but the absolute scale is varied.

---

## 4. Multi-gene tandem

Two co-directional genes on the same DNA. Inter-RNAP steric exclusion is enforced automatically by a soft tanh ramp on each RNAP's velocity (see [Theory → RNAP velocity](theory/dna-mechanics.md#4-rnap-velocity)).

```python
genomic_setup = GenomicSetup(
    chromatin_type='prokaryotic',
    gene_names=['geneA', 'geneB'],
    TSSes=[340.0, 4420.0],          # B starts ~12 kb after A
    gene_lengths=[3400.0, 3400.0],
    gene_directions=[1, 1],
    RNAP_on_rates=[0.02, 0.02],
    promoter_mode='constitutive',
    buffer_length=4420.0,           # extends past the rightmost gene end
    are_multiple_chromosomes_present=False,
    regulatory_network_information=read_regulatory_network(['geneA', 'geneB'], None),
)

model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='global_overall',
    global_supercoiling_relaxation_rate=0.1,
)

model = Model(genomic_setup, model_setup)
sim = SimulationSetupAndState(1, [40, 40], Gillespie_random_seed=1, everything_else_random_seed=1)
simulate_dynamics(model, sim)
```

!!! warning "`clamp_right` is computed from gene 0 only"
    `clamp_right = TSSes[0] + gene_lengths[0] + buffer_length` (or `TSSes[0] + buffer_length` if gene 0 is on the − strand). When simulating multiple genes, either list the rightmost one first or select a `buffer_length` large enough to cover all genes.

---

## 5. Convergent gene pair

Two genes on opposite strands whose 3′ ends face each other. Positive supercoiling accumulates between them as both RNAPs push twist into the shared region.

```text
# convergent_pair.config
gene_R      25300   5300    -1   1.0
gene_F      10000   5300     1   1.0
```

```python
genomic_setup = construct_genomic_setup(
    'convergent_pair.config',
    chromatin_type='prokaryotic',
    promoter_mode='constitutive',
    explicit_RNAP_on_rates=[0.083, 0.083],   # ≈ 5 RNAPs/min each
)

model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.004,
    TOP2_effective_relaxation_rate=0.004,
)

model = Model(genomic_setup, model_setup)
sim = SimulationSetupAndState(1, [20, 20], Gillespie_random_seed=42, everything_else_random_seed=42)
simulate_dynamics(model, sim)
```

`gene_R` is listed first so that `clamp_right` is large enough to cover it. The `topoisomerase_approximated` mode is more realistic here than `global_overall` because TOP2 — which can resolve plectonemes — is the only enzyme that can rescue the highly positively supercoiled segment between the two converging RNAPs.

---

## 6. Divergent gene pair

Two genes on opposite strands whose 5′ ends face each other. Negative supercoiling builds up between them.

```text
# divergent_pair.config
gene_R      15300   5300    -1   1.0
gene_F      20000   5300     1   1.0
```

```python
genomic_setup = construct_genomic_setup(
    'divergent_pair.config',
    chromatin_type='prokaryotic',
    promoter_mode='constitutive',
    explicit_RNAP_on_rates=[0.05, 0.05],
)

model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.01,
    TOP2_effective_relaxation_rate=0.001,    # TOP2 less important when σ is negative
)
model = Model(genomic_setup, model_setup)
sim = SimulationSetupAndState(0, 300.0, Gillespie_random_seed=7, everything_else_random_seed=7)
simulate_dynamics(model, sim)
```

Compare the steady-state supercoiling profile from this run with the convergent case in Tutorial 5 by enabling logging (Tutorial 11).

---

## 7. Comparing relaxation modes

The five implemented modes are summarised in [User Guide → Relaxation modes](user-guide/relaxation-modes.md). The example below sweeps through all of them on the same single-gene setup.

```python
results = {}
configs = [
    ('global_overall',
        {'global_supercoiling_relaxation_rate': 0.008}),
    ('global_per_segment',
        {'global_supercoiling_relaxation_rate': 0.008}),
    ('global_by_type',
        {'local_supercoiling_relaxation_rates': [0.004, 0.004]}),
    ('per_segment_by_type',
        {'local_supercoiling_relaxation_rates': [0.004, 0.004]}),
    ('topoisomerase_approximated',
        {'TOP1_effective_relaxation_rate': 0.004,
         'TOP2_effective_relaxation_rate': 0.004}),
]

for mode, kw in configs:
    genomic_setup = construct_genomic_setup('single_gene.config', 'prokaryotic', 'constitutive')
    model_setup = ModelSetup(supercoiling_relaxation_dynamics_mode=mode, **kw)
    model = Model(genomic_setup, model_setup)
    sim = SimulationSetupAndState(1, [30], Gillespie_random_seed=0, everything_else_random_seed=0)
    simulate_dynamics(model, sim)
    rates = sim.calculate_RNAP_transcription_rates(model)[0]
    results[mode] = sum(rates) / len(rates) if rates else 0.0

print(f'{"Mode":<30}{"Mean transcription rate (bp/s)":>32}')
for mode, r in results.items():
    print(f'{mode:<30}{r:>32.2f}')
```

Down the list, the model becomes progressively more spatially resolved and biologically motivated, and progressively slower per simulation step. Exploratory work should generally start at the top.

---

## 8. Topoisomerase-approximated mode in depth

The `topoisomerase_approximated` mode introduces *two* independent Poisson event types with rates `TOP1_effective_relaxation_rate` and `TOP2_effective_relaxation_rate`:

- **TOP1 event** — pick one segment with probability proportional to its length. If the segment is negatively supercoiled (`σ < 0`), set its `Lk` to `Lk₀`. If it is positively supercoiled (`σ ≥ 0`) **with writhe** (`writhe_fraction > 0`), partially relax it by reducing `Lk` by `σ_s × Lk₀` (where `σ_s` is the plectoneme-formation threshold); if it is positively supercoiled **without writhe**, set its `Lk` to `Lk₀`.
- **TOP2 event** — pick one segment with probability proportional to its length. If the segment **has writhe** (`writhe_fraction > 0`), set its `Lk` to `Lk₀ × (1 + σ_s)`, where `σ_s` is the plectoneme-formation threshold (TOP2 reduces writhe but does not over-relax). Otherwise the event has no effect.

TOP1 therefore relieves *twist* (and only partially relaxes plectonemic positive supercoiling), TOP2 relieves *plectonemes*, and the supercoiling and writhe state of each segment determines how each enzyme acts on it. Reproducing in-vivo behaviour usually requires both rates to be non-zero. A typical exploratory pair is

```python
ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.01,
    TOP2_effective_relaxation_rate=0.005,
)
```

The fully explicit `'topoisomerase_based'` mode (binding/unbinding of individual enzyme molecules with steric interactions) is **not yet implemented**; see [User Guide → Not yet implemented](user-guide/not-yet-implemented.md).

---

## 9. mRNA dynamics

By default, mRNA counts are monotonically non-decreasing: every time an RNAP completes transcription, `model.mRNA_counts[gene_index]` is incremented. To add first-order degradation, enable `mRNA_dynamics_mode`:

```python
model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.008,
    TOP2_effective_relaxation_rate=0.008,
    mRNA_dynamics_mode=1,
    mRNA_degradation_rate=0.05,         # 1/s, per molecule
)
```

The Gillespie loop now includes a per-gene degradation event with rate `mRNA_degradation_rate × mRNA_counts[i]`. To monitor the count over time:

```python
def log_step(model, sim):
    print(f't={sim.curr_simulation_time:7.2f}  mRNA={model.mRNA_counts}')

simulate_dynamics(model, sim, print_at_each_simulation_step=log_step)
```

---

## 10. Boundary conditions: clamped vs free

By default, both ends of the DNA are torsionally clamped, so no twist can escape. Setting one or both to `'free'` allows twist to dissipate at the boundary; the corresponding boundary segment then has its `Lk` driven by RNAP *displacement* rather than by angular-velocity transfer (see [Theory](theory/dna-mechanics.md#6-linking-number-dynamics)).

```python
model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.008,
    TOP2_effective_relaxation_rate=0.008,
    clamps_status=('free', 'free'),
)
```

Transcription rates are typically higher under free boundary conditions because torsional stress is no longer accumulated against a hard wall.

---

## 11. Callbacks and logging

`simulate_dynamics` accepts three optional callables:

| Argument | Signature | Fired when |
|----------|----------|------------|
| `print_at_each_integration_step` | `(model, sim, t, state_vector)` | once at the start of every ODE integration window of length `RNAP_alive_status_check_interval`, throttled by `min_interval_between_calls_to_print_at_each_integration_step` |
| `print_at_each_simulation_step` | `(model, sim)` | once at the top of each iteration of the main Gillespie loop |
| `print_at_end_of_simulation` | `(model, sim)` | once when the loop exits |

The integration callback is appropriate for logging time series of RNAP positions and supercoiling densities, since the `state_vector` it receives is exactly the integrator state (RNAP positions, RNAP angles, segment `Lk`, and the cumulative propensity).

```python
from typing import TextIO

def log_integration(x_file: TextIO, sigma_file: TextIO,
                    model, sim, t: float, state_vector: list):
    RNAP_gene_index, _ = get_state_vectors_from_dicts(model)
    n = len(RNAP_gene_index)
    seg_lengths, seg_sigmas, *_ = calculate_segments_attributes(model, RNAP_gene_index, state_vector)

    x_file.write(str(t) + '\t' + '\t'.join(str(state_vector[i]) for i in range(n)) + '\n')
    sigma_file.write(str(t) + '\t' + '\t'.join(str(s) for s in seg_sigmas) + '\n')

with open('x.log', 'w') as xf, open('sigma.log', 'w') as sf:
    simulate_dynamics(
        model, sim,
        print_at_each_integration_step=lambda m, s, t, sv: log_integration(xf, sf, m, s, t, sv),
    )
```

The number of segments changes whenever an RNAP is recruited or finishes. For a fixed-width file suitable for plotting, post-process the log rather than attempting to align columns at write time.

Every Gillespie iteration also sets `sim.last_event_type` to a string identifying what just happened (e.g. `'RNAP_recruitment_successful'`, `'mRNA_degradation'`, `'promoter_ON'`, `'<protein_name>_binding'`). This makes it straightforward to filter `print_at_each_simulation_step` callbacks by event kind without having to inspect `model` state manually:

```python
def log_on_recruitment(model, sim):
    if sim.last_event_type == 'RNAP_recruitment_successful':
        print(f't={sim.curr_simulation_time:.2f}  RNAP recruited  total={sum(sim.RNAPs_finished_transcription)}')

simulate_dynamics(model, sim, print_at_each_simulation_step=log_on_recruitment)
```

For more detailed per-event introspection (which event fired, propensities, etc.), see [User Guide → Events and propensities](user-guide/events-and-propensities.md).

---

## 12. Generic DNA-binding proteins

Beyond RNAPs and (in eukaryotic mode) nucleosomes, additional generic DNA-binding protein species can be supplied via the `binding_proteins` argument of `Model`.

```python
my_protein = BindingProtein(
    protein_name='ProteinX',
    total_copy_number=10,
    is_steric_barrier_to_RNAPs=False,
    is_topological_barrier=False,
    basal_on_rate=0.001,    # per (s · nm)
    basal_off_rate=0.01,    # per s
)

model = Model(genomic_setup, model_setup, binding_proteins=[my_protein])
```

The basal on-rate is multiplied by segment length and number of unbound copies, so the **per-segment** binding propensity is

```
n_unbound × basal_on_rate × segment_length
```

and the off propensity for each bound molecule is `basal_off_rate`. After binding, the position within the chosen segment is sampled uniformly.

The on-rate may be modulated by a function of `(segment_length, segment_sigma)`, and the off-rate by a function of `(segment_length, segment_sigma, binding_position)` (the third argument is the bound molecule's position along the DNA in nm):

```python
def enhanced_on(segment_length, segment_sigma):
    return 5.0 if segment_sigma < -0.02 else 1.0

sc_sensor = BindingProtein(
    protein_name='SC_sensor',
    total_copy_number=20,
    is_steric_barrier_to_RNAPs=False,
    is_topological_barrier=False,
    basal_on_rate=0.0005,
    basal_off_rate=0.02,
    on_rate_func=enhanced_on,
)
```

The user function multiplies the basal rate; the wrapping by `× segment_length × n_unbound` (for binding) or by `n_bound` (for unbinding) is added by `BindingProtein` automatically. After the run, bound positions for protein `i` are stored in `model.binding_proteins_positions[i]`.

---

## 13. Displaceable proteins at the TSS

When an RNAP is recruited at a TSS, the simulator first checks whether anything is in the way (another RNAP, a nucleosome, or a generic protein). By default, any obstacle blocks recruitment.

Setting `can_be_displaced_at_TSS_by_RNAP=True` on a `BindingProtein` modifies this behaviour: if the *only* obstacle is a displaceable protein, recruitment **succeeds** and the blocking protein is removed from `model.binding_proteins_positions`.

```python
blocker = BindingProtein(
    protein_name='Displaceable',
    total_copy_number=10,
    is_steric_barrier_to_RNAPs=True,
    is_topological_barrier=False,
    basal_on_rate=0.001,
    basal_off_rate=0.01,
    can_be_displaced_at_TSS_by_RNAP=True,
)
```

This is appropriate for promoter-proximal repressors that can be evicted by initiating polymerase. The same flag exists for nucleosomes; see Tutorial 15.

---

## 14. Topological-barrier proteins

Some DNA-binding proteins are believed to act as **topological barriers**: they not only physically obstruct RNAPs but also pin the underlying DNA so that supercoiling cannot diffuse past them. In TWISTED, setting `is_topological_barrier=True` adds the protein to the list of segment boundaries — every bound molecule splits one DNA segment into two with independent linking numbers.

```python
boundary = BindingProtein(
    protein_name='Boundary',
    total_copy_number=4,
    is_steric_barrier_to_RNAPs=True,    # required: must also be a steric barrier
    is_topological_barrier=True,
    basal_on_rate=0.0005,
    basal_off_rate=0.001,               # long residence time
)

model = Model(genomic_setup, model_setup, binding_proteins=[boundary])
```

`is_topological_barrier=True` implies `is_steric_barrier_to_RNAPs=True`; the `BindingProtein` constructor raises an error otherwise. When a topological-barrier protein binds, the local supercoiling density is preserved on both new sub-segments. When it unbinds, the two adjacent segments are merged and their linking numbers are added.

---

## 15. Eukaryotic chromatin with nucleosomes

Switching `chromatin_type` to `'eukaryotic'` does two things:

1. The torque law switches to the chromatin-specific piecewise function (see [Theory → Eukaryotic torque](theory/dna-mechanics.md#3-eukaryotic-torque-law)). The "buffering" plateau in this law widens with nucleosome density `ψ`, modelling the absorption of positive supercoiling by chromatin.
2. A `BindingProtein` named `'nucleosome'` is **automatically added as `model.binding_proteins[0]`**. Its copy number is computed by tiling the DNA with `per_nucleosome_DNA_length + nucleosome_linker_length` spacing.

```python
genomic_setup = construct_genomic_setup(
    'single_gene.config',
    chromatin_type='eukaryotic',
    per_nucleosome_DNA_length=147,            # bp; converted to nm internally
    nucleosome_linker_length=30,              # bp; converted to nm internally
    nucleosomes_are_steric_barriers_to_RNAPs=True,
)
print('Nucleosomes tiled:', genomic_setup.get_total_nucleosome_count())

model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.004,
    TOP2_effective_relaxation_rate=0.004,
)
model = Model(genomic_setup, model_setup)
print('Auto-created nucleosome protein:', model.binding_proteins[0].protein_name)
```

The auto-tiled copy number can be overridden with `nucleosome_count=...`, the binding kinetics customised with `nucleosome_on_rate_func` / `nucleosome_off_rate_func`, and incoming RNAPs may be allowed to evict promoter-proximal nucleosomes with `nucleosomes_can_be_displaced_at_TSS_by_RNAP=True`. All eukaryotic-specific keyword arguments are listed in [User Guide → Genomic setup](user-guide/genomic-setup.md#eukaryotic-keyword-arguments).

Any additional `BindingProtein` objects passed to `Model` appear at indices `1, 2, …` after the auto-created nucleosome entry; index by name rather than by hard-coded position.

---

## 16. Reproducibility and integrator choice

TWISTED uses two independent seeded `random.Random` instances internally: one (`rng_Gillespie`) for event-time sampling and event selection, and one (`rng_everything_else`) for all other stochastic choices (segment selection, binding positions, etc.). Both are stored on the `SimulationSetupAndState` object and seeded from the values passed to its constructor, making each run fully reproducible (modulo SciPy's deterministic adaptive stepping):

```python
sim = SimulationSetupAndState(
    simulation_end_mode=0,
    simulation_end_criterion=300.0,
    Gillespie_random_seed=2026,
    everything_else_random_seed=2026,
)
```

Using separate seeds lets the Gillespie stream and the auxiliary stream be varied independently, which is useful for sensitivity analyses. Both default to `42`. Calling `random.seed()` globally has **no effect** on TWISTED's internal RNGs.

The integrator can be selected via `SimulationSetupAndState`:

```python
sim = SimulationSetupAndState(
    simulation_end_mode=0,
    simulation_end_criterion=300.0,
    integration_method='RK23',           # default; recommended
    integration_time_resolution=0.05,    # finer t_eval for higher-resolution logging
    integration_rtol=1.0e-6,             # relative tolerance passed to solve_ivp
    integration_atol=1.0e-8,             # absolute tolerance passed to solve_ivp
    RNAP_alive_status_check_interval=0.5,
    Gillespie_random_seed=2026,
    everything_else_random_seed=2026,
)
```

`'RK23'` is the default. Selecting any other method (`'RK45'`, `'DOP853'`, `'Radau'`, `'BDF'`, `'LSODA'`) emits a `UserWarning` because higher-order solvers can take steps that produce non-physical intermediate states (RNAPs overlapping, segments shorter than allowed by steric constraints) and crash the simulation. RK23 should be retained unless there is a specific reason to switch and the result has been validated.

`integration_time_resolution` controls only the spacing of `t_eval` points within an integration window; reducing it does **not** change the dynamics, only the temporal resolution at which the integration callback (Tutorial 11) observes the state. `integration_rtol` and `integration_atol` set the relative and absolute error tolerances passed to `solve_ivp` (defaults `1e-6` and `1e-8`); tightening them increases accuracy at the cost of more integration steps. `RNAP_alive_status_check_interval` does affect the dynamics indirectly: it bounds the duration the integrator runs before re-checking which RNAPs have finished and re-computing event rates. Smaller values are safer but slower.

!!! tip "Unstable dynamics"
    If a run shows numerical instability — rapid, non-physical fluctuations in RNAP positions or supercoiling densities — lower (tighten) `integration_rtol` and `integration_atol`, for example back to `1e-8` and `1e-10` or smaller. Tighter tolerances force `solve_ivp` to take smaller steps, which suppresses these fluctuations at the cost of longer runtimes.

---

## 17. Non-constitutive promoters

With `promoter_mode='non-constitutive'` each gene's promoter switches stochastically between OFF (`promoter_status[i] = 0`) and ON (`promoter_status[i] = 1`) as Gillespie events. RNAP recruitment on gene `i` occurs only while its promoter is ON. All promoters start in the OFF state at the beginning of the simulation.

The switching rates are supplied per gene via `TF_on_off_rates`, a list of `(k_on, k_off)` pairs in s⁻¹:

```python
genomic_setup = GenomicSetup(
    chromatin_type='prokaryotic',
    gene_names=['geneA'],
    TSSes=[340.0],
    gene_lengths=[3400.0],
    gene_directions=[1],
    RNAP_on_rates=[0.05],
    promoter_mode='non-constitutive',
    TF_on_off_rates=[(0.02, 0.01)],   # k_on=0.02 s⁻¹, k_off=0.01 s⁻¹
    buffer_length=3400.0,
    are_multiple_chromosomes_present=False,
    regulatory_network_information=read_regulatory_network(['geneA'], None),
)

model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.008,
    TOP2_effective_relaxation_rate=0.008,
)
model = Model(genomic_setup, model_setup)

sim = SimulationSetupAndState(
    simulation_end_mode=0,
    simulation_end_criterion=600.0,
    Gillespie_random_seed=7,
    everything_else_random_seed=7,
)

```

The mean fraction of time the promoter spends ON is `k_on / (k_on + k_off)`. With the values above that is `0.02 / 0.03 ≈ 0.67`, so the effective RNAP recruitment rate is roughly `0.67 × 0.05 ≈ 0.033 s⁻¹`, compared to `0.05 s⁻¹` for a constitutive promoter of the same basal rate.

The `print_at_each_simulation_step` callback is convenient for tracking promoter switching events. Use `sim.last_event_type` to fire the callback only on the relevant transitions:

```python
def log_promoter(model, sim):
    if sim.last_event_type in ('promoter_ON', 'promoter_OFF'):
        status = 'ON' if model.promoter_status[0] == 1 else 'OFF'
        print(f't={sim.curr_simulation_time:7.2f}  promoter={status}  mRNA={model.mRNA_counts[0]}')

simulate_dynamics(model, sim, print_at_each_simulation_step=log_promoter)
print('mRNA produced:', model.mRNA_counts[0])
```

For a multi-gene setup, each gene has its own `(k_on, k_off)` pair. Here gene B is rarely active:

```python
genomic_setup = GenomicSetup(
    chromatin_type='prokaryotic',
    gene_names=['geneA', 'geneB'],
    TSSes=[340.0, 4420.0],
    gene_lengths=[3400.0, 3400.0],
    gene_directions=[1, 1],
    RNAP_on_rates=[0.05, 0.05],
    promoter_mode='non-constitutive',
    TF_on_off_rates=[
        (0.05, 0.01),   # geneA: mostly ON  (~83 % of the time)
        (0.005, 0.05),  # geneB: mostly OFF (~9 % of the time)
    ],
    buffer_length=4420.0,
    are_multiple_chromosomes_present=False,
    regulatory_network_information=read_regulatory_network(['geneA', 'geneB'], None),
)

model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.008,
    TOP2_effective_relaxation_rate=0.008,
)
model = Model(genomic_setup, model_setup)

sim = SimulationSetupAndState(
    simulation_end_mode=0,
    simulation_end_criterion=600.0,
    Gillespie_random_seed=8,
    everything_else_random_seed=8,
)

simulate_dynamics(model, sim)
print('geneA mRNA:', model.mRNA_counts[0])
print('geneB mRNA:', model.mRNA_counts[1])
```

See [User Guide → Genomic setup](user-guide/genomic-setup.md#non-constitutive) and [User Guide → Events and propensities](user-guide/events-and-propensities.md) for the mathematical details of the switching rates.

---

## 18. Multi-chromosome setup

When genes span more than one chromosome, supply a six-column gene file with the chromosome identifier in column 2:

```text
# columns:  name    chromosome   TSS_bp   length_bp   direction   RNAP_on_rate(1/s)
geneA       chrI    10000        5000     1           0.02
geneB       chrI    20000        5000     1           0.02
geneC       chrII   10000        5000    -1           0.02
```

`construct_genomic_setup` detects the extra column, lays the chromosomes end-to-end in reverse order of first appearance (the chromosome appearing last in the file is placed first in the simulation coordinate system), and sets `are_multiple_chromosomes_present=True` along with the computed `chromosomes_end_positions`.

```python
genomic_setup = construct_genomic_setup(
    'multi_chrom.config',
    chromatin_type='prokaryotic',
    promoter_mode='constitutive',
)
```

`Model` then automatically inserts a permanent `'chromosome_barrier'` `BindingProtein` at each internal chromosome boundary. Each barrier has `is_steric_barrier_to_RNAPs=True` and `is_topological_barrier=True` (off-rate 0 — it never unbinds), so it acts as both a physical block to RNAPs and a torsional insulator between chromosomes. The `Lk` list is also re-partitioned: one entry per chromosome rather than one per inter-RNAP segment.

```python
model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.008,
    TOP2_effective_relaxation_rate=0.008,
)
model = Model(genomic_setup, model_setup)

# The last binding-protein entry is the chromosome barrier.
barrier = model.binding_proteins[-1]
print('Barrier name:', barrier.protein_name)               # 'chromosome_barrier'
print('Barrier positions:', model.binding_proteins_positions[-1])

sim = SimulationSetupAndState(
    simulation_end_mode=0,
    simulation_end_criterion=300.0,
    Gillespie_random_seed=19,
    everything_else_random_seed=19,
)
simulate_dynamics(model, sim)
print('Completed transcripts:', sim.RNAPs_finished_transcription)
```

Because the barrier never unbinds, the supercoiling states of the two chromosomes evolve independently: twist generated by RNAPs on chrI cannot propagate to chrII, and vice versa. Use the logging callbacks from Tutorial 11 to inspect the per-chromosome supercoiling profiles.

---

## 19. Gene regulatory network: toggle switch

This tutorial builds on non-constitutive promoters (Tutorial 17) and on the mRNA/protein dynamics enabled by `mRNA_dynamics_mode=1` (Tutorial 9). A **gene regulatory network** couples genes together: the quasi-steady-state protein concentration produced by one gene modulates the promoter ON-rate of another, layered on top of the stochastic TF switching from Tutorial 17. Because the modulation depends on protein levels, `mRNA_dynamics_mode=1` is required so that mRNA and protein concentrations are tracked.

The classic example is a **toggle switch**: two genes that mutually repress each other. Depending on initial conditions and noise, the system settles into one of two stable states — high-A/low-B or low-A/high-B — and can switch stochastically.

**Config file** (`toggle_switch.config`): same two-gene format as Tutorial 4.

**Regulatory network file** (`toggle_switch.topo`): tab-delimited, six columns.

```text
Source	Target	Type	lambda	Theta	n
gene_A	gene_B	2	0.01	2000.0	4.0
gene_B	gene_A	2	0.01	2000.0	4.0
```

`Type=2` (repression). `lambda=0.01` sets a 1 % activity floor. `Theta=2000.0` is the threshold concentration in units of $k_\text{prod}/k_\text{deg}$. `n=4` gives sharp switching (cooperativity).

`construct_genomic_setup` reads both files and wires the network automatically:

```python
import sys
sys.path.append('/path/to/TWISTED/code')

from utilities import construct_genomic_setup
from model_setup import ModelSetup, Model
from simulate_dynamics import SimulationSetupAndState, simulate_dynamics

genomic_setup = construct_genomic_setup(
    'toggle_switch.config',
    chromatin_type='prokaryotic',
    promoter_mode='non-constitutive',
    regulatory_network_file='toggle_switch.topo',
)

model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.008,
    TOP2_effective_relaxation_rate=0.008,
    mRNA_dynamics_mode=1,
    mRNA_degradation_rate=0.05,
    protein_production_rate=1e4,
    protein_degradation_rate=1e2,
)

model = Model(genomic_setup, model_setup)

sim = SimulationSetupAndState(
    simulation_end_mode=0,
    simulation_end_criterion=3000.0,   # 50 minutes
    Gillespie_random_seed=20,
    everything_else_random_seed=20,
)
simulate_dynamics(model, sim)

print('mRNA counts:', model.mRNA_counts)
protein_conc = [
    model.mRNA_counts[i] * model_setup.protein_production_rate / model_setup.protein_degradation_rate
    for i in range(len(genomic_setup.gene_names))
]
print('Protein conc:', protein_conc)
```

The promoter ON-rate for each gene is scaled by the shifted Hill factor:

$$H_i = \prod_j \left[\lambda + (1-\lambda) \frac{1}{1+(p_j^*/\Theta)^n}\right]$$

where $p_j^* = k_\text{prod} \cdot \text{mRNA}_j / k_\text{deg}$ is the quasi-steady-state protein concentration of regulator $j$. Run multiple replicates with different seeds and compare `model.mRNA_counts` between trajectories to observe bimodality.

---

## 20. Saving and resuming a simulation

`model_setup` provides two helpers for checkpointing a run to disk: `save_simulation_state_to_file` serializes the `Model` and the `SimulationSetupAndState` together with `dill`, and `load_simulation_state_from_file` restores them as a `(model, sim)` tuple. This is useful for long runs you want to inspect, share, or continue later.

```python
from model_setup import (
    GenomicSetup, ModelSetup, Model, SimulationSetupAndState,
    save_simulation_state_to_file, load_simulation_state_from_file,
)
from simulate_dynamics import simulate_dynamics

genomic_setup = GenomicSetup(
    chromatin_type='prokaryotic',
    gene_names=['geneA'],
    TSSes=[340.0],
    gene_lengths=[3400.0],
    gene_directions=[1],
    RNAP_on_rates=[0.02],
    promoter_mode='constitutive',
    buffer_length=3400.0,
    are_multiple_chromosomes_present=False,
    regulatory_network_information=read_regulatory_network(['geneA'], None),
)
model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.01,
    TOP2_effective_relaxation_rate=0.005,
)
model = Model(genomic_setup, model_setup)
sim = SimulationSetupAndState(
    simulation_end_mode=1,
    simulation_end_criterion=[50],
    max_RNAPs_to_recruit=[50],
)

# Run to completion, then checkpoint to disk.
simulate_dynamics(model, sim)
save_simulation_state_to_file(model, sim, 'run.pkl')
```

Later — in another script or session — reload the checkpoint and read the results without re-running anything:

```python
model, sim = load_simulation_state_from_file('run.pkl')

print('Completed transcripts:', sim.RNAPs_finished_transcription[0])
print('Total simulated time:', sim.curr_simulation_time)
rates = sim.calculate_RNAP_transcription_rates(model)[0]
print('Mean transcription rate:', sum(rates) / len(rates), 'bp/s')
```

A restored finished run carries `sim.simulation_completed == True`. Calling `simulate_dynamics` on it is a no-op: it prints `This simulation has finished.` and returns immediately, so you cannot accidentally double-count events by re-running a completed checkpoint.

```python
simulate_dynamics(model, sim)   # prints: This simulation has finished.
```

To continue a finished run further, extend the termination criterion with `update_simulation_end_criterion`, which clears the `simulation_completed` flag for you:

```python
model, sim = load_simulation_state_from_file('run.pkl')
sim.update_simulation_end_criterion([100])   # event mode: new per-gene target; or a float for time mode
simulate_dynamics(model, sim)                # resumes and runs until the new criterion is met
```

Because the two random number generators (`rng_Gillespie`, `rng_everything_else`) are stored on `sim`, they are serialized with the checkpoint and restored on load, so a resumed run continues the same random streams rather than restarting them. A run executed in one shot and the same run saved and resumed therefore produce identical trajectories. For time mode, the new end time must not be earlier than `sim.curr_simulation_time`; for event mode, the new counts must require at least one more transcript and stay within `max_RNAPs_to_recruit`.

`SimulationSetupAndState` also accepts `max_wall_time` (default 23 hours), a cap on real (wall-clock) seconds per `simulate_dynamics` call. If this limit is hit before either termination criterion, the loop exits with `sim.simulation_completed` left `False`, leaving the run in the same resumable state as one you interrupted yourself — save it and call `simulate_dynamics` again to pick up where it left off. This is handy for splitting long runs across job-scheduler time limits, independently of the `simulation_end_mode`/`simulation_end_criterion` you are using.

!!! warning
    `load_simulation_state_from_file` uses `dill`, which can execute arbitrary code while deserializing. Only load checkpoint files that you created or otherwise trust.

See [User Guide → Running simulations](user-guide/simulation.md) and the [`model_setup` API reference](api/model-setup.md) for the full signatures.

