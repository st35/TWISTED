# Getting Started

## Installation

TWISTED is a pure-Python package. Clone the repository and install its dependencies:

```bash
git clone https://github.com/st35/TWISTED.git
cd TWISTED/code
pip install scipy numpy
```

No additional installation step is required — simply import the modules directly from the `code/` directory.

### Serving the Documentation Locally

To build and serve the documentation site locally:

```bash
pip install -r docs/requirements.txt
mkdocs serve
```

## Dependencies

| Package | Purpose |
|---------|---------|
| `scipy` | ODE integration (`solve_ivp` / RK45) |
| `numpy` | Numerical array operations |

Both packages are available from PyPI and conda-forge.

---

## Minimal Example

The example below simulates a single constitutively expressed gene for 500 seconds using a simple global supercoiling relaxation model.

```python
from model_setup import GenomicSetup, ModelSetup, Model, SimulationSetupAndState
from simulate_dynamics import simulate_dynamics

# ── 1. Genomic geometry ──────────────────────────────────────────────────────
# All spatial quantities in nm (1 bp = 0.34 nm).
genomic_setup = GenomicSetup(
    chromatin_type='prokaryotic',
    gene_names=['lacZ'],
    TSSes=[340.0],           # TSS at 1000 bp
    gene_lengths=[3400.0],   # 10,000 bp gene
    gene_directions=[1],     # transcribed left → right
    RNAP_on_rates=[0.02],    # recruitment rate (s⁻¹)
    promoter_mode='constitutive',
    buffer_length=3400.0,    # DNA buffer downstream of gene end
)

# ── 2. Model parameters ───────────────────────────────────────────────────────
model_setup = ModelSetup(
    v0=20.0,          # max RNAP velocity (nm/s)
    tau_c=12.0,       # torque sensitivity (pN·nm)
    supercoiling_relaxation_dynamics_mode='global_overall',
    global_supercoiling_relaxation_rate=0.1,
)

# ── 3. Assemble and run ───────────────────────────────────────────────────────
model = Model(genomic_setup, model_setup)
sim = SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode=0,        # 0 = time-based termination
    simulation_end_criterion=500.0,  # run for 500 s
)
simulate_dynamics(model, sim)

# ── 4. Results ────────────────────────────────────────────────────────────────
print('RNAPs finished transcription:', sim.RNAPs_finished_transcription)
print('Recruitment times:', sim.RNAP_recruitment_times)
transcription_rates = sim.calculate_RNAP_transcription_rates(model)
print('Transcription rates (bp/s):', transcription_rates)
```

---

## Event-based Termination

Instead of running for a fixed time, you can stop when a specified number of RNAPs have completed transcription:

```python
sim = SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode=1,              # 1 = event-based termination
    simulation_end_criterion=[50],      # stop after 50 completed RNAPs
)
```

---

## Two-gene Example

```python
genomic_setup = GenomicSetup(
    chromatin_type='prokaryotic',
    gene_names=['geneA', 'geneB'],
    TSSes=[340.0, 4080.0],          # two TSSes
    gene_lengths=[3400.0, 3400.0],
    gene_directions=[1, 1],
    RNAP_on_rates=[0.02, 0.02],
    promoter_mode='constitutive',
    buffer_length=4080.0,            # must extend beyond all genes
)
```

!!! note
    When multiple genes share the same DNA molecule, RNAP steric interactions are automatically enforced. See [Genomic Setup](user-guide/genomic-setup.md) for details on multi-gene configuration constraints.

---

## Next Steps

- [User Guide → Overview](user-guide/overview.md) — Understand the simulation architecture.
- [User Guide → Model Parameters](user-guide/model-setup.md) — Full parameter reference.
- [User Guide → Supercoiling Relaxation Modes](user-guide/relaxation-modes.md) — Choose the right relaxation model for your use case.
- [Theory → DNA Mechanics](theory/dna-mechanics.md) — Physical equations underlying the model.
