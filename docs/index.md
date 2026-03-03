# TWISTED

**TW**in domain model **I**ntegrating **S**upercoiling from **T**ranscription-topoisomeras**E** **D**ynamics

---

TWISTED is a computational biophysics simulation framework for modeling DNA transcription dynamics with explicit treatment of DNA supercoiling and topoisomerase activity. It uses a hybrid continuous-stochastic approach to capture both the fast mechanical dynamics of RNA polymerase (RNAP) elongation and the slow discrete events such as topoisomerase binding.

## Key Features

- **Twin-domain supercoiling model** — DNA supercoiling is tracked explicitly in segments separated by transcribing RNAPs, using a physically realistic torque-twist relationship for prokaryotic DNA.
- **Hybrid ODE/Gillespie simulation** — Continuous integration (RK45) for transcription mechanics, interleaved with stochastic event selection for discrete biological events.
- **Multiple supercoiling relaxation modes** — Six modes ranging from simple global relaxation to explicit topoisomerase binding/unbinding dynamics with steric interactions.
- **Multi-gene support** — Simulate multiple co-transcribed genes on the same DNA molecule with inter-RNAP steric effects.
- **Configurable promoter models** — Constitutive or non-constitutive (TF-regulated) promoters.
- **mRNA dynamics** — Optional mRNA degradation.

## Quick Start

```python
from model_setup import GenomicSetup, ModelSetup, Model, SimulationSetupAndState
from simulate_dynamics import simulate_dynamics

# Define gene geometry (positions in nm)
genomic_setup = GenomicSetup(
    chromatin_type='prokaryotic',
    gene_names=['geneA'],
    TSSes=[340.0],          # 1000 bp from left clamp
    gene_lengths=[3400.0],  # 10,000 bp gene
    gene_directions=[1],    # positive strand
    RNAP_on_rates=[0.01],   # 1 RNAP per 100 s
    promoter_mode='constitutive',
    buffer_length=3400.0,
)

# Configure physical parameters
model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='global_overall',
    global_supercoiling_relaxation_rate=0.1,
)

# Assemble model and run
model = Model(genomic_setup, model_setup)
sim = SimulationSetupAndState(
    genomic_setup,
    simulation_end_mode=0,       # time-based
    simulation_end_criterion=500.0,
)
simulate_dynamics(model, sim)

print('RNAPs finished:', sim.RNAPs_finished_transcription)
```

## Navigation

| Section | Contents |
|---------|----------|
| [Getting Started](getting-started.md) | Installation and a minimal working example |
| [User Guide](user-guide/overview.md) | In-depth explanations of all components |
| [Theory](theory/dna-mechanics.md) | Physical model and equations |
| [API Reference](api/index.md) | Full function and class reference |
