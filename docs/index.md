# TWISTED

**TW**in domain model **I**ntegrating **S**upercoiling from **T**ranscription–topoisomeras**E** **D**ynamics.

TWISTED is a Python framework for simulating transcription on a one-dimensional DNA molecule with **explicit DNA supercoiling**. It tracks where each RNA polymerase (RNAP) is, how fast it elongates, the linking number of the DNA between consecutive RNAPs, and the discrete biochemical events (recruitment, supercoiling relaxation, topoisomerase action, mRNA degradation, protein binding/unbinding) that drive the system.

## What problem does it solve?

Transcribing RNAPs continuously inject positive supercoils ahead of themselves and negative supercoils behind. The resulting torsional stress slows or stalls elongation, modulates RNAP recruitment, and is relaxed by topoisomerases. Capturing this two-way coupling in a tractable simulation requires combining:

- **Continuous mechanics** — RNAP positions, RNAP angular positions, and segment linking numbers evolve through a system of ODEs.
- **Discrete stochastic events** — recruitment, topoisomerase action, mRNA degradation and protein binding are sampled with a Gillespie-style algorithm.

TWISTED interleaves the two: an adaptive ODE solver (RK23 by default) advances the continuous variables while accumulating a *cumulative propensity*; when the propensity crosses an exponentially distributed threshold, integration stops and a discrete event is selected and executed.

## Supported configurations

- Single or multiple co-transcribed genes on the same DNA, in any orientation (tandem, convergent, divergent).
- Five supercoiling-relaxation modes ranging from a single global rate to length-weighted TOP1/TOP2-like activity gated by writhe.
- Generic DNA-binding proteins with user-supplied σ-dependent on/off rates; selected species may act as topological barriers (splitting the DNA into independent supercoiling domains) or be displaced from a TSS by an incoming RNAP.
- Eukaryotic chromatin with dynamically (un)binding nucleosomes and a chromatin-specific torque law.
- Optional first-order mRNA degradation.
- Clamped or free torsional boundary conditions at each end of the DNA.

## Minimal example

```python
from model_setup import GenomicSetup, ModelSetup, Model, SimulationSetupAndState
from simulate_dynamics import simulate_dynamics

genomic_setup = GenomicSetup(
    chromatin_type='prokaryotic',
    gene_names=['geneA'],
    TSSes=[340.0],            # 1000 bp from the left clamp (1 bp = 0.34 nm)
    gene_lengths=[3400.0],    # 10 kb gene
    gene_directions=[1],      # transcribed left → right
    RNAP_on_rates=[0.02],     # 1/s
    promoter_mode='constitutive',
    buffer_length=3400.0,
)

model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='global_overall',
    global_supercoiling_relaxation_rate=0.1,
)

model = Model(genomic_setup, model_setup)
sim = SimulationSetupAndState(
    simulation_end_mode=0,           # time-based
    simulation_end_criterion=500.0,  # seconds
)
simulate_dynamics(model, sim)

print('RNAPs that finished transcription:', sim.RNAPs_finished_transcription)
print('Mean transcription rates (bp/s) per gene:',
      [sum(r)/len(r) if r else None for r in sim.calculate_RNAP_transcription_rates(model)])
```

## Documentation map

| Section | Contents |
|---------|----------|
| [Getting Started](getting-started.md) | Installation, the four-object construction pattern, and a guided minimal example |
| [Tutorials](tutorials.md) | Sixteen worked examples, one per feature, in order of increasing complexity |
| [User Guide](user-guide/overview.md) | Reference descriptions of every component |
| [Theory](theory/dna-mechanics.md) | Physical model, governing equations, and numerical method |
| [API Reference](api/index.md) | Public functions and classes |

## Conventions

- All spatial quantities are in **nanometres**. Conversion: `1 bp = 0.34 nm`. Config-file inputs are in bp and converted by [`construct_genomic_setup`](api/utilities.md).
- All times are in **seconds**, all rates in **s⁻¹**, torque in **pN·nm**, force in **pN**, and energy in **pN·nm** (so `kBT ≈ 4.1 pN·nm` at room temperature).
- Linking number `Lk` is dimensionless (turns); supercoiling density `σ = (Lk − Lk₀) / Lk₀`.
