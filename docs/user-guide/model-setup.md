# Model Parameters

`ModelSetup` holds all tunable physical and biological parameters. Default values correspond to well-characterised prokaryotic transcription parameters from the literature.

---

## Constructor Signature

```python
ModelSetup(
    w0=1.85, chi=0.05, eta=0.0005, alpha=1.5,
    v0=20.0, tau_c=12.0, force=1.0, kBT=4.1,
    TOP1_k0=11.0, TOP1_theta=0.25,
    TOP2_V0=2.6, TOP2_k12=2.0,
    RNAP_diameter=15.0,
    TOPO_diameter=15.0,
    generic_binding_protein_diameter=15.0,
    steric_hindrance_constraint_parameter=2.0,
    clamps_status=('clamped', 'clamped'),
    finite_size_effect_flag=1,
    supercoiling_relaxation_dynamics_mode='global_overall',
    mRNA_dynamics_mode=0,
    model_observation_event_rate=0.5,
    **kwargs
)
```

---

## DNA Parameters

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `w0` | 1.85 | rad/nm | Intrinsic DNA twist rate (relaxed DNA); related to helical pitch `h = 2π/w0 ≈ 3.4 nm/turn` |
| `chi` | 0.05 | pN·nm·s | Rotational drag coefficient of RNAP body |
| `eta` | 0.0005 | pN·nm^(1−α)·s | Pre-factor for distance-dependent drag term |
| `alpha` | 1.5 | — | Exponent for distance-dependent drag |
| `force` | 1.0 | pN | Applied stretching force on DNA |
| `kBT` | 4.1 | pN·nm | Thermal energy at room temperature (~298 K) |

The effective rotational drag experienced by an RNAP at distance `x` from its TSS is:

$$\gamma(x) = \chi + \eta \, x^\alpha$$

---

## RNAP Parameters

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `v0` | 20.0 | nm/s | Maximum RNAP elongation velocity (≈ 59 bp/s) |
| `tau_c` | 12.0 | pN·nm | Torque sensitivity constant — sets the scale for torque-dependent velocity reduction |

RNAP velocity is:

$$v = \frac{v_0}{2}\left(1 - \tanh\!\frac{\tau_f - \tau_b}{\tau_c}\right)$$

where $\tau_f$ and $\tau_b$ are the torques in front of and behind the RNAP (see [DNA Mechanics](../theory/dna-mechanics.md)).

---

## Topoisomerase Parameters

These parameters govern the action of Type I (TOP1) and Type II (TOP2) topoisomerases.

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `TOP1_k0` | 11.0 | s⁻¹ | Intrinsic TOP1 nick-close rate |
| `TOP1_theta` | 0.25 | — | Fraction of the twist-relaxation step done while the DNA is in the energy well |
| `TOP2_V0` | 2.6 | s⁻¹ | Maximum TOP2 strand-passage rate |
| `TOP2_k12` | 2.0 | — | Michaelis-like constant for TOP2 writhe dependence |

---

## Steric Interaction Parameters

Steric interactions between molecular species are modelled as a **soft constraint**. Rather than zeroing an RNAP's velocity at a hard exclusion distance, the velocity is continuously reduced by a smooth factor:

$$f(s) = \frac{1}{2}\left(1 + \tanh\!\frac{s - d}{\lambda}\right)$$

where $s$ is the centre-to-centre separation, $d$ is the exclusion distance `(d1 + d2) / 2`, and $\lambda$ is the `steric_hindrance_constraint_parameter`. The factor approaches 0 when the obstacle is very close ($s \ll d$) and 1 when far away ($s \gg d$). The transition spans approximately $\pm 2\lambda$ around $d$.

Nucleosomes are a special case: their physical extent is `per_nucleosome_DNA_length + nucleosome_linker_length` (from `GenomicSetup`).

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `RNAP_diameter` | 15.0 | nm | RNAP physical diameter; used for RNAP–RNAP, RNAP–protein, and RNAP–nucleosome steric checks |
| `TOPO_diameter` | 15.0 | nm | Topoisomerase physical diameter; used for TSS–TOPO steric checks |
| `generic_binding_protein_diameter` | 15.0 | nm | Default diameter for non-nucleosome binding proteins; used for RNAP–protein and protein–protein steric checks |
| `steric_hindrance_constraint_parameter` | 2.0 | nm | Controls the width of the soft steric transition zone in the tanh ramp |

---

## Boundary Conditions

| Parameter | Default | Values | Description |
|-----------|---------|--------|-------------|
| `clamps_status` | `('clamped', 'clamped')` | `('clamped'` or `'free', 'clamped'` or `'free')` | Whether the left and right ends of the DNA are torsionally clamped or free. Free ends do not resist torsional stress. Internally stored as `left_clamp_status` and `right_clamp_status` (1 = clamped, 0 = free) |

---

## Finite-size Effects

| Parameter | Default | Values | Description |
|-----------|---------|--------|-------------|
| `finite_size_effect_flag` | `1` | `0` or `1` | Enable (`1`) or disable (`0`) finite-size corrections to the plectoneme formation threshold |
| `finite_size_effect_length` | 340 nm | nm (float) | Length scale for the finite-size correction (default: 1000 bp = 340 nm). Only relevant when flag is `1`. Pass via `**kwargs` |

The finite-size correction raises the plectoneme threshold $\sigma_s$ for short DNA segments, reflecting the fact that very short segments cannot form plectonemes as readily:

$$\sigma_s^{\text{eff}} = \sigma_s \left(1 + \left(\frac{L_0}{L}\right)^2\right)$$

where $L_0$ is `finite_size_effect_length` and $L$ is the segment length.

---

## Supercoiling Relaxation Mode

See the dedicated page: [Supercoiling Relaxation Modes](relaxation-modes.md).

---

## mRNA Dynamics

| Parameter | Default | Values | Description |
|-----------|---------|--------|-------------|
| `mRNA_dynamics_mode` | `0` | `0` or `1` | `0`: mRNAs accumulate without degradation; `1`: mRNAs degrade at rate `mRNA_degradation_rate` |
| `mRNA_degradation_rate` | — | s⁻¹ | Required when `mRNA_dynamics_mode=1`; pass via `**kwargs` |

---

## Observation Event Rate

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `model_observation_event_rate` | 0.5 | s⁻¹ | Rate of dummy "observation" events. These events do not change the model state but ensure the simulation clock advances even when all biological event rates are zero. Must be > 0 |

---

## Example: Topoisomerase-based Mode

> **Not yet implemented.** The `topoisomerase_based` mode is planned but not currently available. Use `topoisomerase_approximated` as an alternative.

```python
model_setup = ModelSetup(
    v0=20.0,
    tau_c=12.0,
    supercoiling_relaxation_dynamics_mode='topoisomerase_based',
    topoisomerase_copy_numbers=[10, 5],       # 10 × TOP1, 5 × TOP2
    topoisomerase_on_off_rates=[
        (0.1, 0.05),   # TOP1: (k_on, k_off) in s⁻¹
        (0.05, 0.02),  # TOP2: (k_on, k_off) in s⁻¹
    ],
)
```
