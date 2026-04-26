# Model parameters

`ModelSetup` carries every tunable physical and biological constant. The defaults are reasonable for a prokaryotic transcription system; only parameters that differ from the defaults need be overridden, together with a [supercoiling-relaxation mode](relaxation-modes.md).

```python
from model_setup import ModelSetup

model_setup = ModelSetup(
    supercoiling_relaxation_dynamics_mode='global_overall',
    global_supercoiling_relaxation_rate=0.1,
)
```

---

## Constructor

```python
ModelSetup(
    w0: float = 1.85,
    chi: float = 0.05,
    eta: float = 0.0005,
    alpha: float = 1.5,
    v0: float = 20.0,
    tau_c: float = 12.0,
    force: float = 1.0,
    kBT: float = 4.1,
    TOP1_k0: float = 11.0,
    TOP1_theta: float = 0.25,
    TOP2_V0: float = 2.6,
    TOP2_k12: float = 2.0,
    RNAP_diameter: float = 15.0,
    generic_binding_protein_diameter: float = 15.0,
    steric_hindrance_constraint_parameter: float = 2.0,
    clamps_status: tuple[str, str] = ('clamped', 'clamped'),
    finite_size_effect_flag: int = 1,
    supercoiling_relaxation_dynamics_mode: str = 'global_overall',
    mRNA_dynamics_mode: int = 0,
    model_observation_event_rate: float = 0.5,
    **kwargs,
)
```

---

## DNA elastic parameters

| Parameter | Default | Units | Role |
|-----------|--------|------|------|
| `w0` | 1.85 | rad/nm | Intrinsic DNA twist rate; `h_dna = 2π/w0 ≈ 3.4 nm/turn` is computed from this |
| `force` | 1.0 | pN | Applied stretching force entering the torque law |
| `kBT` | 4.1 | pN·nm | Thermal energy at room temperature |

---

## RNAP rotational drag

The effective rotational drag on RNAP at distance `x` from its TSS is `γ(x) = chi + eta · x^alpha`.

| Parameter | Default | Units |
|-----------|--------|------|
| `chi` | 0.05 | pN·nm·s |
| `eta` | 0.0005 | pN·nm^(1−α)·s |
| `alpha` | 1.5 | — |

---

## RNAP elongation

RNAP elongation velocity is the torque-responsive law

$$v = \tfrac{v_0}{2}\bigl(1 - \tanh\tfrac{\tau_f - \tau_b}{\tau_c}\bigr).$$

| Parameter | Default | Units | Role |
|-----------|--------|------|------|
| `v0` | 20.0 | nm/s | Maximum elongation velocity (≈ 59 bp/s) |
| `tau_c` | 12.0 | pN·nm | Sets the half-width of the velocity ramp |

The sign is flipped for `−1`-strand genes. See [Theory → RNAP velocity](../theory/dna-mechanics.md#4-rnap-velocity).

---

## Topoisomerase rate constants

These constants enter the per-molecule TOP1/TOP2 rate equations defined in [`get_TOP1_effect_on_Lk_dynamics`](../api/biol-methods.md#get_top1_effect_on_lk_dynamics) and [`get_TOP2_effect_on_Lk_dynamics`](../api/biol-methods.md#get_top2_effect_on_lk_dynamics). They are used by the not-yet-implemented `'topoisomerase_based'` mode; the simpler `'topoisomerase_approximated'` mode does **not** use them and uses `TOP1_effective_relaxation_rate` / `TOP2_effective_relaxation_rate` instead.

| Parameter | Default | Units | Role |
|-----------|--------|------|------|
| `TOP1_k0` | 11.0 | s⁻¹ | Intrinsic TOP1 nick-close rate |
| `TOP1_theta` | 0.25 | — | Fraction of the twist-relaxation step that occurs in the energy well |
| `TOP2_V0` | 2.6 | s⁻¹ | Maximum TOP2 strand-passage rate |
| `TOP2_k12` | 2.0 | — | Michaelis-like constant for TOP2 writhe dependence |

---

## Steric interactions

Steric exclusion between RNAPs, nucleosomes and other binding proteins is implemented as a **soft tanh ramp** on RNAP velocity:

$$f(s) = \tfrac{1}{2}\bigl(1 + \tanh\tfrac{s - d}{\lambda}\bigr),$$

where `s` is the centre-to-centre separation, `d` is the half-sum of the two diameters (or `(RNAP_diameter + per_nucleosome_DNA_length + nucleosome_linker_length)/2` for nucleosomes), and `λ = steric_hindrance_constraint_parameter`.

| Parameter | Default | Units |
|-----------|--------|------|
| `RNAP_diameter` | 15.0 | nm |
| `generic_binding_protein_diameter` | 15.0 | nm |
| `steric_hindrance_constraint_parameter` | 2.0 | nm |

The same diameters are used by the **TSS steric-hindrance check** when an RNAP recruitment event fires; if anything is within the relevant exclusion distance of the TSS, recruitment fails (unless the obstacle is a displaceable protein; see [Binding proteins](binding-proteins.md)).

---

## Boundary conditions

| Parameter | Default | Allowed | Role |
|-----------|--------|---------|------|
| `clamps_status` | `('clamped', 'clamped')` | each entry `'clamped'` or `'free'` | Torsional boundary conditions at left and right ends |

Internally stored as `left_clamp_status` and `right_clamp_status` (1 = clamped, 0 = free). A free end allows twist to escape, so the `Lk` of the boundary segment is driven by RNAP *displacement* rather than by angular-velocity transfer; see [`get_segment_Lk_dynamics`](../api/biol-methods.md#get_segment_lk_dynamics).

---

## Finite-size correction

| Parameter | Default | Allowed | Role |
|-----------|--------|---------|------|
| `finite_size_effect_flag` | `1` | `0` or `1` | Toggle the correction |
| `finite_size_effect_length` | 340 nm (= 1000 bp) | float, passed via `**kwargs` | Length scale of the correction; only used when the flag is on |

For short segments the plectoneme-formation threshold `σ_s` is raised by

$$\sigma_s^{\text{eff}} = \sigma_s\left(1 + (L_0/L)^2\right),$$

reflecting the difficulty of forming plectonemes in short loops. In eukaryotic mode the same correction is applied to the positive-twist cutoff of the chromatin torque law.

---

## Supercoiling relaxation mode

```python
ModelSetup(supercoiling_relaxation_dynamics_mode=mode, **mode_kwargs)
```

Allowed values for `mode` and the keyword arguments each one demands:

| Mode | Required `**kwargs` |
|------|--------------------|
| `'global_overall'` | `global_supercoiling_relaxation_rate` |
| `'global_per_segment'` | `global_supercoiling_relaxation_rate` |
| `'global_by_type'` | `local_supercoiling_relaxation_rates` (list of two floats) |
| `'per_segment_by_type'` | `local_supercoiling_relaxation_rates` (list of two floats) |
| `'topoisomerase_approximated'` | `TOP1_effective_relaxation_rate`, `TOP2_effective_relaxation_rate` |
| `'topoisomerase_based'` | *(raises `NotImplementedError`)* |

See [Relaxation modes](relaxation-modes.md) for full semantics.

---

## mRNA dynamics

| Parameter | Default | Allowed | Role |
|-----------|--------|--------|------|
| `mRNA_dynamics_mode` | `0` | `0` or `1` | `0` = mRNA only accumulates; `1` = first-order degradation enabled |
| `mRNA_degradation_rate` | required if mode = 1 | float, via `**kwargs` | per-molecule degradation rate (s⁻¹) |

In `mRNA_dynamics_mode == 1` the Gillespie loop adds a per-gene event with rate `mRNA_degradation_rate × mRNA_counts[i]`. The same field exists on `Model` (`model.mRNA_counts`).

---

## Observation event

| Parameter | Default | Units | Role |
|-----------|--------|------|------|
| `model_observation_event_rate` | 0.5 | s⁻¹ | Rate of a dummy "observation" event that does not modify model state |

The observation event ensures that the simulation clock continues to advance even when every biological event rate is zero (e.g. a system with no RNAPs, no proteins, and no relaxation). It must be strictly positive; the constructor asserts this.

---

## Reading parameters back

After construction, every argument is available as a like-named attribute on the instance, along with several derived attributes:

| Attribute | Definition |
|-----------|-----------|
| `model_setup.h_dna` | `2π / w0`, the helical repeat in nm/turn |
| `model_setup.left_clamp_status`, `model_setup.right_clamp_status` | 1 (clamped) or 0 (free) |
| `model_setup.global_supercoiling_relaxation_rate` | rate set if mode is `'global_overall'` or `'global_per_segment'`, else 0 |
| `model_setup.local_supercoiling_relaxation_rates` | `[rate_pos, rate_neg]` if mode is `'global_by_type'` or `'per_segment_by_type'`, else `[0.0, 0.0]` |
| `model_setup.TOP1_effective_relaxation_rate`, `model_setup.TOP2_effective_relaxation_rate` | rates set if mode is `'topoisomerase_approximated'`, else 0 |
| `model_setup.mRNA_degradation_rate` | rate set if mode = 1, else 0 |
