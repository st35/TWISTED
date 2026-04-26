# `biol_methods`

Pure scalar physics/biology functions. None of them mutate `model`.

---

## `get_RNAP_velocity`

```python
get_RNAP_velocity(
    model: Model,
    gene_index: int,
    left_segment_length: float,
    right_segment_length: float,
    left_torque: float,
    right_torque: float,
) -> float
```

Bare RNAP linear velocity (nm/s) before steric ramp. For `+1` strand: `(v0 / 2) (1 − tanh((τ_right − τ_left) / τ_c))`. For `−1` strand: sign reversed and front/back roles swapped.

---

## `get_RNAP_angular_velocity`

```python
get_RNAP_angular_velocity(
    model: Model,
    gene_index: int,
    x: float,
    dx_dt: float,
    left_torque: float,
    right_torque: float,
) -> float
```

Angular velocity (rad/s). Returns

```
ω0 · dx_dt · χ / (χ + η · |x − TSS|^α) + (τ_right − τ_left) / (χ + η · |x − TSS|^α)
```

with the appropriate signs for the gene direction.

See [Theory → RNAP angular velocity](../theory/dna-mechanics.md#5-rnap-angular-velocity).

---

## `get_segment_Lk_dynamics`

```python
get_segment_Lk_dynamics(
    model: Model,
    dx_dt_front: float, dx_dt_back: float,
    dtheta_dt_front: float, dtheta_dt_back: float,
    left_barrier_type: str, right_barrier_type: str,
    is_rightmost_segment: bool, is_leftmost_segment: bool,
) -> float
```

Per-segment `dLk/dt`. Free clamp ends contribute `(±1/h_dna) · dx_dt` (twist escapes through the boundary); clamped ends contribute the angular-velocity difference of the two end agents divided by `2π`. Topological-barrier proteins contribute 0 to whichever end they pin.

See [Theory → Linking-number dynamics](../theory/dna-mechanics.md#6-linking-number-dynamics).

---

## `get_RNAP_recruitment_rate`

```python
get_RNAP_recruitment_rate(
    model: Model,
    TSS_index: int,
    promoter_status: int,
    TSS_sigma: float,
) -> float
```

Returns `RNAP_on_rates[TSS_index]` if `promoter_status == 1`, else `0`. `TSS_sigma` is currently unused.

---

## `get_TOP1_effect_on_Lk_dynamics`

```python
get_TOP1_effect_on_Lk_dynamics(
    model: Model,
    segment_length: float, segment_sigma: float,
    segment_torque: float, segment_writhe_frac: float,
    bound_TOP1_count: int,
) -> float
```

Continuous per-molecule TOP1 contribution to `dLk/dt`. Returns 0 if `segment_writhe_frac > 0` (TOP1 cannot act on plectonemic DNA). Sign drives `Lk` toward `Lk₀`.

Used only by the planned `'topoisomerase_based'` mode; **not** invoked in `'topoisomerase_approximated'`.

---

## `get_TOP2_effect_on_Lk_dynamics`

```python
get_TOP2_effect_on_Lk_dynamics(
    model: Model,
    segment_length: float, segment_sigma: float,
    segment_torque: float, segment_writhe_frac: float,
    bound_TOP2_count: int,
) -> float
```

Continuous per-molecule TOP2 contribution to `dLk/dt`. Active only when `segment_writhe_frac > 0`. Magnitude follows Michaelis–Menten in writhe; sign matches the sign of `σ`.

---

## `get_mRNA_degradation_rate`

```python
get_mRNA_degradation_rate(model: Model, mRNA_count: int) -> float
```

Returns `mRNA_degradation_rate × mRNA_count` if `mRNA_dynamics_mode == 1`, else 0.

---

## `get_prokaryotic_torque`

```python
get_prokaryotic_torque(
    w0: float, force: float, kBT: float,
    segment_length: float, sigma: float,
    finite_size_effect_flag: int,
    finite_size_effect_length: float,
) -> tuple[float, int, float, float]
```

Returns `(torque, dna_state, writhe_fraction, sigma_s)` for a prokaryotic segment.

| Output | Meaning |
|--------|--------|
| `torque` | pN·nm |
| `dna_state` | 0 = twisted-melted, 1 = melted, 2 = twisted, 5 = positive plectoneme, 6 = twisted plectoneme |
| `writhe_fraction` | 0 (no plectoneme), `(σ−σ_s)/(σ_p−σ_s)` (forming), 1 (full plectoneme) |
| `sigma_s` | Plectoneme-formation threshold (after finite-size rescaling) |

See [Theory → Prokaryotic torque law](../theory/dna-mechanics.md#2-prokaryotic-torque-law-five-state-model).

---

## `get_eukaryotic_torque`

```python
get_eukaryotic_torque(
    force: float, segment_length: float,
    psi: float, sigma: float,
    finite_size_effect_flag: int,
    finite_size_effect_cutoff: float,
) -> tuple[float, int, float, float]
```

Returns `(torque, dna_state, writhe_fraction, pos_twisted_cutoff)` for a eukaryotic segment with nucleosome occupancy fraction `psi`.

| Output | Meaning |
|--------|--------|
| `dna_state` | 1 = melted, 2 = twisted, 3 = buffering, 4 = positive twisted, 5 = positive plectoneme, 6 = twisted plectoneme |
| `pos_twisted_cutoff` | Plays the role of `sigma_s` for TOP2-style relaxation |

See [Theory → Eukaryotic torque law](../theory/dna-mechanics.md#3-eukaryotic-torque-law).
