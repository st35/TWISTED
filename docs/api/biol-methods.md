# `biol_methods` — Biophysical Equations

**File:** `biol_methods.py`
**Imports:** `model_setup`

This module contains the core physical equations for RNAP mechanics, DNA torque, and topoisomerase activity.

---

## `get_RNAP_velocity`

```python
get_RNAP_velocity(
    model: Model,
    gene_index: int,
    left_segment_length: float,
    right_segment_length: float,
    left_torque: float,
    right_torque: float
) -> float
```

Returns the linear velocity of an RNAP in nm/s.

The velocity follows the torque-dependent form:

$$v = \frac{v_0}{2}\left(1 - \tanh\!\frac{\tau_f - \tau_b}{\tau_c}\right)$$

Returns `0.0` if the segment ahead of the RNAP is shorter than `between_RNAPs_steric_effect_cutoff` (steric stalling).

**Parameters:**

| Name | Description |
|------|-------------|
| `gene_index` | Index into `genomic_setup.gene_directions` to resolve which torque is "front" vs. "back" |
| `left_segment_length` | Length (nm) of the segment to the left of the RNAP |
| `right_segment_length` | Length (nm) of the segment to the right of the RNAP |
| `left_torque` | Torque (pN·nm) on the left segment |
| `right_torque` | Torque (pN·nm) on the right segment |

---

## `get_RNAP_angular_velocity`

```python
get_RNAP_angular_velocity(
    model: Model,
    gene_index: int,
    x: float,
    dx_dt: float,
    left_torque: float,
    right_torque: float
) -> float
```

Returns the angular velocity of an RNAP in rad/s:

$$\dot\theta = \frac{\omega_0 \dot{x}\, \chi}{\chi + \eta |x - x_{\mathrm{TSS}}|^\alpha} + \frac{\tau_f - \tau_b}{\chi + \eta |x - x_{\mathrm{TSS}}|^\alpha}$$

**Parameters:**

| Name | Description |
|------|-------------|
| `x` | RNAP position (nm) |
| `dx_dt` | RNAP linear velocity (nm/s) |
| `left_torque` | Torque on the left segment (pN·nm) |
| `right_torque` | Torque on the right segment (pN·nm) |

---

## `get_segment_Lk_dynamics`

```python
get_segment_Lk_dynamics(
    model: Model,
    dtheta_dt_front: float,
    dtheta_dt_back: float
) -> float
```

Returns the rate of change of linking number for a DNA segment (turns/s):

$$\dot{Lk} = \frac{1}{2\pi}(\dot\theta_{\mathrm{front}} - \dot\theta_{\mathrm{back}})$$

---

## `get_RNAP_recruitment_rate`

```python
get_RNAP_recruitment_rate(
    model: Model,
    TSS_index: int,
    promoter_status: int,
    TSS_sigma: float
) -> float
```

Returns the RNAP recruitment rate (s⁻¹) at a given TSS. Returns `0.0` if `promoter_status == 0`. Currently returns `RNAP_on_rates[TSS_index]` when the promoter is ON (supercoiling dependence of recruitment is not yet implemented).

---

## `get_per_TOP1_binding_rate_for_each_segment`

```python
get_per_TOP1_binding_rate_for_each_segment(
    model: Model,
    segments_lengths: list[float],
    segments_sigmas: list[float]
) -> list[float]
```

Returns a list of TOP1 binding rates — one per DNA segment — weighted by segment length:

$$k_{\mathrm{on},j}^{\mathrm{TOP1}} = k_{\mathrm{on}}^{\mathrm{TOP1}} \cdot \frac{L_j}{L_{\mathrm{total}}}$$

---

## `get_per_TOP2_binding_rate_for_each_segment`

```python
get_per_TOP2_binding_rate_for_each_segment(
    model: Model,
    segments_lengths: list[float],
    segments_sigmas: list[float]
) -> list[float]
```

Identical structure to `get_per_TOP1_binding_rate_for_each_segment` but uses `topoisomerase_on_off_rates[1][0]` as the base rate.

---

## `get_TOP1_effect_on_Lk_dynamics`

```python
get_TOP1_effect_on_Lk_dynamics(
    model: Model,
    segment_length: float,
    segment_sigma: float,
    segment_torque: float,
    segment_writhe_frac: float,
    bound_TOP1_count: int
) -> float
```

Returns the contribution of bound TOP1 molecules to $\dot{Lk}$ for a segment. Returns `0.0` if the writhe fraction is non-zero (TOP1 cannot act on plectonemic DNA).

$$\dot{Lk}^{\mathrm{TOP1}} = \begin{cases}
-n \, k_0 e^{\theta\beta\tau}(1 - e^{-2\pi\beta\tau}) & \sigma > 0 \\
+n \, k_0 e^{-\theta\beta\tau}(1 - e^{-2\pi\beta\tau}) & \sigma < 0
\end{cases}$$

---

## `get_TOP2_effect_on_Lk_dynamics`

```python
get_TOP2_effect_on_Lk_dynamics(
    model: Model,
    segment_length: float,
    segment_sigma: float,
    segment_torque: float,
    segment_writhe_frac: float,
    bound_TOP2_count: int
) -> float
```

Returns the contribution of bound TOP2 molecules to $\dot{Lk}$. Returns `0.0` if writhe fraction is zero (TOP2 requires plectonemic DNA):

$$\dot{Lk}^{\mathrm{TOP2}} = \mp n \, V_0 \frac{Wr}{Wr + k_{12}}, \quad Wr = |\sigma| \Phi_w$$

---

## `get_prokaryotic_torque`

```python
get_prokaryotic_torque(
    w0: float,
    force: float,
    kBT: float,
    segment_length: float,
    sigma: float,
    finite_size_effect_flag: int,
    finite_size_effect_length: float
) -> tuple[float, int, float, float]
```

Calculates torque and DNA state for a prokaryotic DNA segment. Returns a 4-tuple:

| Index | Symbol | Type | Description |
|-------|--------|------|-------------|
| 0 | $\tau$ | `float` | Torque (pN·nm) |
| 1 | state | `int` | DNA state code (0=twisted-melted, 1=melted, 2=twisted, 5=pos-plectoneme, 6=twisted-plectoneme) |
| 2 | $\Phi_w$ | `float` | Writhe fraction (0 to 1) |
| 3 | $\sigma_s$ | `float` | Plectoneme-formation threshold |

See [DNA Mechanics](../theory/dna-mechanics.md) for the full torque model derivation.
