# DNA mechanics

This page is the theoretical reference for everything happening **inside the ODE** that TWISTED integrates: how `Lk` translates into torque, how torque slows the polymerase, how the angle dynamics couple translation to twist injection, how the boundary conditions enter, and how the discrete event clock is coupled to the continuous integration.

Equations follow the conventions in `biol_methods.py` and `model_dynamics.py`. Symbols match the source variable names where practical.

---

## 1. Geometry and notation

The DNA is a one-dimensional rod from `clamp_left = 0` (always) to `clamp_right` (in nm), set by `GenomicSetup` from gene 0 and the buffer length. Positions increase left → right. RNAPs on the `+1` strand move in the `+x` direction; RNAPs on the `−1` strand move in the `−x` direction.

RNAPs and topological-barrier proteins partition the molecule into ordered segments, indexed **right to left** in the state vector. For a segment of length `L_seg`:

- relaxed linking number: `Lk₀ = L_seg / h_dna`, where `h_dna = 2π / w₀`;
- supercoiling density: `σ = (Lk − Lk₀) / Lk₀`;
- torque `τ`, DNA state, writhe fraction, and plectoneme threshold `σ_s` are all functions of `(L_seg, σ, force, kBT, …)` returned by [`get_prokaryotic_torque`](../api/biol-methods.md#get_prokaryotic_torque) or [`get_eukaryotic_torque`](../api/biol-methods.md#get_eukaryotic_torque).

`Lk > Lk₀` ⇒ `σ > 0` ⇒ overwound (positively supercoiled) DNA. `Lk < Lk₀` ⇒ `σ < 0` ⇒ underwound.

---

## 2. Prokaryotic torque law (five-state model)

Implemented in [`get_prokaryotic_torque`](../api/biol-methods.md#get_prokaryotic_torque). Material constants (fixed inside the function):

| Symbol | Value | Meaning |
|--------|-------|---------|
| `A` | 50 nm | Bending persistence length |
| `C` | 95 nm | Twist persistence length |
| `P` | 24 nm | Plectoneme supercoil pitch |
| `A_m` | 4 nm | Bending persistence length, melted DNA |
| `C_m` | 1.75 nm | Twist persistence length, melted DNA |
| `e_m` | `6 kBT` | Denaturation free energy per base |
| `σ₀` | −1 | Reference supercoiling for melting |

Effective stiffnesses and the plectoneme thresholds:

$$c_s = c\left(1 - \frac{C}{4A}\sqrt{\frac{k_BT}{AF}}\right), \qquad c = k_BT \cdot C \cdot \omega_0^2,$$

$$g = F - \sqrt{k_BT F / A}, \qquad p = k_BT \cdot P \cdot \omega_0^2,$$

$$\sigma_s = \frac{1}{c_s}\sqrt{\frac{2 p g}{1 - p/c_s}}, \qquad \sigma_p = \frac{1}{p}\sqrt{\frac{2 p g}{1 - p/c_s}}.$$

Two melting cutoffs use the melted-DNA constants and `g_m = 1.2 (F − √(k_BT F / A_m))`, and follow the closed-form expressions in the source.

The torque is then piecewise linear/constant in five regimes:

| Regime | Code | Range of `σ` | Torque |
|--------|------|--------------|-------|
| Twisted–melted | 0 | `σ ≤ σ_m` | `(c_m / w₀)(σ − σ₀)` |
| Melted | 1 | `σ_m < σ ≤ σ_sm` | `(c_m / w₀)(σ_m − σ₀)` (plateau) |
| Twisted | 2 | `σ_sm < σ ≤ σ_s` | `(c_s / w₀) · σ` |
| Positive plectoneme | 5 | `σ_s < σ ≤ σ_p` | `√(2pg / (1 − p/c_s)) / w₀` (plateau) |
| Twisted plectoneme | 6 | `σ > σ_p` | `(p / w₀) · σ`, hard-capped at 40 pN·nm |

The function also reports a writhe fraction:

- 0 in regimes 0–2,
- `(σ − σ_s) / (σ_p − σ_s)` in regime 5 (linear ramp from 0 to 1 across the plectoneme plateau),
- 1 in regime 6.

The plectoneme threshold `σ_s` is what TOP2 (in `'topoisomerase_approximated'` mode) relaxes to.

### Finite-size correction

When `finite_size_effect_flag = 1`, `σ_s` is rescaled per segment:

$$\sigma_s^{\text{eff}} = \sigma_s \left(1 + \left(\frac{L_0}{L_{\text{seg}}}\right)^2\right),$$

with `L₀ = finite_size_effect_length` (default 340 nm = 1000 bp). Short segments resist plectoneme formation more strongly.

---

## 3. Eukaryotic torque law

Implemented in [`get_eukaryotic_torque`](../api/biol-methods.md#get_eukaryotic_torque). The chromatin torque law has six regimes parameterised by the local nucleosome density `ψ` (the fraction of segment length occupied by nucleosomes; computed by [`get_nucleosome_occupied_fraction_per_segment`](../api/utilities.md#get_nucleosome_occupied_fraction_per_segment)).

Constants (fixed inside the function):

| Symbol | Value | Meaning |
|--------|-------|---------|
| `melted_cutoff` | −0.013 | σ at the melted/twisted boundary |
| `twisted_cutoff` | 0.001 | σ at the twisted/buffering boundary |
| `melted_torque` | −10.0026 pN·nm | torque on the melted plateau |
| `twisted_slope` | 763.064 pN·nm | slope of the twisted regime |
| `pos_twisted_slope` | 753.3442 pN·nm | slope of the positive-twisted regime |
| `plectoneme_cutoff` | 0.0772 | σ at the positive plectoneme onset |

The two **ψ-dependent** cutoffs widen the buffering plateau as nucleosome density grows:

| Cutoff | Expression |
|--------|-----------|
| `buffering_cutoff` | `0.0576 ψ + 0.0013` |
| `pos_twisted_cutoff` | `0.0578 ψ + 0.0205` |

The `pos_twisted_cutoff` plays the role of `σ_s` in the prokaryotic law and is rescaled by the same finite-size correction when the flag is on. The slope of the twisted-plectoneme regime is a quartic polynomial in `ψ`:

```
twisted_plectoneme_slope = 1110.5 ψ⁴ − 1373.6 ψ³ + 770.67 ψ² + 37.0125 ψ + 187.22
```

The piecewise torque is:

| Code | Regime | `σ` range | Torque | Writhe |
|------|--------|----------|-------|-------|
| 1 | melted | `σ < melted_cutoff` | `melted_torque` | 0 |
| 2 | twisted | `melted_cutoff ≤ σ < twisted_cutoff` | `melted_torque + (σ − melted_cutoff) × twisted_slope` | 0 |
| 3 | buffering | `twisted_cutoff ≤ σ < buffering_cutoff` | constant `buffering_torque` | 0 |
| 4 | pos. twisted | `buffering_cutoff ≤ σ < pos_twisted_cutoff` | linear with slope `pos_twisted_slope` | 0 |
| 5 | pos. plectoneme | `pos_twisted_cutoff ≤ σ < plectoneme_cutoff` | constant `plectoneme_torque` | linear ramp 0 → 1 |
| 6 | twisted plectoneme | `σ ≥ plectoneme_cutoff` | linear with slope `twisted_plectoneme_slope`, capped at 40 pN·nm | 1 |

The buffering plateau (regime 3) is the physically interesting feature: it absorbs *positive* supercoiling without raising the torque. Its width is `buffering_cutoff − twisted_cutoff = 0.0576 ψ`, so as RNAPs displace nucleosomes (lowering `ψ`) the buffering capacity shrinks and torque begins to build.

---

## 4. RNAP velocity

For RNAP `i` on a `+1`-strand gene, with torques `τ_b` behind (left) and `τ_f` ahead (right) from segments computed at the current state:

$$v_i \;=\; \frac{v_0}{2}\,\Bigl(1 - \tanh\frac{\tau_f - \tau_b}{\tau_c}\Bigr).$$

For a `−1`-strand gene the sign is flipped (`v_i` becomes negative) and `τ_f`, `τ_b` are taken from the segments to the **left** and **right** of the RNAP respectively. Implemented in [`get_RNAP_velocity`](../api/biol-methods.md#get_rnap_velocity).

### Soft steric ramp

The bare velocity is multiplied by a tanh ramp to represent steric hindrance from the nearest obstacle in the RNAP's direction of travel:

$$v_i \to v_i \cdot \tfrac{1}{2}\Bigl(1 + \tanh\tfrac{s - d}{\lambda}\Bigr),$$

where `s` is the centre-to-centre separation, `d` is the half-sum of the two diameters (or the larger nucleosome formula), and `λ = steric_hindrance_constraint_parameter` (default 2 nm). Implemented in [`get_steric_hindrance_factor`](../api/utilities.md#get_steric_hindrance_factor) and applied in [`get_RNAP_velocities`](../api/model-dynamics.md#get_rnap_velocities).

The exclusion distances are:

| Obstacle | Exclusion distance `d` |
|----------|----------------------|
| Another RNAP | `RNAP_diameter` |
| Nucleosome (steric-barrier) | `(RNAP_diameter + per_nucleosome_DNA_length + nucleosome_linker_length) / 2` |
| Generic protein (steric-barrier) | `(RNAP_diameter + generic_binding_protein_diameter) / 2` |

---

## 5. RNAP angular velocity

The rotational degree of freedom is what couples translation to supercoiling injection. For RNAP `i` at distance `x = |x_i − TSS_i|` from its TSS:

$$\dot\theta_i \;=\; \frac{\omega_0\,\dot x_i \cdot \chi}{\chi + \eta\,x^\alpha} \;+\; \frac{\tau_f - \tau_b}{\chi + \eta\,x^\alpha}.$$

The first term is the rotation that *would* be required to track the translational helical motion at full coupling; the prefactor `χ / (χ + η x^α)` discounts that contribution as the (growing) transcript adds drag. The second term is torsional relaxation driven by the torque imbalance across the RNAP.

Here `τ_f` and `τ_b` are the front and back **with respect to segment geometry**: i.e. the right and left segment torques respectively, with no sign flip for `−1`-strand genes. Implemented in [`get_RNAP_angular_velocity`](../api/biol-methods.md#get_rnap_angular_velocity).

---

## 6. Linking-number dynamics

For an interior segment between two RNAPs, the rate of change of `Lk` is the difference between the angular velocities of the RNAPs at its two ends, divided by `2π`:

$$\dot{Lk}_j \;=\; \frac{1}{2\pi}\,\bigl(\dot\theta_{\text{front}} - \dot\theta_{\text{back}}\bigr),$$

where "front" and "back" refer to the right and left ends of the segment (under the right-to-left segment ordering, `front` is the rightward neighbour of the segment and `back` is the leftward one).

A topological-barrier protein at one end of a segment contributes 0 (it is pinned; its angle does not change).

### Free vs clamped boundaries

The leftmost and rightmost segments touch the clamps. The boundary condition flips the dynamics:

- **Clamped end** (`clamps_status = 'clamped'`): the missing neighbour contributes 0, so the formula above applies with one term zero.
- **Free end** (`clamps_status = 'free'`): twist escapes through the boundary. The `Lk` of the boundary segment is then driven by the *displacement* of the only RNAP at its non-clamp end, not by angular velocity transfer:
  - rightmost free segment: `dLk/dt = (1 / h_dna) · (− dx/dt of leftmost RNAP of that segment)`,
  - leftmost free segment: `dLk/dt = (1 / h_dna) · (dx/dt of rightmost RNAP of that segment)`.

Implemented in [`get_segment_Lk_dynamics`](../api/biol-methods.md#get_segment_lk_dynamics).

### Segment merging when an RNAP exits

When integration takes an RNAP past its gene end, that RNAP is removed by [`update_state_vector_to_remove_dead_RNAPs`](../api/model-dynamics.md#update_state_vector_to_remove_dead_rnaps). The two segments adjacent to the departing RNAP are merged into one and their linking numbers are added; supercoiling is conserved across the merge.

### Lk update on recruitment / topological-barrier (un)binding

When a new RNAP or topological-barrier protein appears, the local segment is split; both new sub-segments inherit the parent `σ`, so `Lk` is partitioned in proportion to the new sub-segment lengths ([`update_Lk_vector_after_RNAP_or_protein_recruitment`](../api/model-dynamics.md#update_lk_vector_after_rnap_or_protein_recruitment)). When a topological barrier unbinds, the two adjacent segments are merged and their linking numbers added ([`update_Lk_vector_after_protein_unbinding`](../api/model-dynamics.md#update_lk_vector_after_protein_unbinding)).

---

## 7. Topoisomerase rate equations (used by the not-yet-implemented `'topoisomerase_based'` mode)

For completeness, the per-molecule continuous rate equations contributed by TOP1 and TOP2 to `dLk/dt` in the planned `'topoisomerase_based'` mode are already implemented in `biol_methods.py`.

### TOP1

TOP1 acts only when the segment writhe fraction is zero. With `n_1` bound TOP1 molecules, `β = 1/k_BT`, `θ = TOP1_theta`, `k_0 = TOP1_k0`, and `τ` the segment torque magnitude:

$$\dot{Lk}_{\text{TOP1}} \;=\;
\begin{cases}
- n_1 k_0 \exp(\theta\beta\tau)\,\bigl(1 - \exp(-2\pi\beta\tau)\bigr) & \sigma > 0,\\[4pt]
+ n_1 k_0 \exp(-\theta\beta\tau)\,\bigl(1 - \exp(-2\pi\beta\tau)\bigr) & \sigma < 0.
\end{cases}$$

The sign drives `Lk` toward `Lk₀`. Implemented in [`get_TOP1_effect_on_Lk_dynamics`](../api/biol-methods.md#get_top1_effect_on_lk_dynamics).

### TOP2

TOP2 acts only when the segment writhe fraction is positive. With writhe `Wr = |σ| · writhe_fraction`, `V_0 = TOP2_V0`, `k_{12} = TOP2_k12`:

$$\dot{Lk}_{\text{TOP2}} \;=\; \mp \, n_2 V_0 \,\frac{Wr}{Wr + k_{12}},$$

with the sign matching the sign of `σ`. Implemented in [`get_TOP2_effect_on_Lk_dynamics`](../api/biol-methods.md#get_top2_effect_on_lk_dynamics).

In the implemented `'topoisomerase_approximated'` mode these per-molecule rates are *not* used; instead, a length-weighted Poisson event with a global effective rate fires and snaps `Lk` to either `Lk₀` (TOP1) or `Lk₀ × (1 + σ_s)` (TOP2). See [Relaxation modes](../user-guide/relaxation-modes.md#topoisomerase_approximated).

---

## 8. RNAP recruitment rate

Implemented in [`get_RNAP_recruitment_rate`](../api/biol-methods.md#get_rnap_recruitment_rate). For TSS `i`:

$$r_i \;=\; \begin{cases} \texttt{RNAP\_on\_rates}[i] & \texttt{promoter\_status}[i] = 1,\\ 0 & \texttt{promoter\_status}[i] = 0. \end{cases}$$

Note that the function takes the local TSS supercoiling density `TSS_sigma` as an argument **but does not currently use it**. There is no σ-dependent or torque-dependent modulation of recruitment rate at present. Possible workarounds include scaling `RNAP_on_rates[i]` externally according to σ values measured in a pilot run, or modifying `get_RNAP_recruitment_rate` directly.

The TSS steric check that gates whether a recruitment event *succeeds* is separate (and σ-independent); see [Events and propensities](../user-guide/events-and-propensities.md#block-0-rnap-recruitment).

---

## 9. mRNA degradation rate

Implemented in [`get_mRNA_degradation_rate`](../api/biol-methods.md#get_mrna_degradation_rate). Per gene:

$$r_i \;=\; \texttt{mRNA\_degradation\_rate} \times \texttt{mRNA\_counts}[i].$$

Active only when `mRNA_dynamics_mode == 1`.

---

## 10. Coupling continuous and discrete dynamics

The state vector handed to the ODE solver carries one extra component `A` beyond the physical degrees of freedom: the cumulative propensity since the start of the integration window.

$$\frac{dA}{dt} = a_0(t) = \sum_\mu r_\mu(t).$$

Before each integration window, the solver is given `A(t_start) = 0`. A uniform `p₀ ~ U(0, 1)` is sampled and the threshold

$$T = \ln(1/p_0)$$

is computed. The integrator advances in chunks of `RNAP_alive_status_check_interval`. After each chunk, the solver output is searched for the first `t_eval` point at which `A(t) > T`. If found, that becomes the event time; the event index is then drawn from the *current* propensities `{r_µ}` with probabilities `r_µ / a₀`. If not found in the current chunk, dead RNAPs are removed and the integrator continues.

This is a Gillespie next-reaction sampling that correctly handles time-varying propensities, because the threshold is on the *time integral* of the total propensity rather than on the instantaneous total. Implemented in [`integrate`](../api/model-dynamics.md#integrate) and dispatched in [`simulate_dynamics`](../api/simulate-dynamics.md#simulate_dynamics).

The check `np.diff(sol.y[-1, :]) >= 0.0` enforces that the cumulative propensity is monotonic; if it ever decreases (which would indicate a numerical pathology), the simulation aborts.

---

## 11. Summary of the dynamic state

Putting it all together, the system of ODEs that is integrated in each window is:

$$
\begin{aligned}
\dot x_i      &= \text{velocity of RNAP } i \text{ (Section 4)},\\
\dot\theta_i  &= \text{angular velocity of RNAP } i \text{ (Section 5)},\\
\dot{Lk}_j    &= \text{linking-number dynamics of segment } j \text{ (Section 6)},\\
\dot A        &= \sum_\mu r_\mu \text{ (Section 10)}.
\end{aligned}
$$

When `A` crosses `ln(1/p_0)`, integration stops and one discrete event from the catalogue in [Events and propensities](../user-guide/events-and-propensities.md) fires.
