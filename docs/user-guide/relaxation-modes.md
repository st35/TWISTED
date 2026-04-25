# Supercoiling Relaxation Modes

TWISTED supports six supercoiling relaxation modes controlled by `ModelSetup.supercoiling_relaxation_dynamics_mode`. The choice determines both the physical model of relaxation and whether topoisomerase molecules are explicitly simulated.

---

## Quick Reference

| Mode | Physical interpretation | Requires |
|------|------------------------|---------|
| `global_overall` | One event relaxes the **entire** DNA | `global_supercoiling_relaxation_rate` |
| `global_per_segment` | One event relaxes **one** randomly chosen segment | `global_supercoiling_relaxation_rate` |
| `global_by_type` | Separate events for +σ and −σ across the whole DNA | `local_supercoiling_relaxation_rates` |
| `per_segment_by_type` | Separate events per segment, per supercoiling sign | `local_supercoiling_relaxation_rates` |
| `topoisomerase_approximated` | Effective TOP1/TOP2 rates, segment chosen by length | `TOP1_effective_relaxation_rate`, `TOP2_effective_relaxation_rate` |
| `topoisomerase_based` | Explicit TOP1/TOP2 binding, unbinding, and catalysis with steric interactions | *Not yet implemented* |

---

## Mode Details

### `global_overall`

A single Poisson event at rate `global_supercoiling_relaxation_rate` (s⁻¹) resets the linking number of **every** segment to its relaxed value (`Lk = Lk₀`). This is the simplest model and is appropriate when the identity of the relaxing enzyme is not important.

```python
ModelSetup(
    supercoiling_relaxation_dynamics_mode='global_overall',
    global_supercoiling_relaxation_rate=0.05,
)
```

---

### `global_per_segment`

Same rate as `global_overall`, but each event relaxes a single segment chosen with probability proportional to its fractional length. This captures the length-dependence of topoisomerase accessibility without explicit enzyme tracking.

```python
ModelSetup(
    supercoiling_relaxation_dynamics_mode='global_per_segment',
    global_supercoiling_relaxation_rate=0.05,
)
```

---

### `global_by_type`

Two independent Poisson events: one relaxes all positively supercoiled segments (σ > 0) simultaneously, the other relaxes all negatively supercoiled segments (σ < 0). Rates are specified as a two-element list `[rate_positive, rate_negative]`.

```python
ModelSetup(
    supercoiling_relaxation_dynamics_mode='global_by_type',
    local_supercoiling_relaxation_rates=[0.1, 0.05],
)
```

This mode is useful when positive and negative supercoiling are relaxed by enzymes with different activities (e.g., gyrase vs. TOP1 approximations).

---

### `per_segment_by_type`

Like `global_by_type`, but each event picks a single segment with probability proportional to its length before relaxing. Combines the length-weighting of `global_per_segment` with the sign-specificity of `global_by_type`.

```python
ModelSetup(
    supercoiling_relaxation_dynamics_mode='per_segment_by_type',
    local_supercoiling_relaxation_rates=[0.1, 0.05],
)
```

---

### `topoisomerase_approximated`

Models the *net* effect of topoisomerases as effective relaxation events without tracking individual enzyme copies. Each event:

- **TOP1 event** (rate `TOP1_effective_relaxation_rate`): resets a randomly chosen segment to `Lk = Lk₀` — but **only if** that segment has no writhe (TOP1 cannot act on plectonemic DNA).
- **TOP2 event** (rate `TOP2_effective_relaxation_rate`): resets a plectonemic segment to the plectoneme-formation threshold `Lk = Lk₀(1 + σ_s)` — **only if** that segment has writhe.

```python
ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.5,
    TOP2_effective_relaxation_rate=0.2,
)
```

---

### `topoisomerase_based`

> **Not yet implemented.** This mode is planned but not currently available. Use `topoisomerase_approximated` as an alternative.

---

## Choosing a Mode

| Scenario | Recommended mode |
|----------|-----------------|
| Quick parameter sweeps, minimal complexity | `global_overall` |
| Spatially resolved relaxation, no enzyme tracking | `global_per_segment` or `per_segment_by_type` |
| Biologically motivated TOP1/TOP2 distinction | `topoisomerase_approximated` |
| Full topoisomerase kinetics and steric effects | `topoisomerase_based` |

!!! note "Computational cost"
    `topoisomerase_based` is significantly more expensive per simulation step due to the extra state variables and continuous Lk integration for each bound enzyme. For exploratory work, start with `global_overall` or `topoisomerase_approximated`.
