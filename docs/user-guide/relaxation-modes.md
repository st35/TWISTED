# Supercoiling relaxation modes

`ModelSetup.supercoiling_relaxation_dynamics_mode` selects the mechanism by which supercoiling is relaxed by topoisomerase-like activity. Six values are accepted; **five are implemented** and one is currently a stub.

| Mode | Status | Per-event behaviour | Required `**kwargs` |
|------|--------|---------------------|--------------------|
| `'global_overall'` | implemented | Reset every segment's `Lk` to `Lk₀` | `global_supercoiling_relaxation_rate` |
| `'global_per_segment'` | implemented | Reset one randomly chosen segment's `Lk` to `Lk₀` (length-weighted) | `global_supercoiling_relaxation_rate` |
| `'global_by_type'` | implemented | One event resets all segments with `σ > 0`; another resets all with `σ < 0` | `local_supercoiling_relaxation_rates` (`[rate_pos, rate_neg]`) |
| `'per_segment_by_type'` | implemented | Pick one segment of the right sign (length-weighted) and reset its `Lk` | `local_supercoiling_relaxation_rates` |
| `'topoisomerase_approximated'` | implemented | One TOP1 event clears twist on a writhe-free segment; one TOP2 event clears writhe on a plectonemic segment | `TOP1_effective_relaxation_rate`, `TOP2_effective_relaxation_rate` |
| `'topoisomerase_based'` | **not implemented** | (planned: explicit TOP1/TOP2 binding/unbinding/catalysis with steric interactions) | `NotImplementedError` raised at `ModelSetup` construction |

Across all modes, the rates are interpreted as **Poisson rates per simulation** (not per segment). No segment-length weighting is applied to the rate itself; length weighting only enters in the selection of *which segment* is affected when an event fires.

---

## `'global_overall'`

A single Poisson event at rate `global_supercoiling_relaxation_rate`. When it fires, **every** segment's linking number is set to its relaxed value:

```
model.Lk[j] = segment_lengths[j] / h_dna   for all j
```

This is the simplest mode and an appropriate initial choice for parameter sweeps that require only that `σ` remain bounded. There is no spatial structure.

```python
ModelSetup(
    supercoiling_relaxation_dynamics_mode='global_overall',
    global_supercoiling_relaxation_rate=0.05,
)
```

---

## `'global_per_segment'`

A single Poisson event at the same rate. When it fires, **one** segment is picked with probability proportional to its length, and only that segment's `Lk` is reset:

```
chosen ~ length-weighted
model.Lk[chosen] = segment_lengths[chosen] / h_dna
```

This adds spatial structure (long segments are relaxed more often) without distinguishing positive from negative supercoiling.

```python
ModelSetup(
    supercoiling_relaxation_dynamics_mode='global_per_segment',
    global_supercoiling_relaxation_rate=0.05,
)
```

---

## `'global_by_type'`

Two independent Poisson events with rates `local_supercoiling_relaxation_rates = [rate_pos, rate_neg]`:

- Event 0 fires at rate `rate_pos` and resets every segment with `σ > 0`.
- Event 1 fires at rate `rate_neg` and resets every segment with `σ < 0`.

This mode is appropriate when positive and negative supercoiling should be assigned different bulk relaxation rates, for instance as a coarse stand-in for gyrase versus TOP1 activity.

```python
ModelSetup(
    supercoiling_relaxation_dynamics_mode='global_by_type',
    local_supercoiling_relaxation_rates=[0.1, 0.05],
)
```

---

## `'per_segment_by_type'`

Two independent Poisson events with the same rates as above, but each event picks **one** segment of the appropriate sign with length weighting and resets only that segment:

```
event 0: pick a length-weighted segment with σ > 0 ; reset it
event 1: pick a length-weighted segment with σ < 0 ; reset it
```

This combines length weighting with sign specificity. It is the most physically motivated of the four non-topoisomerase modes.

```python
ModelSetup(
    supercoiling_relaxation_dynamics_mode='per_segment_by_type',
    local_supercoiling_relaxation_rates=[0.1, 0.05],
)
```

---

## `'topoisomerase_approximated'`

Two independent Poisson events that approximate TOP1 and TOP2 activity without tracking individual enzyme molecules:

- **TOP1 event** at rate `TOP1_effective_relaxation_rate`. Picks a length-weighted segment. **If that segment has no writhe**, set its `Lk` to `Lk₀`. Otherwise the event has no effect (TOP1 cannot resolve plectonemes).
- **TOP2 event** at rate `TOP2_effective_relaxation_rate`. Picks a length-weighted segment. **If that segment has writhe**, set its `Lk` to `Lk₀ × (1 + σ_s)`, where `σ_s` is the plectoneme-formation threshold of that segment (TOP2 reduces writhe but does not over-relax: it stops at the threshold beyond which plectonemes start to form).

The writhe and plectoneme threshold come from the torque-state classification done by [`get_prokaryotic_torque`](../api/biol-methods.md#get_prokaryotic_torque) or [`get_eukaryotic_torque`](../api/biol-methods.md#get_eukaryotic_torque).

```python
ModelSetup(
    supercoiling_relaxation_dynamics_mode='topoisomerase_approximated',
    TOP1_effective_relaxation_rate=0.5,
    TOP2_effective_relaxation_rate=0.2,
)
```

This is the recommended starting point for biologically motivated runs: the TOP1/TOP2 distinction is preserved without the cost of explicit enzyme tracking.

---

## `'topoisomerase_based'`: not yet implemented

The `ModelSetup` constructor raises `NotImplementedError` when this mode is selected. The intended semantics are:

- A discrete copy number of TOP1 and TOP2 molecules.
- Each enzyme binds and unbinds segments stochastically.
- While bound, each enzyme contributes a continuous `dLk/dt` term to its host segment given by the per-molecule rate equations [`get_TOP1_effect_on_Lk_dynamics`](../api/biol-methods.md#get_top1_effect_on_lk_dynamics) and [`get_TOP2_effect_on_Lk_dynamics`](../api/biol-methods.md#get_top2_effect_on_lk_dynamics).
- Steric interactions between bound enzymes and other DNA-bound species are honoured.

See [Not yet implemented](not-yet-implemented.md). Use `'topoisomerase_approximated'` as a substitute.

---

## Choosing a mode

| Requirement | Suggested mode |
|-------------|---------------|
| Simplest possible relaxation; quick parameter sweep | `'global_overall'` |
| Spatial structure without enzyme distinction | `'global_per_segment'` |
| Distinguish positive vs negative supercoiling at the bulk level | `'global_by_type'` |
| Distinguish positive vs negative supercoiling per segment | `'per_segment_by_type'` |
| Distinguish twist (TOP1) vs writhe (TOP2); recommended physical default | `'topoisomerase_approximated'` |
| Full enzyme kinetics with binding sterics | (not yet available) |
