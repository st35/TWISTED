# Events and propensities

When debugging, extending the simulator, or writing detailed callbacks, an exhaustive list of the discrete events that the Gillespie loop can fire is useful. This page provides such a reference.

The event rates are assembled by [`get_events_rates`](../api/model-dynamics.md#get_events_rates), which returns a flat list `rates_vector` together with a list of cumulative split points `events_indices`. The dispatch in [`simulate_dynamics`](../api/simulate-dynamics.md#simulate_dynamics) walks the propensities, samples one event index, and routes it to the corresponding handler.

---

## Order of event blocks in the rate vector

```
rates_vector = [
    RNAP_recruitment_rates           per gene,                 # block 0
    model_observation_event_rate,                              # block 1 (single rate)
    global_supercoiling_relaxation_rate,                       # block 2 (single rate)
    local_supercoiling_relaxation_rates  [pos, neg],           # block 3 (two rates)
    [TOP1_effective_relaxation_rate, TOP2_effective_relaxation_rate],   # block 4
    mRNA_degradation_rates           per gene,                 # block 5
    binding_proteins_on_rates        per species (summed),     # block 6
    binding_proteins_off_rates       per species (summed),     # block 7
    promoter_on_rates                per gene,                 # block 8
    promoter_off_rates               per gene,                 # block 9
]
```

Inactive modes/features contribute zero-rate entries in the appropriate block (e.g. `local_supercoiling_relaxation_rates = [0.0, 0.0]` if the mode does not use them). The block boundaries are the cumulative lengths returned in `events_indices[0..9]`.

---

## Block 0: RNAP recruitment

One propensity per gene. The rate for gene `i` is

```
RNAP_on_rates[i]    if promoter_status[i] == 1
0                   if promoter_status[i] == 0
```

When this event fires for gene `i`:

1. Steric hindrance at `TSSes[i]` is checked against all alive RNAPs, then all bound nucleosomes that are steric barriers, then all bound generic steric-barrier proteins.
2. If the only obstacle is another RNAP, recruitment **fails** (no model state change).
3. If the only obstacle is a non-displaceable bound protein, recruitment **fails**.
4. If `max_RNAPs_to_recruit` is set and the per-gene cap is reached, recruitment **fails**.
5. Otherwise, an RNAP is appended at `(x = TSSes[i], θ = 0)`. The local segment's `Lk` is split into two new sub-segments preserving `σ`. The recruitment time is logged in `sim.RNAP_recruitment_times[i]`. If a displaceable steric-barrier protein blocked the TSS, it is removed from `model.binding_proteins_positions`.

---

## Block 1: observation event

A dummy event fired at `model_observation_event_rate` (default 0.5 s⁻¹). It modifies no state. Its sole purpose is to keep the simulation clock advancing in degenerate situations where every other propensity is zero.

---

## Block 2: global supercoiling relaxation

Only non-zero when `supercoiling_relaxation_dynamics_mode` is `'global_overall'` or `'global_per_segment'`. Rate = `global_supercoiling_relaxation_rate`.

- `'global_overall'`: every segment's `Lk` is set to its relaxed value `Lk₀ = L / h_dna`. Sets `last_event_type = 'global_supercoiling_relaxation'`.
- `'global_per_segment'`: one segment is picked with length weighting and its `Lk` is reset. Sets `last_event_type = 'global_supercoiling_relaxation_per_segment'`, with `'_rightmost_segment'` or `'_leftmost_segment'` appended if the chosen segment is the rightmost (`index 0`) or leftmost (`index -1`) segment.

---

## Block 3: type-specific supercoiling relaxation

Two propensities, `local_supercoiling_relaxation_rates = [rate_pos, rate_neg]`. Non-zero only for modes `'global_by_type'` and `'per_segment_by_type'`.

- `'global_by_type'`, event 0: every segment with `σ > 0` is reset to `Lk₀`. Sets `last_event_type = 'per_segment_supercoiling_relaxation_positive_only'`.
- `'global_by_type'`, event 1: every segment with `σ < 0` is reset to `Lk₀`. Sets `last_event_type = 'per_segment_supercoiling_relaxation_negative_only'`.
- `'per_segment_by_type'`, event 0: pick a length-weighted segment **of positive sign** and reset its `Lk`. Sets `last_event_type = 'per_segment_supercoiling_relaxation_positive_only'`.
- `'per_segment_by_type'`, event 1: pick a length-weighted segment **of negative sign** and reset its `Lk`. Sets `last_event_type = 'per_segment_supercoiling_relaxation_negative_only'`.

For the `'per_segment_by_type'` sub-cases, `'_rightmost_segment'` or `'_leftmost_segment'` is appended to `last_event_type` if the chosen segment is the rightmost (`index 0`) or leftmost (`index -1`) segment.

If no segment of the right sign exists, no state change occurs.

---

## Block 4: topoisomerase-approximated activity

Two propensities, `[TOP1_effective_relaxation_rate, TOP2_effective_relaxation_rate]`. Non-zero only when mode is `'topoisomerase_approximated'`.

- **Event 0 (TOP1)**: pick a length-weighted segment. If its writhe fraction is **zero**, set its `Lk` to `Lk₀`. Otherwise no change (TOP1 cannot act on plectonemic DNA). Sets `last_event_type = 'TOP1_supercoiling_relaxation_approximation'`.
- **Event 1 (TOP2)**: pick a length-weighted segment. If its writhe fraction is **positive**, set its `Lk` to `Lk₀ × (1 + σ_s)`, where `σ_s` is the plectoneme-formation threshold of that segment. Otherwise no change. Sets `last_event_type = 'TOP2_supercoiling_relaxation_approximation'`.

For both TOP1 and TOP2, `'_rightmost_segment'` or `'_leftmost_segment'` is appended to `last_event_type` if the chosen segment is the rightmost (`index 0`) or leftmost (`index -1`) segment.

The writhe fraction and plectoneme threshold come from [`get_prokaryotic_torque`](../api/biol-methods.md#get_prokaryotic_torque) (or the eukaryotic equivalent).

---

## Block 5: mRNA degradation

One propensity per gene. Rate for gene `i` = `mRNA_degradation_rate × mRNA_counts[i]`. Non-zero only when `mRNA_dynamics_mode == 1`. When fired, decrement `model.mRNA_counts[i]` by 1.

A degradation event fired on a gene with zero mRNA raises a `ValueError`; this condition should not arise when rates are computed correctly.

---

## Block 6: binding-protein binding (one event type per species)

One aggregated propensity per `BindingProtein` species. The aggregated rate is the sum, over all segments, of

```
n_unbound × basal_on_rate × segment_length × user_on_rate_func(L, σ)
```

When this event fires for species `k`:

1. A segment is chosen with the per-segment on-rates as weights.
2. A position is sampled uniformly inside the chosen segment.
3. The TSS-independent steric check `is_protein_binding_blocked` is run; if blocked, the event is a no-op.
4. Otherwise, the position is appended to `model.binding_proteins_positions[k]`. If the species is a topological barrier, the segment is split into two and the linking number bookkeeping is updated.

---

## Block 7: binding-protein unbinding (one event type per species)

One aggregated propensity per species, equal to the sum over bound molecules of `basal_off_rate * user_off_rate_func(L, σ)` evaluated at each molecule's host segment.

When this event fires for species `k`:

1. A bound molecule is selected with weights equal to its per-molecule off-rate.
2. Its position is removed from `model.binding_proteins_positions[k]`.
3. If the species is a topological barrier, the segments on either side of the unbinding position are merged and their linking numbers are added.

---

## Block 8: promoter ON (one event type per gene)

One propensity per gene, computed by [`get_promoter_on_rate`](../api/biol-methods.md#get_promoter_on_rate). The rate for gene `i` is

```
TF_on_off_rates[i][0]   if promoter_status[i] == 0
0                       if promoter_status[i] == 1
```

When the event fires for gene `i`, `model.promoter_status[i]` is set to 1. Selecting this event for an already-ON promoter raises a `ValueError` (the rate should have been zero).

In `'constitutive'` mode `TF_on_off_rates[i][0]` is `0.0`, so this block contributes nothing.

---

## Block 9: promoter OFF (one event type per gene)

One propensity per gene, computed by [`get_promoter_off_rate`](../api/biol-methods.md#get_promoter_off_rate). The rate for gene `i` is

```
TF_on_off_rates[i][1]   if promoter_status[i] == 1
0                       if promoter_status[i] == 0
```

When the event fires for gene `i`, `model.promoter_status[i]` is set to 0. Selecting this event for an already-OFF promoter raises a `ValueError`.

In `'constitutive'` mode `TF_on_off_rates[i][1]` is `0.0`, so this block contributes nothing.

---

## Inspecting per-block propensities

A common debugging pattern is to print the dominant event family at each step:

```python
def report(model, sim):
    RNAP_gene_index, sv = get_state_vectors_from_dicts(model)
    rates, idx = get_events_rates(model, RNAP_gene_index, sv)
    block_totals = [sum(rates[start:end]) for start, end in zip([0] + idx[:-1], idx)]
    labels = ['recruitment', 'observation', 'global SC', 'by-type SC',
              'topoisomerase', 'mRNA degr.', 'protein on', 'protein off',
              'promoter on', 'promoter off']
    for label, total in zip(labels, block_totals):
        print(f'  {label:>14}: {total:.4g}')

simulate_dynamics(model, sim, print_at_each_simulation_step=report)
```
