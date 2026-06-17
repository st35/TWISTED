# `utilities`

Helper functions used across the package.

---

## File and setup helpers

### `read_genes_information`

```python
read_genes_information(filename: str) -> tuple[list[str], list[float], list[float], list[int], list[float]]
```

Read a tab-delimited gene file with columns `gene_name`, `TSS_bp`, `gene_length_bp`, `direction`, `RNAP_on_rate`. Returns `(gene_names, TSSes_nm, gene_lengths_nm, gene_directions, RNAP_on_rates)`. Positions and lengths are converted from bp to nm.

### `construct_genomic_setup`

```python
construct_genomic_setup(
    filename: str,
    chromatin_type: str,
    promoter_mode: str = 'constitutive',
    buffer_length: float = 10000.0 * 0.34,
    **kwargs,
) -> GenomicSetup
```

Convenience: read a gene file and construct a `GenomicSetup` in one call.

---

## Sampling

### `select_event_based_on_propensities`

```python
select_event_based_on_propensities(rates_vector: list[float], p: float) -> int | None
```

Pick the smallest index `k` such that `cumulative_rates[k] >= p × sum(rates_vector)`. Returns `None` only on a degenerate empty-rates input.

---

## Spatial lookups

### `get_spot_segment_index`

```python
get_spot_segment_index(spot: float, segments_lengths: list[float]) -> int
```

Returns the right-to-left segment index containing `spot`.

### `get_nucleosome_occupied_fraction_per_segment`

```python
get_nucleosome_occupied_fraction_per_segment(model, segments_lengths, segment_index) -> float
```

Fraction of segment length covered by bound nucleosomes (range `[0, 1]`). Used by the eukaryotic torque law.

### `get_ordering_of_RNAPs_and_proteins`

```python
get_ordering_of_RNAPs_and_proteins(model, RNAP_gene_index, state_vector, topological_barrier_proteins_only: bool = False) -> tuple[list[float], list[str]]
```

Returns `(positions, ids)` sorted left to right. `ids` are `'RNAP'` or the species index as a string.

### `get_nearest_steric_obstacles_for_RNAP`

```python
get_nearest_steric_obstacles_for_RNAP(RNAP_index, RNAP_pos, sorted_positions, sorted_ids) -> tuple[float|None, str|None, float|None, str|None]
```

Returns `(left_pos, left_id, right_pos, right_id)` for the nearest obstacles on either side; `None` if no obstacle on that side.

---

## Steric checks

### `get_TSS_steric_hindrance_status`

```python
get_TSS_steric_hindrance_status(model, TSS_position, RNAP_gene_index, state_vector) -> tuple[int, float, int]
```

Returns `(status, blocking_position, blocking_id)`. `status` is 1 if any RNAP, nucleosome, or generic steric-barrier protein is within the TSS exclusion zone, else 0. `blocking_id` is `-1` for an RNAP, otherwise the binding-protein species index. Used by the recruitment handler.

### `is_protein_binding_blocked`

```python
is_protein_binding_blocked(model, RNAP_gene_index, state_vector, protein_index, binding_position) -> int
```

Returns 1 if a candidate binding event for species `protein_index` at `binding_position` would overlap an existing RNAP, nucleosome, or other protein.

### `get_steric_hindrance_factor`

```python
get_steric_hindrance_factor(model, separation: float, steric_hindrance_distance: float) -> float
```

Soft tanh ramp `0.5 (1 + tanh((separation − steric_hindrance_distance) / λ))` with `λ = steric_hindrance_constraint_parameter`.

### `check_separation_between_nucleosomes`

```python
check_separation_between_nucleosomes(model: Model) -> None
```

Asserts that bound nucleosomes are pairwise separated by at least `per_nucleosome_DNA_length + nucleosome_linker_length`.

### `check_separation_between_nucleosomes_and_RNAPs`

```python
check_separation_between_nucleosomes_and_RNAPs(model, RNAP_gene_index, state_vector) -> None
```

Asserts that no RNAP overlaps a bound nucleosome.

---

## Misc

### `find_and_remove_from_list`

```python
find_and_remove_from_list(lst: list[float], value: float, tolerance: float = 1e-6) -> None
```

Remove the first occurrence of `value` (within tolerance) from `lst`. Raises if not found.

### `print_list`

```python
print_list(name: str, lst: list[float]) -> None
```

Tabbed printout for debugging.
