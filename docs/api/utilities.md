# `utilities` — I/O and Helper Functions

**File:** `utilities.py`
**Imports:** `model_setup`

---

## File I/O

### `read_genes_information`

```python
read_genes_information(filename: str) -> tuple[
    list[str],    # gene_names
    list[float],  # TSSes (nm)
    list[float],  # gene_lengths (nm)
    list[int],    # gene_directions
    list[float],  # RNAP_on_rates (s⁻¹)
]
```

Reads gene data from a **tab-delimited file** (no header). Each line defines one gene:

```
<name>  <TSS_bp>  <length_bp>  <direction>  <RNAP_on_rate>
```

Columns 2 and 3 (TSS and length) are multiplied by `0.34` to convert bp → nm. Direction must be `+1` or `-1`.

**Example file:**
```
lacZ    1000    10000   1   0.02
lacA    12000   3000    1   0.01
```

**Example usage:**
```python
names, TSSes, lengths, dirs, rates = read_genes_information('genes.tsv')
```

---

### `construct_genomic_setup`

```python
construct_genomic_setup(
    filename: str,
    chromatin_type: str,
    promoter_mode: str = 'constitutive',
    buffer_length: float = 3400.0,   # 10,000 bp in nm
    **kwargs
) -> GenomicSetup
```

Factory function that reads a gene file and returns a configured `GenomicSetup`.

**Keyword arguments:**

| Keyword | Type | Description |
|---------|------|-------------|
| `explicit_RNAP_on_rates` | `list[float]` | Multipliers applied element-wise to the file's RNAP on-rates |
| `TF_on_off_rates` | `list[tuple[float,float]]` | Required when `promoter_mode='non-constitutive'` |

All eukaryotic keyword arguments accepted by `GenomicSetup` can also be passed through this function (see [Eukaryotic Keyword Arguments](../user-guide/genomic-setup.md#eukaryotic-keyword-arguments)):

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `per_nucleosome_DNA_length` | `float` | 147 (bp) | DNA wrapped per nucleosome (converted to nm internally) |
| `nucleosome_linker_length` | `float` | 30 (bp) | Linker DNA between nucleosomes (converted to nm internally) |
| `nucleosomes_are_steric_barriers_to_RNAPs` | `bool` | `True` | Whether nucleosomes block RNAP passage |
| `nucleosome_count` | `int` | auto | Explicit nucleosome count; if omitted, computed by tiling |
| `nucleosomes_can_be_displaced_at_TSS_by_RNAP` | `bool` | `False` | Whether nucleosomes blocking a TSS can be displaced by an incoming RNAP |

**Example:**
```python
from utilities import construct_genomic_setup

genomic_setup = construct_genomic_setup(
    'genes.tsv',
    chromatin_type='prokaryotic',
    promoter_mode='constitutive',
    buffer_length=3400.0,
    explicit_RNAP_on_rates=[2.0],   # double the rate from the file
)
```

---

## Random Numbers

### `uniform_random_in_interval`

```python
uniform_random_in_interval(start: float, end: float) -> float
```

Returns a uniform random float in the half-open interval `[start, end)`.

---

## Geometry

### `get_spot_segment_index`

```python
get_spot_segment_index(spot: float, segments_lengths: list[float]) -> int
```

Returns the index of the DNA segment containing `spot` (given as an absolute position measured from the **left** clamp). Segments are ordered right-to-left.

Returns `-1` if the spot falls outside all segments (which should not occur in a correctly configured simulation).

---

## Steric Hindrance

### `get_TSS_steric_hindrance_status`

```python
get_TSS_steric_hindrance_status(
    model: Model,
    TSS_position: float,
    RNAP_gene_index: list[int],
    state_vector: list[float]
) -> tuple[int, float, int]
```

Returns a 3-tuple `(status, position, entity_id)`:

- `status`: `1` if the TSS at `TSS_position` is sterically blocked, `0` otherwise.
- `position`: the position (nm) of the blocking entity, or `0.0` if unblocked.
- `entity_id`: identifies the blocking entity type — `-1` for an RNAP, `0` to `len(model.binding_proteins) - 1` for a bound protein (index into `model.binding_proteins`), `len(model.binding_proteins)` for a topoisomerase, or `-1` when unblocked.

Blocking is triggered by:

- Any active RNAP within `RNAP_diameter` nm of the TSS.
- Any bound topoisomerase within `(RNAP_diameter + TOPO_diameter) / 2` nm of the TSS (in `topoisomerase_based` mode only).
- Any bound nucleosome (with `is_steric_barrier_to_RNAPs=True`) within `(RNAP_diameter + per_nucleosome_DNA_length + nucleosome_linker_length) / 2` nm of the TSS.
- Any other bound protein (with `is_steric_barrier_to_RNAPs=True`) within `(RNAP_diameter + generic_binding_protein_diameter) / 2` nm of the TSS.

### `is_protein_binding_blocked`

```python
is_protein_binding_blocked(
    model: Model,
    RNAP_gene_index: list[int],
    state_vector: list[float],
    protein_index: int,
    binding_position: float
) -> int
```

Returns `1` if a binding protein of the given type cannot bind at `binding_position`, `0` otherwise.

Blocking is triggered by:

- Any RNAP within `(RNAP_diameter + size) / 2` nm, where `size` is `per_nucleosome_DNA_length + nucleosome_linker_length` for nucleosomes or `generic_binding_protein_diameter` for other proteins.
- For nucleosomes: any other bound nucleosome within `per_nucleosome_DNA_length + nucleosome_linker_length` nm.
- For non-nucleosome proteins: any bound protein within `generic_binding_protein_diameter` nm.

### `is_TOPO_binding_blocked`

> *Not yet implemented.* This function is defined but the `topoisomerase_based` mode that uses it is not currently available.

```python
is_TOPO_binding_blocked(
    model: Model,
    RNAP_gene_index: list[int],
    state_vector: list[float],
    binding_position: float
) -> int
```

Returns `1` if a topoisomerase cannot bind at `binding_position` because an RNAP is within `(RNAP_diameter + TOPO_diameter) / 2` nm, `0` otherwise.

### `get_nucleosome_occupied_fraction_per_segment`

```python
get_nucleosome_occupied_fraction_per_segment(
    model: Model,
    segments_lengths: list[float],
    segment_index: int
) -> float
```

Returns the fraction of a DNA segment occupied by nucleosomes (0.0 to 1.0). Computes geometric overlap between each bound nucleosome's footprint (`per_nucleosome_DNA_length + nucleosome_linker_length`) and the segment boundaries, then divides total occupied length by segment length.

### `check_separation_between_nucleosomes`

```python
check_separation_between_nucleosomes(model: Model) -> None
```

Validation function that raises `ValueError` if any two bound nucleosomes are closer than `per_nucleosome_DNA_length + nucleosome_linker_length` nm (centre-to-centre).

### `check_separation_between_nucleosomes_and_RNAPs`

```python
check_separation_between_nucleosomes_and_RNAPs(
    model: Model,
    RNAP_gene_index: list[int],
    state_vector: list[float]
) -> None
```

Validation function that raises `ValueError` if any RNAP is closer than `(RNAP_diameter + per_nucleosome_DNA_length) / 2` nm to a bound nucleosome (centre-to-centre). This uses a deliberately weaker threshold (core footprint only, no linker) compared to the steric enforcement checks.

### `get_ordering_of_RNAPs_and_proteins`

```python
get_ordering_of_RNAPs_and_proteins(
    model: Model,
    RNAP_gene_index: list[int],
    state_vector: list[float]
) -> tuple[list[float], list[str]]
```

Merges all RNAP positions and bound steric-barrier protein positions into a single list sorted by position. Returns `(sorted_positions, sorted_ids)` where each ID is either `'RNAP_<index>'` for an RNAP or `'<binding_protein_index>'` for a bound protein (matching the index into `model.binding_proteins`).

### `get_nearest_steric_obstacles_for_RNAP`

```python
get_nearest_steric_obstacles_for_RNAP(
    RNAP_index: int,
    RNAP_pos: float,
    sorted_positions: list[float],
    sorted_ids: list[str]
) -> tuple[float | None, str | None, float | None, str | None]
```

Given the sorted ordering from `get_ordering_of_RNAPs_and_proteins`, finds the nearest obstacle to the left and right of the specified RNAP. Returns `(left_position, left_id, right_position, right_id)`. Any value is `None` if no obstacle exists on that side.

### `get_steric_hindrance_factor`

```python
get_steric_hindrance_factor(
    model: Model,
    separation: float,
    steric_hindrance_distance: float
) -> float
```

Computes the soft steric velocity reduction factor:

$$f(s) = \frac{1}{2}\left(1 + \tanh\!\frac{s - d}{\lambda}\right)$$

where $s$ is `separation`, $d$ is `steric_hindrance_distance`, and $\lambda$ is `model.model_setup.steric_hindrance_constraint_parameter`. Returns a value in $(0, 1)$ that approaches 0 near the exclusion distance and 1 far from it.

---

## List Manipulation

### `find_and_remove_from_list`

```python
find_and_remove_from_list(
    lst: list[float],
    value: float,
    tolerance: float = 1e-6
) -> None
```

Finds the first element in `lst` within `tolerance` of `value` and removes it in-place. Raises `ValueError` if no matching element is found.

---

## Event Selection

### `select_event_based_on_propensities`

```python
select_event_based_on_propensities(
    rates_vector: list[float],
    p: float
) -> int | None
```

Gillespie event-type selection. Returns the index of the selected event in `rates_vector` given a uniform random number `p ∈ [0, 1)`.

Returns `None` if the total rate is zero (no events possible).

The i-th event is selected when:

$$\sum_{j=0}^{i-1} \frac{r_j}{a_0} \leq p < \sum_{j=0}^{i} \frac{r_j}{a_0}$$

---

## Debugging

### `print_list`

```python
print_list(name: str, lst: list[float]) -> None
```

Prints `name: val1 val2 val3 ...` on a single line. Convenience function for debugging state vectors.
