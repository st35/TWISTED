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
) -> int
```

Returns `1` if the TSS at `TSS_position` is sterically blocked, `0` otherwise.

Blocking is triggered by:
- Any active RNAP within `between_RNAPs_steric_effect_cutoff` nm of the TSS.
- Any bound topoisomerase within `RNAP_TOPO_steric_effect_cutoff` nm of the TSS (in `topoisomerase_based` mode only).

### `is_TOPO_binding_blocked`

```python
is_TOPO_binding_blocked(
    model: Model,
    RNAP_gene_index: list[int],
    state_vector: list[float],
    binding_position: float
) -> int
```

Returns `1` if a topoisomerase cannot bind at `binding_position` because an RNAP is within `RNAP_TOPO_steric_effect_cutoff` nm, `0` otherwise.

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
