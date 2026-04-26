# Genomic setup

`GenomicSetup` describes the *static* geometry of the DNA construct: total length, gene count, TSS positions, transcription direction per gene, and (in eukaryotic mode) the nucleosome geometry. It is the first object constructed in a typical workflow.

```python
from model_setup import GenomicSetup

genomic_setup = GenomicSetup(
    chromatin_type='prokaryotic',
    gene_names=['geneA'],
    TSSes=[340.0],
    gene_lengths=[3400.0],
    gene_directions=[1],
    RNAP_on_rates=[0.02],
    promoter_mode='constitutive',
    buffer_length=3400.0,
)
```

---

## Constructor

```python
GenomicSetup(
    chromatin_type: str,
    gene_names: list[str],
    TSSes: list[float],
    gene_lengths: list[float],
    gene_directions: list[int],
    RNAP_on_rates: list[float],
    promoter_mode: str,
    buffer_length: float,
    **kwargs,
)
```

### Required positional/keyword arguments

| Name | Type | Meaning |
|------|------|--------|
| `chromatin_type` | `'prokaryotic' \| 'eukaryotic'` | Selects the torque law and (for eukaryotic) auto-populates a nucleosome `BindingProtein` |
| `gene_names` | `list[str]` | Identifier strings, one per gene |
| `TSSes` | `list[float]` | Transcription start sites in **nm** (not bp) |
| `gene_lengths` | `list[float]` | Distance an RNAP must travel from the TSS in order to count as completed (nm) |
| `gene_directions` | `list[int]` | `+1` (transcribed left → right) or `−1` (transcribed right → left); other values raise |
| `RNAP_on_rates` | `list[float]` | Recruitment rate (s⁻¹) used when the promoter is ON |
| `promoter_mode` | `'constitutive' \| 'non-constitutive'` | See [Promoter modes](#promoter-modes) below. Only `'constitutive'` is currently supported |
| `buffer_length` | `float` | Extra DNA past the gene on the right end (nm); contributes to `clamp_right` |

All five list arguments must have the same length, equal to the number of genes.

### Eukaryotic keyword arguments

These keywords are accepted only when `chromatin_type == 'eukaryotic'`. All inputs in **bp** are converted to nm internally by `× 0.34`.

| Name | Default | Meaning |
|------|--------|--------|
| `per_nucleosome_DNA_length` | 147 bp | DNA wrapped per nucleosome |
| `nucleosome_linker_length` | 30 bp | Linker DNA between adjacent nucleosomes |
| `nucleosomes_are_steric_barriers_to_RNAPs` | `True` | If `True`, nucleosomes block RNAP elongation through the soft tanh ramp |
| `nucleosome_count` | auto | Override the auto-tiled count returned by `get_total_nucleosome_count()` |
| `nucleosome_on_rate_func` | `None` | Optional callable `(L, σ) → float` multiplying the basal nucleosome on-rate |
| `nucleosome_off_rate_func` | `None` | Optional callable `(L, σ) → float` multiplying the basal nucleosome off-rate |
| `nucleosomes_can_be_displaced_at_TSS_by_RNAP` | `False` | If `True`, an incoming RNAP can evict a blocking nucleosome at a TSS |

These values are forwarded to the auto-created nucleosome `BindingProtein` (which lives at `model.binding_proteins[0]` after `Model` is constructed).

---

## Unit conventions

All spatial quantities passed to `GenomicSetup` are in **nm**. Conversion is `1 bp = 0.34 nm`. To work in bp, multiply explicitly or use the [`construct_genomic_setup`](../api/utilities.md#construct_genomic_setup) helper, which reads a tab-delimited bp-denominated config file and converts on input.

---

## DNA boundaries

The DNA spans from `clamp_left = 0.0` (always) to `clamp_right`, which is computed at construction from gene 0 only:

| Direction of gene 0 | Formula |
|---------------------|---------|
| `+1` | `clamp_right = TSSes[0] + gene_lengths[0] + buffer_length` |
| `−1` | `clamp_right = TSSes[0] + buffer_length` |

Other genes are not consulted. **For multi-gene setups, list the rightmost gene first**, or pad `buffer_length` such that `clamp_right` extends past every gene.

The torsional boundary condition at each end is set on `ModelSetup` via `clamps_status` (`'clamped'` vs `'free'`), not on `GenomicSetup`.

---

## Promoter modes

### Constitutive

Every promoter is initialised to ON (`promoter_status[i] = 1`) and stays ON for the entire run. RNAP recruitment proceeds at `RNAP_on_rates[i]` whenever the TSS is unoccupied.

### Non-constitutive

!!! warning "Not yet implemented"
    Passing `promoter_mode='non-constitutive'` raises `NotImplementedError` from both `GenomicSetup` and the helper `construct_genomic_setup`. See [Not yet implemented](not-yet-implemented.md).

---

## Multi-gene constraints (recap)

```python
genomic_setup = GenomicSetup(
    chromatin_type='prokaryotic',
    gene_names=['geneA', 'geneB'],
    TSSes=[340.0, 4420.0],
    gene_lengths=[3400.0, 3400.0],
    gene_directions=[1, 1],
    RNAP_on_rates=[0.02, 0.02],
    promoter_mode='constitutive',
    buffer_length=4420.0,         # makes clamp_right cover both genes
)
```

| Rule | Reason |
|------|--------|
| All list arguments must have the same length | enforced by `assert` |
| `gene_directions[i]` must be `+1` or `−1` | enforced by `assert` |
| `clamp_right` derives from gene 0 only | list the rightmost gene first or increase `buffer_length` |

---

## Helper methods

### `get_total_nucleosome_count() -> int`

For prokaryotic setups, returns 0. For eukaryotic setups, returns either the explicit `nucleosome_count` if supplied, or the count obtained by tiling the DNA from `clamp_left + linker_length/2` rightward in steps of `per_nucleosome_DNA_length + nucleosome_linker_length` until the next nucleosome would overhang `clamp_right`.

### `print_genomic_setup() -> None`

Prints a formatted summary table to stdout. Useful for sanity-checking prior to long runs.

```
Chromatin type: Prokaryotic
Promoter mode: Constitutive
========================================
Name    TSS (nm)    Length (nm)    Direction    RNAP on-rate (1 / s)
----------------------------------------
geneA   340.0       3400.0         1            0.02
========================================
```

---

## Loading from a tab-delimited file

The companion helper [`construct_genomic_setup`](../api/utilities.md#construct_genomic_setup) reads a bp-denominated five-column file and forwards everything else to the `GenomicSetup` constructor:

```text
# columns:  name   TSS_bp   length_bp   direction   RNAP_on_rate(1/s)
geneA       10000  5300     1           1.0
```

```python
from utilities import construct_genomic_setup

genomic_setup = construct_genomic_setup(
    'genes.config',
    chromatin_type='prokaryotic',
    promoter_mode='constitutive',
    explicit_RNAP_on_rates=[0.05],   # multiplies file rate, optional
)
```

Any further keyword arguments (e.g. `per_nucleosome_DNA_length`, `nucleosomes_can_be_displaced_at_TSS_by_RNAP`) are forwarded as-is.
