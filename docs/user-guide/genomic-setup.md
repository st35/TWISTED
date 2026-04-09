# Genomic Setup

`GenomicSetup` stores the static geometric and biological properties of the DNA construct being simulated. It is the first object you create when setting up a simulation.

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
    **kwargs
)
```

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `chromatin_type` | `str` | `'prokaryotic'` or `'eukaryotic'` |
| `gene_names` | `list[str]` | Identifier strings for each gene |
| `TSSes` | `list[float]` | Transcription start site positions in **nm** |
| `gene_lengths` | `list[float]` | Length of each gene in **nm** — RNAP exits when it has traveled this distance from the TSS |
| `gene_directions` | `list[int]` | `+1` for the positive (left→right) strand; `-1` for the negative (right→left) strand |
| `RNAP_on_rates` | `list[float]` | Basal RNAP recruitment rate for each gene in s⁻¹ |
| `promoter_mode` | `str` | `'constitutive'` or `'non-constitutive'` |
| `buffer_length` | `float` | Extra DNA added beyond the last gene end (nm); extends the right clamp position |

All list parameters must have the same length.

### Keyword Arguments

| Keyword | Required when | Type | Description |
|---------|--------------|------|-------------|
| `TF_on_off_rates` | `promoter_mode == 'non-constitutive'` | `list[tuple[float, float]]` | List of `(k_on, k_off)` pairs for each gene's transcription factor |

### Eukaryotic Keyword Arguments

These keywords are used only when `chromatin_type == 'eukaryotic'`:

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `per_nucleosome_DNA_length` | `float` | 147 (bp) | DNA wrapped per nucleosome (passed in bp; converted to nm internally) |
| `nucleosome_linker_length` | `float` | 30 (bp) | Linker DNA between nucleosomes (passed in bp; converted to nm internally) |
| `nucleosomes_are_steric_barriers_to_RNAPs` | `bool` | `True` | Whether nucleosomes block RNAP passage |
| `nucleosome_count` | `int` | auto | Explicit nucleosome count; if omitted, computed by tiling the domain |
| `nucleosome_on_rate_func` | `callable or None` | `None` | Optional function `(segment_length, segment_sigma) → float` that modulates the nucleosome binding rate. Passed as `on_rate_func` to the auto-created nucleosome `BindingProtein` (see [BindingProtein](../api/model-setup.md#bindingprotein)) |
| `nucleosome_off_rate_func` | `callable or None` | `None` | Optional function `(segment_length, segment_sigma) → float` that modulates the nucleosome unbinding rate. Passed as `off_rate_func` to the auto-created nucleosome `BindingProtein` |
| `nucleosomes_can_be_displaced_at_TSS_by_RNAP` | `bool` | `False` | Whether nucleosomes blocking a TSS can be displaced (removed) by an incoming RNAP. Forwarded as `can_be_displaced_at_TSS_by_RNAP` to the auto-created nucleosome `BindingProtein` |

---

## Unit Conventions

**All spatial quantities are in nanometres (nm).** The conversion is:

```
1 bp = 0.34 nm
```

When reading gene coordinates from experimental data (usually in base pairs), multiply by `0.34` before passing to `GenomicSetup`. The helper function `read_genes_information` (see [utilities](../api/utilities.md)) performs this conversion automatically when reading from a file.

---

## DNA Boundaries

The simulated DNA runs from `clamp_left = 0.0 nm` to a `clamp_right` position computed automatically:

- For a positive-strand gene: `clamp_right = TSS + gene_length + buffer_length`
- For a negative-strand gene: `clamp_right = TSS + buffer_length`

The torsional boundary condition at each end is set by `ModelSetup.clamps_status`. A `'clamped'` end fixes the DNA angle (no twist can escape past the boundary). A `'free'` end allows twist to dissipate, so $\dot{Lk}$ for the terminal segment is driven by RNAP displacement rather than angular velocity transfer — see [`get_segment_Lk_dynamics`](../api/biol-methods.md#get_segment_lk_dynamics).

---

## Promoter Modes

### Constitutive

All promoters start in the `ON` state (`promoter_status = 1`) and remain there. RNAP recruitment proceeds at rate `RNAP_on_rates[i]` whenever the TSS is unoccupied.

### Non-constitutive

!!! warning "Not yet implemented"
    The non-constitutive promoter mode is defined in `GenomicSetup` but stochastic promoter toggling is not yet implemented in the simulation loop. Promoters initialised in the `OFF` state will remain `OFF` for the entire simulation, and vice versa. Full support will be added in a future release.

When implemented, promoter state will toggle stochastically between `ON` (1) and `OFF` (0) with rates provided by `TF_on_off_rates`. RNAP will be recruited only when the promoter is `ON`.

```python
genomic_setup = GenomicSetup(
    chromatin_type='prokaryotic',
    gene_names=['geneA'],
    TSSes=[340.0],
    gene_lengths=[3400.0],
    gene_directions=[1],
    RNAP_on_rates=[0.05],
    promoter_mode='non-constitutive',
    buffer_length=3400.0,
    TF_on_off_rates=[(0.1, 0.05)],   # (k_on, k_off) in s⁻¹
)
```

---

## Multi-gene Example

```python
genomic_setup = GenomicSetup(
    chromatin_type='prokaryotic',
    gene_names=['geneA', 'geneB'],
    TSSes=[340.0, 4420.0],           # gene B starts 12,000 bp after gene A's TSS
    gene_lengths=[3400.0, 3400.0],   # both 10,000 bp
    gene_directions=[1, 1],
    RNAP_on_rates=[0.02, 0.02],
    promoter_mode='constitutive',
    buffer_length=4420.0,            # must extend beyond all genes
)
```

!!! warning "Multi-gene constraints"
    The right boundary (`clamp_right`) is computed from the first gene listed only. When simulating multiple genes, either list the gene that defines the rightmost extent of the DNA first, or increase `buffer_length` so that `clamp_right` extends beyond all genes.

---

## Utility Method

### `print_genomic_setup()`

Prints a formatted table of all gene properties to stdout. Useful for verifying setup before long simulations.

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

## Loading from File

Use `utilities.construct_genomic_setup` to build a `GenomicSetup` from a tab-delimited gene file:

```
# gene_file.tsv (columns: name, TSS_bp, length_bp, direction, RNAP_on_rate)
geneA   1000    10000   1   0.02
```

```python
from utilities import construct_genomic_setup

genomic_setup = construct_genomic_setup(
    filename='gene_file.tsv',
    chromatin_type='prokaryotic',
)
```

See [utilities API](../api/utilities.md) for full details.
