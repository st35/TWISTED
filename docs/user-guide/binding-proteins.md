# Binding proteins

`BindingProtein` describes any DNA-binding protein species other than RNAP. The same class powers nucleosomes (auto-created in eukaryotic mode) and arbitrary user-defined species (passed to the `Model` constructor via `binding_proteins=[...]`).

```python
from model_setup import BindingProtein

protein = BindingProtein(
    protein_name='ProteinX',
    total_copy_number=10,
    is_steric_barrier_to_RNAPs=False,
    is_topological_barrier=False,
    basal_on_rate=0.001,    # per (s · nm)
    basal_off_rate=0.01,    # per s
)
```

---

## Constructor

```python
BindingProtein(
    protein_name: str,
    total_copy_number: int,
    is_steric_barrier_to_RNAPs: bool,
    is_topological_barrier: bool,
    basal_on_rate: float,
    basal_off_rate: float,
    on_rate_func: Callable | None = None,
    off_rate_func: Callable | None = None,
    is_a_nucleosome: bool = False,
    can_be_displaced_at_TSS_by_RNAP: bool = False,
)
```

| Argument | Meaning |
|----------|--------|
| `protein_name` | Identifier string |
| `total_copy_number` | Total number of molecules of this species (bound + unbound) |
| `is_steric_barrier_to_RNAPs` | If `True`, this protein blocks RNAP elongation through the soft tanh ramp and contributes to the TSS steric-hindrance check during recruitment |
| `is_topological_barrier` | If `True`, each bound molecule splits one DNA segment into two with independent linking numbers. **Implies `is_steric_barrier_to_RNAPs=True`**; the constructor raises otherwise |
| `basal_on_rate` | Per-(s·nm) intrinsic on-rate (see semantics below) |
| `basal_off_rate` | Per-(s) intrinsic off-rate, per bound molecule |
| `on_rate_func` | Optional `(segment_length, segment_sigma) → float` multiplier on the basal on-rate |
| `off_rate_func` | Optional `(segment_length, segment_sigma) → float` multiplier on the basal off-rate |
| `is_a_nucleosome` | Set `True` only by the framework when auto-creating the nucleosome species in eukaryotic mode; nucleosomes use a different exclusion distance in steric checks |
| `can_be_displaced_at_TSS_by_RNAP` | If `True`, an incoming RNAP can evict a bound molecule of this species at the TSS even though it is otherwise a steric barrier |

---

## Rate semantics

The wrapping that turns the basal rate into a Gillespie propensity is added by the constructor itself. The actual per-segment expressions used inside the simulator are:

### Per-segment on-propensity

```
n_unbound × basal_on_rate × segment_length × user_on_rate_func(L, σ)
```

- `n_unbound = total_copy_number − sum(len(positions))`.
- `user_on_rate_func` defaults to `1` when not supplied.
- The `× segment_length` factor causes longer segments to attract proportionally more binding, so `basal_on_rate` is denominated in **per (s · nm)**.
- The chosen segment is sampled with weight equal to its on-propensity. Within the chosen segment, the position is sampled uniformly.

### Per-bound-molecule off-propensity

```
basal_off_rate × user_off_rate_func(L, σ)
```

- `L` and `σ` are the length and supercoiling density of the segment in which the molecule is bound.
- The constructor wraps any user-supplied `off_rate_func` so that it multiplies `basal_off_rate`. When unset, the off-rate equals `basal_off_rate` per bound molecule.

### Total propensities reported to the Gillespie loop

The simulator sums these per-segment on-rates and per-bound off-rates and exposes a single *aggregated* propensity per species per direction (one for binding, one for unbinding). When the species is selected, segment selection (for binding) and per-molecule selection (for unbinding) is performed by a second random draw weighted by the per-segment / per-molecule propensities.

---

## Custom σ-sensitive kinetics

A common use case is biasing binding toward negatively supercoiled regions, or accelerating unbinding under positive torque:

```python
def fast_on_when_negative(segment_length, segment_sigma):
    return 5.0 if segment_sigma < -0.02 else 1.0

def fast_off_when_positive(segment_length, segment_sigma):
    return 10.0 if segment_sigma > 0.02 else 1.0

sc_sensor = BindingProtein(
    protein_name='SC_sensor',
    total_copy_number=15,
    is_steric_barrier_to_RNAPs=False,
    is_topological_barrier=False,
    basal_on_rate=0.001,
    basal_off_rate=0.01,
    on_rate_func=fast_on_when_negative,
    off_rate_func=fast_off_when_positive,
)
```

Both user functions receive `(segment_length, segment_sigma)`; the underlying call site does not pass any additional arguments.

---

## Steric-barrier behaviour

When `is_steric_barrier_to_RNAPs=True`, the protein contributes to two distinct steric checks:

1. **Elongation steric ramp**: for any approaching RNAP, the velocity is multiplied by the soft tanh factor with exclusion distance
   - `(RNAP_diameter + per_nucleosome_DNA_length + nucleosome_linker_length) / 2` if the protein is a nucleosome,
   - `(RNAP_diameter + generic_binding_protein_diameter) / 2` otherwise.
2. **TSS recruitment block**: when an RNAP recruitment event fires at a TSS, the simulator iterates over all active RNAPs, then over all bound nucleosomes (steric-barrier ones), then over all bound generic steric-barrier proteins. If anything is within the relevant exclusion distance of the TSS, recruitment is blocked. If the *only* obstacle is a protein with `can_be_displaceable_at_TSS_by_RNAP=True`, recruitment succeeds and that protein is removed from `model.binding_proteins_positions`.

A protein with `is_steric_barrier_to_RNAPs=False` is invisible to RNAPs in both senses; it can still bind, unbind, and act as a topological barrier.

---

## Topological barriers

`is_topological_barrier=True` causes each bound molecule to act as an independent boundary in the segmentation of the DNA, on top of the boundaries created by RNAPs. Concretely:

- When the protein **binds**, the simulator inserts a boundary at the binding position. The two new sub-segments inherit the supercoiling density of the parent segment (i.e. their `Lk` are set so that `σ` is preserved on both sides).
- When the protein **unbinds**, the two adjacent segments are merged and their linking numbers are added.
- The bookkeeping is implemented by [`update_Lk_vector_after_RNAP_or_protein_recruitment`](../api/model-dynamics.md#update_lk_vector_after_rnap_or_protein_recruitment) and [`update_Lk_vector_after_protein_unbinding`](../api/model-dynamics.md#update_lk_vector_after_protein_unbinding).

The constructor enforces `is_topological_barrier ⇒ is_steric_barrier_to_RNAPs`. A topological barrier that did not also block RNAPs would be inconsistent with the segmentation algorithm.

---

## Nucleosomes (eukaryotic mode)

When `chromatin_type == 'eukaryotic'`, `Model.__init__` automatically constructs a nucleosome `BindingProtein` with:

- `protein_name='nucleosome'`
- `total_copy_number = genomic_setup.get_total_nucleosome_count()`
- `is_a_nucleosome=True`
- `is_steric_barrier_to_RNAPs = genomic_setup.nucleosomes_are_steric_barriers_to_RNAPs`
- `is_topological_barrier=False`
- `basal_on_rate = 1.2 / (clamp_right − clamp_left)`
- `basal_off_rate = 0.4`
- `on_rate_func = genomic_setup.nucleosome_on_rate_func`
- `off_rate_func = genomic_setup.nucleosome_off_rate_func`
- `can_be_displaced_at_TSS_by_RNAP = genomic_setup.nucleosomes_can_be_displaced_at_TSS_by_RNAP`

This species is inserted at index 0 of `model.binding_proteins`, **before** any user-supplied species. The eukaryotic indexing pattern is therefore:

```python
model = Model(genomic_setup, model_setup, binding_proteins=[my_protein])
# model.binding_proteins[0]  -> auto-created nucleosome
# model.binding_proteins[1]  -> my_protein
```

To skip the auto-created nucleosomes entirely, set `chromatin_type='prokaryotic'` and supply a nucleosome-like `BindingProtein` explicitly.

---

## Inspecting bound state

After (or during) a run, the bound positions are accessible at:

```python
model.binding_proteins_positions[k]    # list of nm positions for species k
```

The list of `BindingProtein` objects is at `model.binding_proteins[k]`. To pretty-print all species:

```python
model.print_model_setup()
```
