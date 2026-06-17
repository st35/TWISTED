# `model_setup`

Data classes the user constructs.

---

## class `GenomicSetup`

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

Holds the genomic layout: chromatin type, gene table, promoter mode, and DNA extents. The right-hand clamp position is computed automatically as `TSSes[0] + gene_lengths[0] + buffer_length` for a `+1` first gene, or `TSSes[0] + buffer_length` for a `−1` first gene.

`chromatin_type` ∈ `{'prokaryotic', 'eukaryotic'}`. `promoter_mode` ∈ `{'constitutive', 'non-constitutive'}`. For `'non-constitutive'`, the kwarg `TF_on_off_rates: list[tuple[float, float]]` is required — one `(TF_on_rate, TF_off_rate)` pair per gene (s⁻¹). For `'constitutive'`, `TF_on_off_rates` defaults to `[(0.0, 0.0), ...]` and is unused.

Eukaryotic kwargs: `per_nucleosome_DNA_length` (bp, default 147), `nucleosome_linker_length` (bp, default 30), `nucleosomes_are_steric_barriers_to_RNAPs` (bool, default `True`), `nucleosomes_can_be_displaced_at_TSS_by_RNAP` (bool, default `False`), `nucleosome_count` (int, override the auto-computed count), `nucleosome_on_rate_func`, `nucleosome_off_rate_func` (callables `(L, σ) → factor`).

### `get_total_nucleosome_count() -> int`

Returns 0 for prokaryotic, otherwise either the user-supplied `nucleosome_count` or the count obtained by densely tiling `clamp_left → clamp_right` with `per_nucleosome_DNA_length + nucleosome_linker_length`.

### `print_genomic_setup() -> None`

Pretty-print the gene table.

See also: [User guide → Genomic setup](../user-guide/genomic-setup.md).

---

## class `ModelSetup`

```python
ModelSetup(
    w0: float = 1.85,
    chi: float = 0.05,
    eta: float = 0.0005,
    alpha: float = 1.5,
    v0: float = 20.0,
    tau_c: float = 12.0,
    force: float = 1.0,
    kBT: float = 4.1,
    TOP1_k0: float = 11.0,
    TOP1_theta: float = 0.25,
    TOP2_V0: float = 2.6,
    TOP2_k12: float = 2.0,
    RNAP_diameter: float = 15.0,
    generic_binding_protein_diameter: float = 15.0,
    steric_hindrance_constraint_parameter: float = 2.0,
    clamps_status: tuple[str, str] = ('clamped', 'clamped'),
    finite_size_effect_flag: int = 1,
    supercoiling_relaxation_dynamics_mode: str = 'global_overall',
    mRNA_dynamics_mode: int = 0,
    model_observation_event_rate: float = 0.5,
    **kwargs,
)
```

Holds DNA mechanics constants, RNAP coupling constants, topoisomerase rate constants, geometric constants, boundary-condition flags, and the global "what kind of supercoiling relaxation" switch.

`h_dna = 2π / w0` is computed in the constructor.

Required kwargs depending on mode:

| Selected mode | Required kwargs |
|--------------|----------------|
| `'global_overall'`, `'global_per_segment'` | `global_supercoiling_relaxation_rate` |
| `'global_by_type'`, `'per_segment_by_type'` | `local_supercoiling_relaxation_rates` (list `[pos, neg]`) |
| `'topoisomerase_approximated'` | `TOP1_effective_relaxation_rate`, `TOP2_effective_relaxation_rate` |
| `'topoisomerase_based'` | raises `NotImplementedError` |

If `mRNA_dynamics_mode == 1`, `mRNA_degradation_rate` is required. If `finite_size_effect_flag == 1`, `finite_size_effect_length` defaults to 340 nm.

See also: [User guide → Model parameters](../user-guide/model-setup.md).

---

## class `BindingProtein`

```python
BindingProtein(
    protein_name: str,
    total_copy_number: int,
    is_steric_barrier_to_RNAPs: bool,
    is_topological_barrier: bool,
    basal_on_rate: float,
    basal_off_rate: float,
    on_rate_func: callable = None,
    off_rate_func: callable = None,
    is_a_nucleosome: bool = False,
    can_be_displaced_at_TSS_by_RNAP: bool = False,
)
```

Description of one DNA-binding species. The aggregated per-segment on-rate is

```
n_unbound × basal_on_rate × segment_length × on_rate_func(L, σ)
```

(if `on_rate_func` is `None`, it defaults to the constant `1`, giving the documented `basal_on_rate × L` form).

The aggregated per-molecule off-rate is

```
basal_off_rate × off_rate_func(L, σ)
```

A `is_topological_barrier=True` species must also be `is_steric_barrier_to_RNAPs=True` (enforced at construction).

See also: [User guide → Binding proteins](../user-guide/binding-proteins.md).

---

## class `Model`

```python
Model(
    genomic_setup: GenomicSetup,
    model_setup: ModelSetup,
    binding_proteins: list[BindingProtein] = None,
)
```

Combines `GenomicSetup` and `ModelSetup` and initialises mutable state: `x_dict`, `theta_dict`, `Lk` (single segment spanning the whole DNA), `promoter_status`, `mRNA_counts`, `binding_proteins`, `binding_proteins_positions`. Eukaryotic models have a synthetic `nucleosome` `BindingProtein` prepended automatically.

### `print_model_setup() -> None`

Print the genomic setup followed by the binding-protein table.

---

## class `SimulationSetupAndState`

```python
SimulationSetupAndState(
    simulation_end_mode: int,
    simulation_end_criterion: float | list[int],
    integration_method: str = 'RK23',
    integration_time_resolution: float = 0.1,
    RNAP_alive_status_check_interval: float = 1.0,
    max_RNAPs_to_recruit: list[int] = None,
)
```

Termination policy + integrator settings + result accumulators in one object. The constructor no longer takes a `GenomicSetup`; per-gene state is allocated later by `setup_simulation_state` (called automatically by `simulate_dynamics`).

State attributes (filled during `simulate_dynamics`):

- `RNAPs_finished_transcription[i]`: int
- `RNAP_recruitment_times[i]`, `RNAP_exit_times[i]`, `RNAPs_exit_positions[i]`: lists
- `curr_simulation_time`: float
- `last_event_index`: int — index of the most recently selected Gillespie event; `-1` before any event has occurred
- `last_event_type`: str | None — human-readable label of the most recently dispatched event type; `None` before any event has occurred (see [dispatch table](simulate-dynamics.md#event-dispatch-table) for possible values)
- `simulation_completed`: bool

### `setup_simulation_state(genomic_setup) -> None`

Allocates the per-gene result lists (`RNAPs_finished_transcription`, `RNAPs_exit_positions`, `RNAP_recruitment_times`, `RNAP_exit_times`) sized to `len(genomic_setup.gene_names)`, and validates that `simulation_end_criterion` and `max_RNAPs_to_recruit` (if provided) have the matching length and that `simulation_end_event_counts[i] <= max_RNAPs_to_recruit[i]` in event-count mode. Called automatically at the start of `simulate_dynamics`; users normally do not call it directly.

### `calculate_RNAP_transcription_rates(model) -> list[list[float]]`

Per-gene per-RNAP transcription rate (bp/s), computed only for RNAPs that completed transcription. Asserts that the distance covered is at least `gene_length`.

See also: [User guide → Running simulations](../user-guide/simulation.md).
