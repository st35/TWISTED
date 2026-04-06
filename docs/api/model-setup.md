# `model_setup` — Core Data Classes

**File:** `model_setup.py`
**Imports:** `utilities`

This module defines the four principal data classes used throughout TWISTED.

---

## `GenomicSetup`

Holds the static genomic and geometric properties of the simulated DNA construct.

```python
class GenomicSetup:
    def __init__(
        self,
        chromatin_type: str,
        gene_names: list[str],
        TSSes: list[float],
        gene_lengths: list[float],
        gene_directions: list[int],
        RNAP_on_rates: list[float],
        promoter_mode: str,
        buffer_length: float,
        **kwargs
    ) -> None
```

See [Genomic Setup](../user-guide/genomic-setup.md) for full parameter descriptions and examples.

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `chromatin_type` | `str` | `'prokaryotic'` or `'eukaryotic'` |
| `gene_names` | `list[str]` | Gene identifiers |
| `TSSes` | `list[float]` | TSS positions (nm) |
| `gene_lengths` | `list[float]` | Gene lengths (nm) |
| `gene_directions` | `list[int]` | `+1` or `-1` per gene |
| `RNAP_on_rates` | `list[float]` | Basal recruitment rates (s⁻¹) |
| `promoter_mode` | `str` | `'constitutive'` or `'non-constitutive'` |
| `TF_on_off_rates` | `list[tuple[float,float]]` | TF toggle rates (non-constitutive only) |
| `clamp_left` | `float` | Left DNA boundary (always 0.0 nm) |
| `clamp_right` | `float` | Right DNA boundary (nm) |

Additional attributes present only when `chromatin_type == 'eukaryotic'`:

| Attribute | Type | Description |
|-----------|------|-------------|
| `per_nucleosome_DNA_length` | `float` | DNA wrapped per nucleosome (nm); default 147 bp × 0.34 = 49.98 nm |
| `nucleosome_linker_length` | `float` | Linker DNA between nucleosomes (nm); default 30 bp × 0.34 = 10.2 nm |
| `nucleosomes_are_steric_barriers_to_RNAPs` | `bool` | Whether nucleosomes block RNAP passage; default `True` |
| `explicit_nucleosome_count` | `int or None` | User-specified nucleosome count; if `None`, computed automatically by tiling |

These can be set via `**kwargs` in the constructor (pass values in **bp** for lengths; they are converted to nm internally).

### Methods

#### `get_total_nucleosome_count() → int`

Returns the total number of nucleosomes for the construct. Returns `0` for prokaryotic setups. If `explicit_nucleosome_count` is set, returns that value. Otherwise, tiles nucleosomes across the domain starting from `clamp_left + linker/2`, placing one nucleosome per `per_nucleosome_DNA_length + nucleosome_linker_length` interval.

#### `print_genomic_setup() → None`

Prints a formatted summary table of gene properties to stdout.

---

## `ModelSetup`

Holds all tunable physical and biological parameters.

```python
class ModelSetup:
    def __init__(
        self,
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
        TOPO_diameter: float = 15.0,
        generic_binding_protein_diameter: float = 15.0,
        steric_hindrance_constraint_parameter: float = 2.0,
        clamps_status: tuple[str, str] = ('clamped', 'clamped'),
        finite_size_effect_flag: int = 1,
        supercoiling_relaxation_dynamics_mode: str = 'global_overall',
        mRNA_dynamics_mode: int = 0,
        model_observation_event_rate: float = 0.5,
        **kwargs
    ) -> None
```

See [Model Parameters](../user-guide/model-setup.md) for full parameter descriptions.

### Key Derived Attributes

| Attribute | Description |
|-----------|-------------|
| `h_dna` | Helical repeat `= 2π / w0` (nm/turn) |
| `left_clamp_status` | `1` = torsionally clamped, `0` = free (derived from `clamps_status[0]`) |
| `right_clamp_status` | `1` = torsionally clamped, `0` = free (derived from `clamps_status[1]`) |
| `global_supercoiling_relaxation_rate` | Set from `**kwargs` when mode is `global_overall` or `global_per_segment` |
| `local_supercoiling_relaxation_rates` | `[rate_pos, rate_neg]` when mode is `global_by_type` or `per_segment_by_type` |
| `TOP1_effective_relaxation_rate` | Effective TOP1 rate when mode is `topoisomerase_approximated` |
| `TOP2_effective_relaxation_rate` | Effective TOP2 rate when mode is `topoisomerase_approximated` |
| `topoisomerase_copy_numbers` | `[n_TOP1, n_TOP2]` when mode is `topoisomerase_based` |
| `topoisomerase_on_off_rates` | `[(k_on1, k_off1), (k_on2, k_off2)]` when mode is `topoisomerase_based` |
| `mRNA_degradation_rate` | mRNA degradation rate (s⁻¹) when `mRNA_dynamics_mode=1` |
| `finite_size_effect_length` | Length scale for finite-size corrections (nm) |
| `supercoiling_relaxation_dynamics_modes_with_no_steric_hindrance` | List of mode names that do not model explicit topoisomerase binding |

---

## `BindingProtein`

Represents a type of DNA-binding protein that can bind and unbind from DNA segments during the simulation.

```python
class BindingProtein:
    def __init__(
        self,
        protein_name: str,
        total_copy_number: int,
        is_steric_barrier_to_RNAPs: bool,
        is_topological_barrier: bool,
        basal_on_rate: float,
        basal_off_rate: float,
        on_rate_func: callable = None,
        off_rate_func: callable = None,
        is_a_nucleosome: bool = False
    ) -> None
```

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `protein_name` | `str` | Identifier for this protein type |
| `total_copy_number` | `int` | Total number of molecules of this protein |
| `is_steric_barrier_to_RNAPs` | `bool` | Whether bound proteins block RNAP passage |
| `is_topological_barrier` | `bool` | Whether bound proteins act as barriers to supercoiling diffusion (not yet implemented) |
| `basal_on_rate` | `float` | Basal binding rate (s⁻¹ nm⁻¹); the per-segment on-rate is `basal_on_rate × segment_length` |
| `basal_off_rate` | `float` | Basal unbinding rate (s⁻¹) |
| `on_rate_func` | `callable or None` | Optional function `(segment_length, segment_sigma) → float` that multiplies the basal on-rate × segment_length. Defaults to `1.0` (no modulation) |
| `off_rate_func` | `callable or None` | Optional function `(segment_length, segment_sigma) → float` that multiplies the basal off-rate. Defaults to `1.0` (no modulation) |
| `is_a_nucleosome` | `bool` | If `True`, this protein is treated as a nucleosome with special steric handling (uses `per_nucleosome_DNA_length + nucleosome_linker_length` as its physical extent instead of `generic_binding_protein_diameter`). Default `False` |

The effective per-segment on-rate for unbound proteins is:

$$r_{\text{on}} = n_{\text{unbound}} \times \texttt{basal\_on\_rate} \times L_{\text{segment}} \times f_{\text{on}}(L, \sigma)$$

The effective off-rate for a bound protein on a given segment is:

$$r_{\text{off}} = \texttt{basal\_off\_rate} \times f_{\text{off}}(L, \sigma)$$

!!! note "Steric barriers"
    `is_steric_barrier_to_RNAPs` is enforced: bound proteins with this flag set to `True` will stall nearby RNAPs and block RNAP recruitment at nearby TSSes. Nucleosomes (`is_a_nucleosome=True`) use `per_nucleosome_DNA_length + nucleosome_linker_length` as their physical extent for steric checks; other proteins use `generic_binding_protein_diameter`.

!!! warning "Not yet implemented"
    `is_topological_barrier` is stored but not yet enforced during the simulation. Support will be added in a future release.

---

## `Model`

Container for the complete dynamic state of a simulation.

```python
class Model:
    def __init__(
        self,
        genomic_setup: GenomicSetup,
        model_setup: ModelSetup,
        binding_proteins: list[BindingProtein] = None
    ) -> None
```

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `genomic_setup` | `GenomicSetup` | Reference to the genomic setup |
| `model_setup` | `ModelSetup` | Reference to the model setup |
| `x_dict` | `list[list[float]]` | RNAP positions per gene (nm). Index: `[gene_index][rnap_index]` |
| `theta_dict` | `list[list[float]]` | RNAP angular positions per gene (rad). Same indexing as `x_dict` |
| `Lk` | `list[float]` | Linking number of each DNA segment (length = n_RNAPs + 1). Ordered right-to-left |
| `promoter_status` | `list[int]` | `1` (ON) or `0` (OFF) per gene |
| `mRNA_counts` | `list[int]` | mRNA copy numbers per gene |
| `binding_proteins` | `list[BindingProtein]` | List of binding protein types. For eukaryotic setups, a nucleosome `BindingProtein` is automatically prepended at index 0 |
| `binding_proteins_positions` | `list[list[float]]` | Positions (nm) of bound proteins, indexed `[protein_type][bound_index]` |

Additional attributes present only when `supercoiling_relaxation_dynamics_mode == 'topoisomerase_based'`:

| Attribute | Type | Description |
|-----------|------|-------------|
| `topoisomerase_type` | `list[int]` | `0` = TOP1, `1` = TOP2 for each topoisomerase copy |
| `topoisomerase_positions` | `list[float]` | Current position (nm); `-1.0` if unbound |
| `topoisomerase_segment_indices` | `list[int]` | DNA segment index occupied; `-1` if unbound |
| `topoisomerase_status` | `list[int]` | `0` = unbound, `1` = bound |

### Initial State

- `Lk` is initialised to a single element corresponding to the fully relaxed linking number of the complete DNA.
- All RNAPs start absent (`x_dict` and `theta_dict` are empty lists per gene).
- Constitutive promoters start `ON`; non-constitutive promoters start with a random binary state.
- For eukaryotic setups, a nucleosome `BindingProtein` is automatically created and prepended to `binding_proteins`. Its `total_copy_number` is determined by `GenomicSetup.get_total_nucleosome_count()`, its `basal_on_rate` is normalised by total DNA length (`1.2 / (clamp_right - clamp_left)`), and `is_steric_barrier_to_RNAPs` is set from `GenomicSetup.nucleosomes_are_steric_barriers_to_RNAPs`. User-supplied binding proteins follow at subsequent indices.

!!! warning "Non-constitutive promoters"
    Stochastic promoter toggling is not yet implemented. Non-constitutive promoters retain their initial random state for the entire simulation. Full support is planned for a future release.

---

## `SimulationSetupAndState`

Controls simulation termination and accumulates results.

```python
class SimulationSetupAndState:
    def __init__(
        self,
        genomic_setup: GenomicSetup,
        simulation_end_mode: int,
        simulation_end_criterion: Union[float, list[int]],
        integration_method: str = 'RK23',
        integration_time_resolution: float = 0.1,
        RNAP_alive_status_check_interval: float = 1.0,
        max_RNAPs_to_recruit: list[int] = None
    ) -> None
```

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `simulation_end_mode` | `int` | `0` = time-based; `1` = event-based |
| `simulation_end_time` | `float` | Target simulation time (mode 0 only) |
| `simulation_end_event_counts` | `list[int]` | Per-gene RNAP completion targets (mode 1 only) |
| `integration_method` | `str` | ODE solver method passed to `solve_ivp`. One of `RK23`, `RK45`, `DOP853`, `Radau`, `BDF`, `LSODA`. Solvers other than `RK23` may crash due to non-physical intermediate states violating steric constraints |
| `integration_time_resolution` | `float` | ODE evaluation time step (s) |
| `RNAP_alive_status_check_interval` | `float` | RNAP status check interval (s) |
| `max_RNAPs_to_recruit` | `list[int] or None` | Cap on total recruits per gene |
| `RNAPs_finished_transcription` | `list[int]` | Count of completed RNAPs per gene |
| `RNAPs_exit_positions` | `list[list[float]]` | Exit positions (nm) per gene |
| `RNAP_recruitment_times` | `list[list[float]]` | Recruitment times (s) per gene |
| `RNAP_exit_times` | `list[list[float]]` | Completion times (s) per gene |
| `curr_simulation_time` | `float` | Current simulation clock (s) |
| `simulation_completed` | `bool` | True when termination criterion is met |

### Methods

#### `calculate_RNAP_transcription_rates(model: Model) → list[list[float]]`

Returns the effective transcription rate in **bp/s** for each RNAP that completed transcription, grouped by gene.

```python
rates = sim.calculate_RNAP_transcription_rates(model)
# rates[gene_index][rnap_index] → float (bp/s)
```

The rate is computed as the distance from the TSS to the exit position (converted to bp) divided by the elapsed time from recruitment to exit.
