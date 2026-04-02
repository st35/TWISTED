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

### Methods

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
        between_RNAPs_steric_effect_cutoff: float = 15.0,
        RNAP_TOPO_steric_effect_cutoff: float = 15.0,
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

## `Model`

Container for the complete dynamic state of a simulation.

```python
class Model:
    def __init__(
        self,
        genomic_setup: GenomicSetup,
        model_setup: ModelSetup
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
