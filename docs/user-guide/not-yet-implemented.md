# Not yet implemented

Two features are accepted by the API surface but **not active in the simulation loop**. They are documented in a single location to avoid duplication elsewhere.

---

## 1. `'topoisomerase_based'` supercoiling-relaxation mode

```python
ModelSetup(supercoiling_relaxation_dynamics_mode='topoisomerase_based', ...)
```

The `ModelSetup` constructor recognises this mode value but immediately raises:

```
NotImplementedError: supercoiling_relaxation_dynamics_mode "topoisomerase_based" is not yet implemented.
```

The intended semantics, when implemented, are:

- An explicit copy number of TOP1 and TOP2 enzymes that bind and unbind segments stochastically.
- While bound, each enzyme contributes the per-molecule continuous `dLk/dt` term defined in [`get_TOP1_effect_on_Lk_dynamics`](../api/biol-methods.md#get_top1_effect_on_lk_dynamics) and [`get_TOP2_effect_on_Lk_dynamics`](../api/biol-methods.md#get_top2_effect_on_lk_dynamics).
- Steric exclusion between bound enzymes and other DNA-bound species.

**Use `'topoisomerase_approximated'` as a substitute.** That mode captures the most important biological effect (the TOP1 vs TOP2 distinction, gated by writhe) without the per-enzyme bookkeeping.

---

## 2. `'non-constitutive'` promoter mode

```python
GenomicSetup(..., promoter_mode='non-constitutive', ...)
```

Both `GenomicSetup.__init__` and the helper `construct_genomic_setup` raise:

```
NotImplementedError: Promoter mode "non-constitutive" is not yet implemented.
```

The intended semantics are per-gene transcription-factor on/off kinetics, with the Gillespie loop adding promoter on/off events at user-supplied `TF_on_off_rates`. None of this is wired up at present.

Until promoter switching is implemented, use `promoter_mode='constitutive'` (which initialises every promoter to ON) and emulate transcription-factor regulation by either:

- scaling `RNAP_on_rates` to reflect a steady-state ON probability, or
- wrapping `simulate_dynamics` in an outer loop that toggles `model.promoter_status[i]` between calls.

---

For updates on either feature, see the [GitHub issue tracker](https://github.com/st35/TWISTED/issues).
