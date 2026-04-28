# Not yet implemented

One feature is accepted by the API surface but **not active in the simulation loop**. It is documented here to avoid duplication elsewhere.

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

For updates on this feature, see the [GitHub issue tracker](https://github.com/st35/TWISTED/issues).
