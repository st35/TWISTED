# Regulatory network

A **gene regulatory network** couples genes together: the amount of protein produced by one gene modulates the promoter switching rate of another. This lets you build circuits such as toggle switches, feed-forward loops, and repressilators on top of the supercoiling dynamics.

The network is optional. When no network is supplied, every gene's promoter switches independently and the model behaves exactly as described in [Genomic setup](genomic-setup.md).

## Requirements

Regulation acts on the **promoter ON-rate** — the rate at which a promoter switches from OFF to ON. Two conditions must be met for a network to have any effect:

- `promoter_mode='non-constitutive'` on `GenomicSetup`. Constitutive promoters are always ON, so there is no switching rate to modulate.
- `mRNA_dynamics_mode=1` on `ModelSetup`. Protein concentrations are derived from mRNA counts, which are only tracked in this mode. The protein production and degradation rates (`protein_production_rate`, `protein_degradation_rate`) are also read from `ModelSetup` (see [Model parameters](model-setup.md)).

If either condition is not met, the network is stored but never influences the dynamics.

## Network file format

The network is loaded from a tab-delimited file with a header row and one regulatory edge per subsequent line:

```text
Source	Target	Type	lambda	Theta	n
gene_A	gene_B	2	0.01	2000.0	4.0
gene_B	gene_A	2	0.01	2000.0	4.0
```

| Column | Meaning |
|--------|---------|
| `Source` | Regulator gene name (must match a gene in the config file) |
| `Target` | Regulated gene name (must match a gene in the config file) |
| `Type` | Regulation type: `1` = activation, `2` = repression (used as a label; see note below) |
| `lambda` | Basal activity floor $\lambda \in [0,1]$; the minimum promoter ON-rate as a fraction of its unregulated value |
| `Theta` | Threshold protein concentration (in units of $k_\text{prod}/k_\text{deg}$ molecules) |
| `n` | Hill coefficient (cooperativity) |

Gene names not present in the config file raise `ValueError`. Each `(Source, Target)` pair defines one edge; a gene may appear as both a source and a target.

## Mathematics

Each regulated gene's promoter ON-rate is scaled by a modulation factor $H$ built from the shifted Hill functions of all its regulators:

$$H = \prod_{\text{regulators } j} \left[\lambda_j + (1-\lambda_j)\frac{1}{1+(p_j^*/\Theta_j)^{n_j}}\right]$$

where the quasi-steady-state protein concentration of regulator $j$ is

$$p_j^* = k_\text{prod} \cdot \text{mRNA}_j / k_\text{deg}.$$

The effective ON-rate for a target gene is then `TF_on_rate × H`, where `TF_on_rate` is the first entry of that gene's `TF_on_off_rates` pair. When $H = 1$ (no active regulators) the ON-rate is unchanged; as a regulator's protein level rises past $\Theta$, $H$ falls toward $\lambda$, throttling the target's ON-rate to a fraction $\lambda$ of its baseline.

!!! note "Activation vs. repression"
    The `Type` column is recorded and shown by `print_genomic_setup`, but the same shifted Hill form is applied to every edge. The function is repressive in shape (decreasing in $p^*$), so the current implementation models repression directly. Use $\lambda$, $\Theta$, and $n$ to tune the strength and sharpness of each interaction.

## Attaching a network

The simplest path is `construct_genomic_setup`, which reads the network file and wires it in automatically via the `regulatory_network_file` keyword:

```python
from utilities import construct_genomic_setup

genomic_setup = construct_genomic_setup(
    'genes.config',
    chromatin_type='prokaryotic',
    promoter_mode='non-constitutive',
    regulatory_network_file='toggle_switch.topo',
)
```

When constructing `GenomicSetup` directly, pass the 4-tuple returned by `read_regulatory_network` as `regulatory_network_information`:

```python
from model_setup import GenomicSetup
from utilities import read_regulatory_network

gene_names = ['gene_A', 'gene_B']
genomic_setup = GenomicSetup(
    chromatin_type='prokaryotic',
    gene_names=gene_names,
    TSSes=[340.0, 4420.0],
    gene_lengths=[3400.0, 3400.0],
    gene_directions=[1, 1],
    RNAP_on_rates=[0.05, 0.05],
    promoter_mode='non-constitutive',
    TF_on_off_rates=[(0.05, 0.01), (0.05, 0.01)],
    buffer_length=4420.0,
    are_multiple_chromosomes_present=False,
    regulatory_network_information=read_regulatory_network(gene_names, 'toggle_switch.topo'),
)
```

For a no-regulation setup, pass `filename=None` to get an empty (all-zero) network:

```python
regulatory_network_information=read_regulatory_network(gene_names, None)
```

## Inspecting the network

`print_genomic_setup` lists every active edge with its interaction type:

```python
genomic_setup.print_genomic_setup()
# ...
# Source	Target	Interaction Type
# gene_A	gene_B	Repression
# gene_B	gene_A	Repression
```

## See also

- [Tutorial 19: Gene regulatory network: toggle switch](../tutorials.md#19-gene-regulatory-network-toggle-switch) — a complete worked example.
- [Model parameters](model-setup.md) — `protein_production_rate` and `protein_degradation_rate`.
- [`utilities` API reference](../api/utilities.md) — `read_regulatory_network` signature.
