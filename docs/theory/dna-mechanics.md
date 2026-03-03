# DNA Mechanics

This page describes the physical model underlying TWISTED: supercoiling mechanics, torque-twist relations, RNAP velocity, and topoisomerase activity.

---

## Coordinate System

The DNA is modelled as a one-dimensional elastic rod spanning positions $[0, L]$ nm. RNA polymerases (RNAPs) divide the molecule into segments. Each segment has its own linking number $Lk$, supercoiling density $\sigma$, and torque $\tau$.

Positions increase left-to-right; RNAPs on the positive strand move in the $+x$ direction, and those on the negative strand move in the $-x$ direction.

---

## Supercoiling Density

For a segment of length $L_{\text{seg}}$ with relaxed linking number $Lk_0 = L_{\text{seg}} / h_{\text{dna}}$ where $h_{\text{dna}} = 2\pi / \omega_0$ is the helical repeat:

$$\sigma = \frac{Lk - Lk_0}{Lk_0}$$

Positive $\sigma$ corresponds to overwound (positively supercoiled) DNA; negative $\sigma$ to underwound (negatively supercoiled) DNA.

---

## Prokaryotic Torque–Twist Relation

TWISTED implements the **five-state torque model** for prokaryotic DNA (Marko 2007; Brutzer et al. 2010; Sevier & Levine 2017). The DNA can be in one of five states depending on $\sigma$:

| State | Code | Condition | Torque |
|-------|------|-----------|--------|
| Twisted–melted | 0 | $\sigma \leq \sigma_m$ | $\tau = (c_m / \omega_0)(\sigma - \sigma_0)$ |
| Melted | 1 | $\sigma_m < \sigma \leq \sigma_{sm}$ | $\tau = (c_m / \omega_0)(\sigma_m - \sigma_0)$ |
| Twisted | 2 | $\sigma_{sm} < \sigma \leq \sigma_s$ | $\tau = (c_s / \omega_0)\sigma$ |
| Positive plectoneme | 5 | $\sigma_s < \sigma \leq \sigma_p$ | $\tau = \tau_{\max}$ (torque plateau) |
| Twisted plectoneme | 6 | $\sigma > \sigma_p$ | $\tau = (p / \omega_0)\sigma$ (capped at 40 pN·nm) |

**Material parameters** (fixed):

| Symbol | Value | Description |
|--------|-------|-------------|
| $A$ | 50 nm | Bending persistence length |
| $C$ | 95 nm | Twist persistence length |
| $P$ | 24 nm | Plectoneme supercoil pitch |
| $A_m$ | 4 nm | Bending persistence length (melted) |
| $C_m$ | 1.75 nm | Twist persistence length (melted) |
| $e_m$ | $6\,k_BT$ | Denaturation free energy per base |
| $\sigma_0$ | −1 | Reference supercoiling for melting |

The effective twist stiffness at force $F$:

$$c_s = c \left(1 - \frac{C}{4A}\sqrt{\frac{k_BT}{AF}}\right), \quad c = k_BT C \omega_0^2$$

The plectoneme-formation threshold:

$$\sigma_s = \frac{\sqrt{2pg/(1 - p/c_s)}}{c_s}, \quad g = F - \sqrt{k_BT F / A}$$

**Finite-size correction** (enabled by default for short segments):

$$\sigma_s^{\text{eff}} = \sigma_s \left(1 + \left(\frac{L_0}{L_{\text{seg}}}\right)^2\right)$$

where $L_0 = 340\,\text{nm}$ (1000 bp) by default.

---

## RNAP Velocity

Translational velocity of RNAP $i$ (on a positive-strand gene):

$$v_i = \frac{v_0}{2}\left(1 - \tanh\!\frac{\tau_f - \tau_b}{\tau_c}\right)$$

where $\tau_f$ is the torque in the segment ahead (right) and $\tau_b$ behind (left). For negative-strand genes the sign is flipped.

**Stalling conditions:**

- The RNAP is set to $v = 0$ if the gap to the next RNAP ahead is less than `between_RNAPs_steric_effect_cutoff`.
- In `topoisomerase_based` mode, the RNAP is also stalled when a bound topoisomerase is within `RNAP_TOPO_steric_effect_cutoff` nm ahead.

---

## RNAP Angular Velocity

The rotational degree of freedom couples supercoiling injection to translational motion. For RNAP $i$ at distance $x_i$ from its TSS:

$$\dot{\theta}_i = \frac{\omega_0 \dot{x}_i \chi}{\chi + \eta x_i^\alpha} + \frac{\tau_f - \tau_b}{\chi + \eta x_i^\alpha}$$

The first term reflects the natural rotation imposed by helical translation; the second is the torsional relaxation driven by the torque imbalance.

---

## Linking Number Dynamics

The rate of change of linking number for segment $j$ between two consecutive RNAPs:

$$\dot{Lk}_j = \frac{1}{2\pi}\left(\dot{\theta}_{j-1} - \dot{\theta}_j\right)$$

where $\dot{\theta}_{j-1}$ is the angular velocity of the RNAP on the right (front) side and $\dot{\theta}_j$ is that of the RNAP on the left (back) side. For the rightmost and leftmost segments the boundary terms are zero (clamped ends).

---

## Topoisomerase Activity (`topoisomerase_based` mode)

### Type I Topoisomerase (TOP1)

TOP1 relaxes torsional stress without changing writhe (plectoneme fraction $\Phi_w = 0$ required). The Lk relaxation rate contributed by $n_1$ bound TOP1 molecules on a segment:

$$\dot{Lk}^{\mathrm{TOP1}} = \begin{cases}
-n_1 \, k_0 \, e^{\theta \beta \tau}(1 - e^{-2\pi\beta\tau}) & \sigma > 0 \\
+n_1 \, k_0 \, e^{-\theta \beta \tau}(1 - e^{-2\pi\beta\tau}) & \sigma < 0
\end{cases}$$

where $\beta = 1/k_BT$ and the sign convention ensures relaxation toward $\sigma = 0$.

### Type II Topoisomerase (TOP2)

TOP2 relaxes positive plectonemic supercoiling (segment writhe $\Phi_w > 0$) by passing one DNA double helix through another, reducing writhe by 2:

$$\dot{Lk}^{\mathrm{TOP2}} = \mp n_2 \, V_0 \frac{Wr}{Wr + k_{12}}$$

where $Wr = |\sigma| \Phi_w$ is the writhe and the sign matches the sign of $\sigma$.

---

## Gillespie Event Selection

Between ODE integration intervals, one discrete event is selected. The waiting time $\Delta t$ until the next event satisfies:

$$\int_{t}^{t + \Delta t} a_0(t') \, dt' = \ln(1/p_0), \quad p_0 \sim \mathcal{U}(0,1)$$

where $a_0(t') = \sum_i r_i(t')$ is the total propensity. This integral is tracked via the auxiliary state variable (last element of the state vector) that is accumulated during ODE integration.

The event type $\mu$ is then selected from the normalised rates:

$$P(\text{event} = \mu) = \frac{r_\mu}{a_0}$$
