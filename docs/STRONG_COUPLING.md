# Strongly Coupled Free-Flight UVLM-MuJoCo Solver: Mathematical Framework

This document specifies the implicit (strong) coupling between the free-flight UVLM aerodynamics and the MuJoCo rigid-body dynamics. Where explicit (loose) coupling advances each subsystem once per step on the other's lagged, previous-step data, strong coupling drives the two to mutual consistency within the step. That within-step solve is a fixed-point iteration accelerated with Aitken's method.

---

## 1. Tagging conventions

1. Snapshot quantities carry the subscript of their snapshot. Snapshots are instants $t_n$, indexed from zero (e.g., $t_0 = 0.0\,\mathrm{s}$, $t_1 = 0.1\,\mathrm{s}$), with $\Delta t = t_{n+1} - t_n$.
2. Interval quantities (loads applied over $\mathopen{}\mathclose{\left(t_n, t_{n+1}\right]}$, state increments, the newly shed wake row) carry the subscript of the snapshot that ends the interval.
3. Trial values inside the solve for snapshot $n{+}1$ carry a parenthesized iteration superscript $\mathopen{}\mathclose{\left(k\right)}$ in addition to the subscript $n{+}1$. The iteration axis $k$ is orthogonal to the time axis $n$.
4. Constants are unsubscripted by construction: $\Delta t$, the diagonal weighting matrix $\mathbf{D}$ (Section 7), the initial relaxation factor $\lambda_{\text{init}}$, and the relative and absolute convergence tolerances $\varepsilon_{\text{rel}}$ and $\varepsilon_{\text{abs}}$.
5. Functions $\operatorname{G}$ (bound geometry), $\operatorname{S}$ (structural propagator), $\operatorname{P}$ (wake propagator), $\operatorname{C}$ (circulation solve), $\operatorname{L}$ (load evaluation), and $\operatorname{A}$ (composite aerodynamic evaluation) are time-generic: the snapshot is carried entirely by their tagged arguments, and the prescribed flapping kinematics are evaluated at the tag of the state argument. $\operatorname{G}$ maps a body state and time to the bound ring vortex positions (Section 3), so the bound geometry is a derived quantity rather than independent state. $\operatorname{P}$ spans the interval: it reads the snapshot-$n$ bound geometry, strengths, and wake for its induced transport and shed strengths, the snapshot-$n$ state $\mathbf{x}_n$ for the surface image in that transport (which drops out in free air), the trial bound geometry $\mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$ for its shed leading edges, and the trial interval-end state $\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$ for its frame terms (Section 3).

## 2. State vector and snapshot data

The body state is the 12-vector:

$$\mathbf{x}_n = \mathopen{}\mathclose{\left(\mathbf{p}_n, \boldsymbol{\theta}_n, \mathbf{v}_n, \boldsymbol{\omega}_n\right)} \in \mathbb{R}^{12}$$

Its blocks are the project-tagged quantities:

| Block                   | Tagged quantity           | Meaning                                                                                                    |
|-------------------------|---------------------------|------------------------------------------------------------------------------------------------------------|
| $\mathbf{p}_n$          | `CgP1_E_Eo`               | Position of the first `Airplane`'s CG, in Earth axes, relative to the Earth origin (m)                     |
| $\boldsymbol{\theta}_n$ | `anglesRad_E_to_BP1_izyx` | Intrinsic zyx Euler angles from Earth axes to the first `Airplane`'s body axes (rad)                       |
| $\mathbf{v}_n$          | `vCgP1_E__E`              | Velocity of the first `Airplane`'s CG, in Earth axes, observed from the Earth frame (m/s)                  |
| $\boldsymbol{\omega}_n$ | `omegasP1Rad_E__E`        | Angular velocity of the first `Airplane`'s body axes, in Earth axes, observed from the Earth frame (rad/s) |

All four blocks are components in a single iterate-independent axes system (Earth axes), with angular quantities in radians, so the residuals, updates, and dot products of Section 5 are plain vector arithmetic. The state is a conversion of MuJoCo's internal representation (position plus quaternion in `qpos`, velocities in `qvel`, where MuJoCo world axes equal Earth axes), not MuJoCo's variables directly, and the conversion runs one way: every $\operatorname{S}$ evaluation ends with an extraction into these tags, while restores go through the saved native MuJoCo snapshot. The `OperatingPoint` seam inside $\operatorname{A}$ converts per evaluation: degrees for the angle block, and a rotation into body axes plus degrees per second for `omegas_BP1__E`.

Known before solving for snapshot $n{+}1$:

| Quantity                | Meaning                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|-------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| $\mathbf{x}_n$          | Body state at $t_n$                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| $\mathbf{B}_n$          | Bound ring vortex positions at $t_n$ (GP1 axes, the first `Airplane`'s geometry axes), $\mathbf{B}_n = \operatorname{G}\mathopen{}\mathclose{\left(\mathbf{x}_n, t_n\right)}$: the prescribed flapping geometry set by $t_n$, plus the trailing-edge back-leg offset set by the freestream (hence by the velocity and attitude blocks of $\mathbf{x}_n$, Section 3). The initial core radius is a per-wing geometric constant (as for the wake). Bound vortices carry zero age. |
| $\boldsymbol{\Gamma}_n$ | Bound ring vortex strengths at $t_n$: the column vector produced by the circulation solve $\operatorname{C}$.                                                                                                                                                                                                                                                                                                                                           |
| $W_n$                   | Wake state object at $t_n$: positions (GP1 axes) and strengths of all rows shed through interval $n$. The viscous core-growth model also consumes each wake vortex's age and initial core radius, but neither is independent state: the age is rederived from the vortex's position in the wake, and the initial core radius is a per-wing geometric constant.                                                                                          |
| $\boldsymbol{\ell}_n$   | Official snapshot-$n$ load: $\boldsymbol{\ell}_n = \operatorname{L}\mathopen{}\mathclose{\left(\mathbf{x}_n, \mathbf{B}_n, \boldsymbol{\Gamma}_n, \boldsymbol{\Gamma}_{n-1}, W_n\right)}$ with $\boldsymbol{\Gamma}_n = \operatorname{C}\mathopen{}\mathclose{\left(\mathbf{x}_n, \mathbf{B}_n, W_n\right)}$ (functions in Section 3), equivalently $\boldsymbol{\ell}_n = \operatorname{A}\mathopen{}\mathclose{\left(\mathbf{x}_n, \mathbf{B}_n, \boldsymbol{\Gamma}_{n-1}, W_n\right)}$                                      |

The wake $W_n$, the bound geometry $\mathbf{B}_n$, and the strengths $\boldsymbol{\Gamma}_n$ are structured aggregates the functions pass and return, distinct from the iterated state vector $\mathbf{x}_n$ and not subject to its dot-product arithmetic (Section 5). The bound geometry is not independent: $\mathbf{B}_n = \operatorname{G}\mathopen{}\mathclose{\left(\mathbf{x}_n, t_n\right)}$.

### Initialization

At the first snapshot ($n = 0$): $\boldsymbol{\Gamma}_{-1} := 0$ and $W_0 := \emptyset$. The first official solve is therefore $\boldsymbol{\Gamma}_0 = \operatorname{C}\mathopen{}\mathclose{\left(\mathbf{x}_0, \mathbf{B}_0, \emptyset\right)}$ and $\boldsymbol{\ell}_0 = \operatorname{L}\mathopen{}\mathclose{\left(\mathbf{x}_0, \mathbf{B}_0, \boldsymbol{\Gamma}_0, 0, \emptyset\right)}$ with $\mathbf{B}_0 = \operatorname{G}\mathopen{}\mathclose{\left(\mathbf{x}_0, t_0\right)}$, whose unsteady term $\dot{\boldsymbol{\Gamma}}_0 \approx \boldsymbol{\Gamma}_0 / \Delta t$ is the impulsive-start convention.

## 3. Functions

### Bound geometry function

$$\mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} = \operatorname{G}\mathopen{}\mathclose{\left(\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, t_{n+1}\right)}$$

$\operatorname{G}$ returns the bound ring vortex positions in GP1 axes. The time tag sets the prescribed flapping geometry: the quarter-chord front legs and the interior back legs are functions of $t_{n+1}$ alone. The state enters only through the trailing-edge rings, whose back legs are set back from the prescribed trailing edge by a quarter of the trailing edge's per-step travel (a freestream-dependent offset, after Lambert, for drag accuracy). That freestream depends on the velocity and attitude blocks of $\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$. So $\operatorname{G}$ is iterate-dependent only in its trailing-edge rows, and for a rigid wing in a steady freestream it collapses to a pure function of time. The function is shown at its sub-iteration instantiation. The snapshot-$n$ geometry is $\mathbf{B}_n = \operatorname{G}\mathopen{}\mathclose{\left(\mathbf{x}_n, t_n\right)}$.

### Structural propagator function

$$\mathbf{x}_{n+1} = \operatorname{S}\mathopen{}\mathclose{\left(\mathbf{x}_n, \hat{\boldsymbol{\ell}}_{n+1}\right)}$$

$\operatorname{S}$ integrates the rigid-body dynamics in MuJoCo over $\mathopen{}\mathclose{\left(t_n, t_{n+1}\right]}$ with the applied wrench held constant at the interval load $\hat{\boldsymbol{\ell}}_{n+1} \in \mathbb{R}^6$ (Earth axes, gravity and external loads included). The hat distinguishes the interval load from the official snapshot load $\boldsymbol{\ell}_{n+1}$ of Section 6: distinct objects sharing the same time tag.

### Wake propagator function

$$W^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} = \operatorname{P}\mathopen{}\mathclose{\left(\mathbf{x}_n, \mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \mathbf{B}_n, \boldsymbol{\Gamma}_n, \mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, W_n\right)}$$

The wake is stored in GP1 axes (the first `Airplane`'s geometry axes, which move with the body), matching the solver. $\operatorname{P}$ has two parts. Convection: displace the pre-existing rings of $W_n$ over the interval by their transport velocity, which splits by what it depends on. An induced part (the bound and wake Biot-Savart induction, plus the image system when a surface is active) is evaluated from the snapshot-$n$ bound geometry and strengths and the wake, $\mathopen{}\mathclose{\left(\mathbf{B}_n, \boldsymbol{\Gamma}_n, W_n\right)}$, the standard explicit fluid transport. The snapshot-$n$ state $\mathbf{x}_n$ enters this part only through the image, setting the snapshot-$n$ surface reflection. In free air the induced velocities are computed in GP1 from relative positions alone, and $\mathbf{x}_n$ drops out. A frame part (the freestream re-expressed in GP1 and the apparent velocity $-\boldsymbol{\omega}\times\mathbf{r}$) is evaluated at the trial interval-end state $\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$, because the wake lives in body-attached axes: holding each ring fixed in the Earth frame means removing the body's own motion across the interval, and that motion is the one ending at $\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$, the state the step converges to and against which the loads are evaluated. Shed: append the interval-$\mathopen{}\mathclose{\left(n{+}1\right)}$ trailing-edge row, with strengths $\boldsymbol{\Gamma}_n$ from the snapshot-$n$ trailing-edge bound rings (Kelvin condition) and leading edges at the trailing-edge rows of $\mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$, which themselves carry the trial freestream offset (the geometry function $\operatorname{G}$).

$\operatorname{P}$ is therefore iterate-dependent through both its frame part and the trial bound geometry $\mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$ that sets its shed leading edges, so $W^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$ moves with the trial state and is re-evaluated every sub-iteration rather than frozen once. The split keeps that cheap: the induced transport depends only on the frozen snapshot-$n$ data $\mathopen{}\mathclose{\left(\mathbf{x}_n, \mathbf{B}_n, \boldsymbol{\Gamma}_n, W_n\right)}$ (with $\mathbf{x}_n$ used only by the surface image), so it is computed once at freeze (Section 5) and reused unchanged, while the trial-dependent pieces (the frame displacement, a freestream rotation and one cross product per wake point, and the trailing-edge shed leading edges from $\operatorname{G}$) are cheap to recompute. Among the shed row, only the strengths $\boldsymbol{\Gamma}_n$ are iterate-independent. Its leading-edge geometry moves with the trial freestream.

Frame handling: expressing the wake in GP1 at the trial interval-end state keeps the wake-body relative geometry consistent with the state at which the circulation solve and load evaluation run, so at convergence the wake is positioned relative to the accepted $\mathbf{x}_{n+1}$ rather than lagging a frame behind. This consistency is why the frame part takes the end state. The convection itself remains a single explicit displacement per step.

### Circulation solve

$$\boldsymbol{\Gamma}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} = \operatorname{C}\mathopen{}\mathclose{\left(\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, W^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}\right)}$$

Build the snapshot-$\mathopen{}\mathclose{\left(n{+}1\right)}$ operating point from $\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$. Assemble the matrix of aerodynamic influence coefficients (AIC) at step $n{+}1$ and trial $k$, $\mathbf{A}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$, from the bound geometry $\mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$, and set the right-hand side by the kinematic normal wash and the induction of the wake $W^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$. Solve for the bound strengths $\boldsymbol{\Gamma}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$. The function is shown at its sub-iteration instantiation. The official snapshot solve (Section 6) calls the same function at the accepted state.

### Load evaluation

$$\boldsymbol{\ell}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} = \operatorname{L}\mathopen{}\mathclose{\left(\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \boldsymbol{\Gamma}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \boldsymbol{\Gamma}_n, W^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}\right)}$$

The Kutta-Joukowski term uses the current snapshot's strengths $\boldsymbol{\Gamma}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$ directly (local velocities include the kinematic velocity from $\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$ plus bound and wake induction). The unsteady Bernoulli term uses $\dot{\boldsymbol{\Gamma}}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} \approx \mathopen{}\mathclose{\left(\boldsymbol{\Gamma}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} - \boldsymbol{\Gamma}_n\right)}/\Delta t$. Add any external loads from `external_loads_fn`, returned in wind axes and evaluated at the same operating point. Transform to Earth axes. Add gravity. Note the asymmetry: $\boldsymbol{\Gamma}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$ enters both terms, while $\boldsymbol{\Gamma}_n$ enters only through the finite difference.

### Composite aerodynamic evaluation

For the fixed-point algebra it is convenient to compose the circulation solve and the load evaluation into one function, keeping the bound geometry explicit:

$$\begin{aligned}
\boldsymbol{\ell}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} &= \operatorname{A}\mathopen{}\mathclose{\left(\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \boldsymbol{\Gamma}_n, W^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}\right)} \\
&:= \operatorname{L}\mathopen{}\mathclose{\left[\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \operatorname{C}\mathopen{}\mathclose{\left(\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, W^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}\right)}, \boldsymbol{\Gamma}_n, W^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}\right]}
\end{aligned}$$

$\boldsymbol{\Gamma}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$ is thus an internal product of $\operatorname{A}$: a deterministic function of $\mathopen{}\mathclose{\left(\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, W^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}\right)}$, not an independent input. The bound geometry and wake are internal products as well, $\mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} = \operatorname{G}\mathopen{}\mathclose{\left(\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, t_{n+1}\right)}$ and $W^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} = \operatorname{P}\mathopen{}\mathclose{\left(\mathbf{x}_n, \mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \mathbf{B}_n, \boldsymbol{\Gamma}_n, \mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, W_n\right)}$. They are kept as explicit arguments for the algorithm but are not free.

To confirm that $\operatorname{A}$ collapses to a function of the trial state alone, expand those internal products. Substitute $\mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} = \operatorname{G}\mathopen{}\mathclose{\left(\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, t_{n+1}\right)}$ in every slot it occupies (the explicit one, inside $\operatorname{C}$, and inside $\operatorname{P}$) and $W^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$ in both of its occurrences (one shared object), giving the fully expanded load:

$$\begin{aligned}
\boldsymbol{\ell}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} = \operatorname{L}\Bigl( &\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \operatorname{G}\mathopen{}\mathclose{\left(\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, t_{n+1}\right)}, \operatorname{C}\Bigl\{\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \operatorname{G}\mathopen{}\mathclose{\left(\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, t_{n+1}\right)}, \operatorname{P}\Bigl[\mathbf{x}_n, \mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \mathbf{B}_n, \boldsymbol{\Gamma}_n, \dots \\
&\operatorname{G}\mathopen{}\mathclose{\left(\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, t_{n+1}\right)}, W_n\Bigr]\Bigr\}, \boldsymbol{\Gamma}_n, \operatorname{P}\mathopen{}\mathclose{\left[\mathbf{x}_n, \mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \mathbf{B}_n, \boldsymbol{\Gamma}_n, \operatorname{G}\mathopen{}\mathclose{\left(\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, t_{n+1}\right)}, W_n\right]}\Bigr)
\end{aligned}$$

Every leaf is now either the trial state $\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$ or one of the frozen step quantities $\mathbf{x}_n, \mathbf{B}_n, \boldsymbol{\Gamma}_n, W_n, t_{n+1}$, which are constant across the sub-iteration (Section 4). Dropping those constants leaves the trial state as the only varying argument:

$$\boldsymbol{\ell}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} = \operatorname{A}_n\mathopen{}\mathclose{\left(\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}\right)}$$

This is a single-variable function over the sub-iteration, with the subscript $n$ carrying the frozen step data exactly as $\operatorname{H}_n$ does (Section 4). Two frozen leaves are worth calling out: $\mathbf{x}_n$ appears (inside $\operatorname{P}$, for the surface image only), but as a fixed parameter rather than a second free variable, so it does not compromise the single-argument reduction. $\boldsymbol{\Gamma}_n$ appears twice, as $\operatorname{L}$'s explicit fourth argument and again as $\operatorname{P}$'s shed strength.

The official snapshot solve (Section 6) is the same composition evaluated at the accepted state.

### Hat convention

Each trial evaluates $\operatorname{A}$ against its own bound geometry $\mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$ and wake $W^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$, and the official snapshot solve (Section 6) rebuilds both at the accepted state. The only hatted object is the interval load $\hat{\boldsymbol{\ell}}_{n+1}$ (the load the body is advanced under), which matches the official $\boldsymbol{\ell}_{n+1}$ of Section 6 to within the convergence tolerance.

## 4. The coupled system and its fixed point

The discrete coupled equations for the step from $t_n$ to $t_{n+1}$ are:

$$\begin{aligned}
\mathbf{B}_{n+1} &= \operatorname{G}\mathopen{}\mathclose{\left(\mathbf{x}_{n+1}, t_{n+1}\right)} \\
W_{n+1} &= \operatorname{P}\mathopen{}\mathclose{\left(\mathbf{x}_n, \mathbf{x}_{n+1}, \mathbf{B}_n, \boldsymbol{\Gamma}_n, \mathbf{B}_{n+1}, W_n\right)} \\
\mathbf{x}_{n+1} &= \operatorname{S}\mathopen{}\mathclose{\left(\mathbf{x}_n, \hat{\boldsymbol{\ell}}_{n+1}\right)} \\
\hat{\boldsymbol{\ell}}_{n+1} &= \operatorname{A}\mathopen{}\mathclose{\left(\mathbf{x}_{n+1}, \mathbf{B}_{n+1}, \boldsymbol{\Gamma}_n, W_{n+1}\right)}
\end{aligned}$$

The first two equations are no longer free of unknowns: the bound geometry and the wake both depend on $\mathbf{x}_{n+1}$. Only the snapshot-$n$ data $\mathopen{}\mathclose{\left(\mathbf{x}_n, \mathbf{B}_n, \boldsymbol{\Gamma}_n, W_n\right)}$ is frozen for the step. Substituting through yields one implicit equation in the state:

$$\begin{aligned}
\mathbf{x}_{n+1} &= \operatorname{H}_n\mathopen{}\mathclose{\left(\mathbf{x}_{n+1}\right)} \\
\operatorname{H}_n\mathopen{}\mathclose{\left(\mathbf{x}_{n+1}\right)} &:= \operatorname{S}\mathopen{}\mathclose{\left[\mathbf{x}_n, \operatorname{A}\mathopen{}\mathclose{\left(\mathbf{x}_{n+1}, \mathbf{B}_{n+1}, \boldsymbol{\Gamma}_n, W_{n+1}\right)}\right]}
\end{aligned}$$

Here $\mathbf{B}_{n+1} = \operatorname{G}\mathopen{}\mathclose{\left(\mathbf{x}_{n+1}, t_{n+1}\right)}$ and $W_{n+1} = \operatorname{P}\mathopen{}\mathclose{\left(\mathbf{x}_n, \mathbf{x}_{n+1}, \mathbf{B}_n, \boldsymbol{\Gamma}_n, \mathbf{B}_{n+1}, W_n\right)}$ are themselves functions of $\mathbf{x}_{n+1}$. $\operatorname{H}_n: \mathbb{R}^{12} \to \mathbb{R}^{12}$ is a self-map, and the snapshot-$\mathopen{}\mathclose{\left(n{+}1\right)}$ state is its fixed point. The subscript on $\operatorname{H}_n$ records that this step's frozen data $\mathopen{}\mathclose{\left(\mathbf{x}_n, \mathbf{B}_n, \boldsymbol{\Gamma}_n, W_n\right)}$ is baked into the self-map. The unknown is $\mathbf{x}_{n+1}$. The trial bound geometry, strengths, wake, and loads are produced as byproducts, with the official $\boldsymbol{\Gamma}_{n+1}$ and $\boldsymbol{\ell}_{n+1}$ finalized at the accepted state in Section 6.

The stored and relaxed variable is the increment:

$$\Delta \mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} := \mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} - \mathbf{x}_n$$

The absolute state is reconstructed only inside calls to $\operatorname{S}$ and $\operatorname{A}$. Because the sub-iteration algorithm touches only differences (residuals and inter-iterate updates), the increment and absolute parameterizations generate identical iterate sequences. The increment is preferred for floating-point hygiene and to keep Euler-angle arithmetic in a small, wrap-free neighborhood.

## 5. Sub-iteration algorithm (solve for snapshot $n{+}1$)

### Freeze

Save the MuJoCo state at $\mathbf{x}_n$. The step's frozen data is the snapshot-$n$ state $\mathbf{x}_n$ (used only by the surface image), bound geometry $\mathbf{B}_n$, strengths $\boldsymbol{\Gamma}_n$, and wake $W_n$. The expensive, iterate-independent part of $\operatorname{P}$ (the induced Biot-Savart transport, which reads only that frozen data) is precomputed once here, while the trial-dependent parts of $\operatorname{P}$ (the frame displacement and the $\mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$-set shed leading edges) recompute each iteration. The simulation clock advances only on acceptance.

### Seed (k = 0)

$$\mathbf{x}^{\mathopen{}\mathclose{\left(0\right)}}_{n+1} = \operatorname{S}\mathopen{}\mathclose{\left(\mathbf{x}_n, \boldsymbol{\ell}_n\right)}$$

The snapshot-$n$ load is recycled as the zeroth trial for the interval load: precisely the loose scheme's mis-tagging, used once as a warm start. The seed affects only the path of the iteration, not its fixed point.

### Loop

For $k = 0, 1, 2, \dots, k_{\text{max}}$:

1. At the trial state, build the trial bound geometry and wake, $\mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} = \operatorname{G}\mathopen{}\mathclose{\left(\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, t_{n+1}\right)}$ and $W^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} = \operatorname{P}\mathopen{}\mathclose{\left(\mathbf{x}_n, \mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \mathbf{B}_n, \boldsymbol{\Gamma}_n, \mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, W_n\right)}$, then evaluate the load:

   $$\boldsymbol{\ell}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} = \operatorname{A}\mathopen{}\mathclose{\left(\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \boldsymbol{\Gamma}_n, W^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}\right)}$$

   This yields the trial strengths $\boldsymbol{\Gamma}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} = \operatorname{C}\mathopen{}\mathclose{\left(\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, W^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}\right)}$ internally.
2. Dynamics under the trial load (MuJoCo restored to $\mathbf{x}_n$ first):

   $$\tilde{\mathbf{x}}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} = \operatorname{S}\mathopen{}\mathclose{\left(\mathbf{x}_n, \boldsymbol{\ell}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}\right)}$$
3. Residual (two candidates for the same snapshot):

   $$\mathbf{r}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} = \tilde{\mathbf{x}}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} - \mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$$
4. Test for convergence. Accept at $K_{n+1} := k$ and exit once the weighted residual satisfies the bound:

   $$\mathopen{}\mathclose{\left\lVert \mathbf{D} \mathbf{r}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} \right\rVert} \le \varepsilon_{\text{rel}} \mathopen{}\mathclose{\left\lVert \mathbf{D} \Delta \mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} \right\rVert} + \varepsilon_{\text{abs}}$$
5. Compute the relaxation factor. $\lambda^{\mathopen{}\mathclose{\left(0\right)}}_{n+1} = \lambda_{\text{init}}$. For $k \ge 1$, Aitken $\Delta^2$ in the $\mathbf{D}$-weighted dot product:

   $$\lambda^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} = -\lambda^{\mathopen{}\mathclose{\left(k-1\right)}}_{n+1} \frac{\mathopen{}\mathclose{\left( \mathbf{D} \mathbf{r}^{\mathopen{}\mathclose{\left(k-1\right)}}_{n+1} \right)} \cdot \mathopen{}\mathclose{\left[ \mathbf{D}\mathopen{}\mathclose{\left(\mathbf{r}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} - \mathbf{r}^{\mathopen{}\mathclose{\left(k-1\right)}}_{n+1}\right)} \right]}}{\mathopen{}\mathclose{\left\lVert \mathbf{D}\mathopen{}\mathclose{\left(\mathbf{r}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} - \mathbf{r}^{\mathopen{}\mathclose{\left(k-1\right)}}_{n+1}\right)} \right\rVert}^2}$$

   The factor falls back to $\lambda_{\text{init}}$ if the denominator is small relative to the residual, $\mathopen{}\mathclose{\left\lVert \mathbf{D}\mathopen{}\mathclose{\left(\mathbf{r}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} - \mathbf{r}^{\mathopen{}\mathclose{\left(k-1\right)}}_{n+1}\right)} \right\rVert}^2 \le \varepsilon_{\text{div}} \mathopen{}\mathclose{\left\lVert \mathbf{D} \mathbf{r}^{\mathopen{}\mathclose{\left(k-1\right)}}_{n+1} \right\rVert}^2$ (Section 8). All differences are taken across $k$ at fixed tag $n{+}1$.
6. Update:

   $$\mathbf{x}^{\mathopen{}\mathclose{\left(k+1\right)}}_{n+1} = \mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} + \lambda^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} \mathbf{r}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$$

Log $\mathopen{}\mathclose{\left\lVert \mathbf{D} \mathbf{r}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} \right\rVert}$ and $\lambda^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$ at every iteration for diagnostics.

## 6. Acceptance and handoff to the outer loop

At $k = K_{n+1}$:

- Interval load: $\hat{\boldsymbol{\ell}}_{n+1} := \boldsymbol{\ell}^{\mathopen{}\mathclose{\left(K_{n+1}\right)}}_{n+1}$.
- State: $\mathbf{x}_{n+1} := \tilde{\mathbf{x}}^{\mathopen{}\mathclose{\left(K_{n+1}\right)}}_{n+1}$, so that $\mathbf{x}_{n+1} = \operatorname{S}\mathopen{}\mathclose{\left(\mathbf{x}_n, \hat{\boldsymbol{\ell}}_{n+1}\right)}$ holds exactly (the accepted state is a genuine MuJoCo trajectory output), and $\hat{\boldsymbol{\ell}}_{n+1} = \operatorname{A}\mathopen{}\mathclose{\left(\mathbf{x}_{n+1}, \mathbf{B}_{n+1}, \boldsymbol{\Gamma}_n, W_{n+1}\right)}$ holds to tolerance.

The outer chain then proceeds with the official snapshot solve, rebuilding the bound geometry and wake at the accepted state:

$$\begin{aligned}
\mathbf{B}_{n+1} &= \operatorname{G}\mathopen{}\mathclose{\left(\mathbf{x}_{n+1}, t_{n+1}\right)} \\
W_{n+1} &= \operatorname{P}\mathopen{}\mathclose{\left(\mathbf{x}_n, \mathbf{x}_{n+1}, \mathbf{B}_n, \boldsymbol{\Gamma}_n, \mathbf{B}_{n+1}, W_n\right)} \\
\boldsymbol{\Gamma}_{n+1} &= \operatorname{C}\mathopen{}\mathclose{\left(\mathbf{x}_{n+1}, \mathbf{B}_{n+1}, W_{n+1}\right)} \\[6pt]
\boldsymbol{\ell}_{n+1} &= \operatorname{L}\mathopen{}\mathclose{\left(\mathbf{x}_{n+1}, \mathbf{B}_{n+1}, \boldsymbol{\Gamma}_{n+1}, \boldsymbol{\Gamma}_n, W_{n+1}\right)} \\
&= \operatorname{A}\mathopen{}\mathclose{\left(\mathbf{x}_{n+1}, \mathbf{B}_{n+1}, \boldsymbol{\Gamma}_n, W_{n+1}\right)}
\end{aligned}$$

This $\boldsymbol{\ell}_{n+1}$ is recorded in the load history and seeds the next sub-iteration, $\mathbf{x}^{\mathopen{}\mathclose{\left(0\right)}}_{n+2} = \operatorname{S}\mathopen{}\mathclose{\left(\mathbf{x}_{n+1}, \boldsymbol{\ell}_{n+1}\right)}$. The solve is necessary, not a redundant re-derivation: alongside $\boldsymbol{\ell}_{n+1}$ it produces the canonical bound geometry $\mathbf{B}_{n+1}$ and strengths $\boldsymbol{\Gamma}_{n+1}$ at the accepted state $\mathbf{x}_{n+1}$, which the next step's wake advance and finite difference consume as $\mathbf{B}_n$ and $\boldsymbol{\Gamma}_n$. It would reproduce the final trial evaluation only at the exact fixed point. Because the loop stops within tolerance, the accepted state $\mathbf{x}_{n+1} = \tilde{\mathbf{x}}^{\mathopen{}\mathclose{\left(K_{n+1}\right)}}_{n+1}$ differs (within tolerance) from the trial state $\mathbf{x}^{\mathopen{}\mathclose{\left(K_{n+1}\right)}}_{n+1}$ that evaluation ran at. The bound geometry and wake are rebuilt at the accepted state too, since both move with the state, so this is a genuine re-evaluation rather than a reuse of frozen objects. That the baseline run loop already performs this solve each step is a convenience (no new code seam), not the reason it is computed.

Nothing constrains the step increment $\Delta \mathbf{x}_{n+1} = \mathbf{x}_{n+1} - \mathbf{x}_n$ itself: the residual compares two snapshot-$\mathopen{}\mathclose{\left(n{+}1\right)}$ candidates, and the increment converges to whatever the physics dictates.

### Convergence mismatch

The interval load $\hat{\boldsymbol{\ell}}_{n+1}$ is evaluated at the final trial state and the official $\boldsymbol{\ell}_{n+1}$ at the accepted state. The two states differ within the convergence tolerance, and that gap propagates through the bound geometry $\mathbf{B}$ and the wake $W$ (both functions of the state) as well as directly, so $\mathopen{}\mathclose{\left\lVert \hat{\boldsymbol{\ell}}_{n+1} - \boldsymbol{\ell}_{n+1} \right\rVert}$ is bounded by the convergence tolerance times the local sensitivity of $\operatorname{A}$. This gap is a consequence of stopping the iteration within tolerance, not a modeling error. It grows only when a step fails to converge, for example when the $k_{\text{max}}$ cap is hit, which the per-step $k_{\text{max}}$ warning (below) already flags.

### Non-convergence handling

If the loop reaches $k = k_{\text{max}}$ without passing the convergence test, the step is accepted anyway, with $\hat{\boldsymbol{\ell}}_{n+1}$ and $\mathbf{x}_{n+1}$ taken from the final iterate as above. Accepting an unconverged step is safe rather than corrupting. The accepted state is still a genuine MuJoCo trajectory output (it satisfies $\mathbf{x}_{n+1} = \operatorname{S}\mathopen{}\mathclose{\left(\mathbf{x}_n, \hat{\boldsymbol{\ell}}_{n+1}\right)}$ exactly), so the rigid-body state is never unphysical. The only quantity left imperfect is the aero-structure consistency, by exactly the residual at the cap. The seed is the loose (explicit) step, and every relaxed iterate moves the estimate closer to the fixed point, so a capped step is at least as accurate as the loose coupling the solver tolerated before strong coupling existed, and usually much closer. Each step also re-solves the fixed point from scratch on the genuine accepted state, with a fresh seed and the next wake, so the solver carries no unconverged residual forward and a brief stiff transient self-heals once it passes.

The cap therefore continues the simulation rather than aborting it. The only hard stop in the step is the orientation validation that rejects out-of-range Euler angles (Section 10, item 3). Continuing is the right default, but it must not be silent: every step that reaches $k_{\text{max}}$ emits a warning. A single cap during a transient is benign, while sustained caps are dangerous (in an added-mass-dominated regime they are the loose-coupling instability returning, since a capped step is only a partially-relaxed loose step), and the two are indistinguishable from a single step but separable by their frequency over the run, which per-step warnings expose. $k_{\text{max}}$ is a user-tunable parameter (default 20). Raising it trades iterations for a tighter per-step solve when convergence is merely slow, as it can be for very light vehicles (large $m_a/m$, Section 8).

## 7. Construction of the weighting matrix D

The weighting matrix $\mathbf{D}$ enters the convergence test (Section 5, step 4) and the Aitken dot products (Section 5, step 5). It exists to make those operations meaningful across a 12-vector whose blocks carry different units (m, rad, m/s, rad/s) and, under load perturbations, different powers of $\Delta t$. It is diagonal and block-constant, ordered as the state $\mathopen{}\mathclose{\left(\mathbf{p}, \boldsymbol{\theta}, \mathbf{v}, \boldsymbol{\omega}\right)}$:

$$\mathbf{D} = \operatorname{diag}\mathopen{}\mathclose{\left( \frac{1}{L_{\text{ref}}}\mathbf{I}_3, \frac{1}{\theta_{\text{ref}}}\mathbf{I}_3, \frac{\Delta t}{L_{\text{ref}}}\mathbf{I}_3, \frac{\Delta t}{\theta_{\text{ref}}}\mathbf{I}_3 \right)}$$

Here $\mathbf{I}_3$ is the 3x3 identity matrix, so each block's scalar weight is repeated on the three diagonal entries of that block.

Two ideas determine it: a $\Delta t$ factor on the derivative blocks (the time-scaling fix) and the reference scales $L_{\text{ref}}$ and $\theta_{\text{ref}}$ (the units fix).

### Time scaling: the delta-t factor on the derivative blocks

The convergence test and the Aitken dot products combine the four blocks directly, so any systematic disparity in their magnitudes biases both the relaxation factor and the stopping criterion toward whichever block is largest. Under a load perturbation the per-step response splits by power of $\Delta t$: each rate block moves as (generalized force / generalized inertia) $\times\Delta t$ (the force-over-mass acceleration for the velocity block, the torque-over-inertia angular acceleration for the angular-velocity block), while the corresponding configuration block moves as that same acceleration times $\Delta t^2$ (a second integration). An unweighted norm is therefore velocity-dominated, and increasingly so as the step is refined.

Multiplying each derivative block by $\Delta t$ converts a rate into the increment it produces over one step, $\mathbf{v} \mapsto \mathbf{v}\Delta t$ and $\boldsymbol{\omega} \mapsto \boldsymbol{\omega}\Delta t$, lifting both derivative blocks to $O\mathopen{}\mathclose{\left(\Delta t^2\right)}$ to match the configuration blocks. This leading-order scaling is integrator-independent, so it holds for MuJoCo's RK4 the same as for a single Euler step (Section 10, item 1 notes only that the velocity update itself is not a simple function of the displacement increment). It removes the magnitude disparity, which is a discretization artifact, without touching the gain disparity between blocks, which is genuine physics. The relaxation factor then tracks the stiff velocity subspace because it is stiff, not because of a step-size accident.

### Reference length

Use the `Airplane`'s reference chord, $L_{\text{ref}} := c_{\text{ref}}$, which defaults to the main `Wing`'s mean aerodynamic chord. It is the aerodynamically meaningful length of the unsteady problem, it is already an `Airplane` attribute, and it scales with the vehicle, so the same $\mathbf{D}$ serves a full-size flapper and a centimeter-scale device without retuning. It is a characteristic scale, not a coordinate range. That distinction drives the angle choice below.

### Reference angle

Choose $\theta_{\text{ref}}$ so the attitude block balances the position block. A load perturbation applies a force $\mathbf{f}$ at a characteristic aerodynamic moment arm of order $c_{\text{ref}}$, hence a torque of order $c_{\text{ref}}\mathopen{}\mathclose{\left\lVert\mathbf{f}\right\rVert}$. The $\Delta t$-scaled, per-step residual blocks are then of order $\mathopen{}\mathclose{\left(\mathopen{}\mathclose{\left\lVert\mathbf{f}\right\rVert}/m\right)}\Delta t^2$ for position and $\mathopen{}\mathclose{\left(c_{\text{ref}}\mathopen{}\mathclose{\left\lVert\mathbf{f}\right\rVert}/\bar{I}\right)}\Delta t^2$ for attitude, where $\bar{I}$ is a scalar inertia. Weighting by $1/L_{\text{ref}}$ and $1/\theta_{\text{ref}}$ and equating the two gives:

$$\begin{aligned}
\theta_{\text{ref}} &:= \frac{mc_{\text{ref}}^2}{\bar{I}} \\
\bar{I} &:= \frac{1}{3}\operatorname{tr}\mathopen{}\mathclose{\left(\mathbf{J}\right)}
\end{aligned}$$

This is the squared ratio of the reference chord to the radius of gyration $\sqrt{\bar{I}/m}$. Both the first `Airplane`'s mass $m$ and its body inertia $\mathbf{J}$ (the tag `I_BP1_CgP1`, in body axes about the CG) are already supplied to the MuJoCo model, so $\theta_{\text{ref}}$ introduces no new input. For a typical vehicle the gyration radius is a few chords, so $\theta_{\text{ref}}$ lands well below 1 rad, weighting attitude more heavily than position. That is consistent with attitude being the aerodynamically sensitive direction, since the $\boldsymbol{\theta}$ block sets the angle of attack.

The scalar $\bar{I}$ is the mean principal moment of inertia, a rotation invariant, deliberately chosen over the per-axis moments of inertia. The trace of a matrix equals the sum of its eigenvalues, and the eigenvalues of the symmetric inertia tensor are its principal moments of inertia, so the identity holds in any axes regardless of the products of inertia. The inertia `I_BP1_CgP1` is constant in body axes, but the state's $\boldsymbol{\theta}$ and $\boldsymbol{\omega}$ blocks live in Earth axes, so a per-axis angular weight would rotate with the body and make $\mathbf{D}$ iterate-dependent, violating its constant status (Section 1, point 4). A scalar invariant keeps $\mathbf{D}$ frozen for the run, at the cost of collapsing the roll, pitch, and yaw anisotropy. That is an acceptable trade for a norm whose only job is to prevent gross block imbalance.

### Properties and assumptions

1. $\mathbf{D}$ is constant per run but not universal, assembled once from $\Delta t$, $c_{\text{ref}}$, $m$, and $\bar{I}$, all fixed for a run, so it satisfies the constant status of Section 1, point 4. It is deliberately discretization-aware (it depends on $\Delta t$), because $\mathbf{D}$ exists precisely to cancel a discretization-induced imbalance. The construction assumes a fixed $\Delta t$. Adaptive stepping would make $\mathbf{D}$ vary step to step and would require revisiting.
2. The $\Delta t$ factors admit an equivalent reference-speed reading: they are the same as nondimensionalizing the rate blocks by a step-dependent reference speed $V_{\text{ref}} = L_{\text{ref}}/\Delta t$ and reference rate $\Omega_{\text{ref}} = \theta_{\text{ref}}/\Delta t$. A fixed physical reference speed would fix units but leave the velocity residual $O\mathopen{}\mathclose{\left(\Delta t\right)}$ against the position residual $O\mathopen{}\mathclose{\left(\Delta t^2\right)}$, reintroducing velocity dominance under refinement. Tying the speed to $\Delta t$ is what removes it.
3. Because $\boldsymbol{\theta}$ holds intrinsic zyx Euler angles rather than a rotation vector, the inertia-to-angle correspondence holds only through the configuration-dependent Euler kinematic map, so $\theta_{\text{ref}}$ is a characteristic balance, not an exact preconditioner. This is consistent with the small-angle, away-from-singularity regime already assumed in Section 10, item 3.

## 8. Convergence tolerances

The convergence test (Section 5, step 4):

$$\mathopen{}\mathclose{\left\lVert \mathbf{D} \mathbf{r}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} \right\rVert} \le \varepsilon_{\text{rel}} \mathopen{}\mathclose{\left\lVert \mathbf{D} \Delta \mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} \right\rVert} + \varepsilon_{\text{abs}}$$

It is the standard mixed local-error form of the ODE-integration literature: a relative term scaled by the iterate plus an absolute floor. Its one nonstandard feature is that $\mathbf{D}$ has already nondimensionalized every block (Section 7), so $\mathbf{D}\mathbf{r}$ and $\mathbf{D}\Delta\mathbf{x}$ are pure numbers, and the two tolerances are therefore universal dimensionless constants rather than quantities to be rescaled per vehicle or per unit system. The same pair serves a full-size flapper and a centimeter-scale device unchanged. The recommended values are $\varepsilon_{\text{rel}} = 10^{-6}$ and $\varepsilon_{\text{abs}} = 10^{-10}$, derived below.

### Relative tolerance

$\varepsilon_{\text{rel}}$ is bounded below by floating-point noise and above by the truncation error of the discretization. Stability imposes no bound: the strongly coupled step is stable whenever the sub-iteration converges, which is a condition on the relaxation factor (Section 9). Both bounds are dimensionless groups, so the admissible interval can be checked for the whole class of vehicles the solver targets.

The lower bound is set by floating-point noise. A single evaluation of $\operatorname{H}_n$ is $\operatorname{S}$ composed with $\operatorname{A}$, and the accuracy-limiting step inside $\operatorname{A}$ is the dense bound-circulation solve in $\operatorname{C}$, whose relative error is of order $\operatorname{\kappa}\mathopen{}\mathclose{\left(\mathbf{A}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}\right)}\varepsilon_{\text{mach}}$, with $\varepsilon_{\text{mach}} \approx 2.2 \times 10^{-16}$ the interval machine epsilon and $\operatorname{\kappa}\mathopen{}\mathclose{\left(\mathbf{A}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}\right)}$ the AIC condition number. The bound-influence matrix is well conditioned: preliminary testing puts $\operatorname{\kappa}\mathopen{}\mathclose{\left(\mathbf{A}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}\right)}$ at order $10$, with the caveat that finer meshes may carry higher values. MuJoCo's RK4 in double precision contributes only a few multiples of $\varepsilon_{\text{mach}}$, negligible by comparison. The residual $\mathbf{r} = \tilde{\mathbf{x}}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} - \mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$ is a difference of two such evaluations, so near the fixed point its $\mathbf{D}$-weighted magnitude bottoms out at a small multiple of $\operatorname{\kappa}\mathopen{}\mathclose{\left(\mathbf{A}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}\right)}\varepsilon_{\text{mach}}$, on the order of $10^{-15}$. Write $h := V\Delta t / c_{\text{ref}}$ for the per-step travel measured in reference chords. A representative $\mathbf{D}$-weighted increment $\mathopen{}\mathclose{\left\lVert \mathbf{D}\Delta\mathbf{x} \right\rVert}$ is on the order of $h$, which is on the order of $10^{-1}$ since the unsteady solver's time step is conventionally near one panel chord of free-stream travel and a `Wing` typically carries on the order of ten chordwise panels. The stopping level $\varepsilon_{\text{rel}}\mathopen{}\mathclose{\left\lVert\mathbf{D}\Delta\mathbf{x}\right\rVert}$ must sit well above the floor, so

$$\varepsilon_{\text{rel}} \gg \frac{\operatorname{\kappa}\mathopen{}\mathclose{\left(\mathbf{A}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}\right)}\varepsilon_{\text{mach}}}{h}$$

which is on the order of $10^{-14}$ for the representative values. The chosen $\varepsilon_{\text{rel}} = 10^{-6}$ clears it by eight orders.

The upper bound is set by accuracy relative to the discretization. The sub-iteration removes only the splitting (partition) error, while the step still carries the truncation error of the time discretization and of the vortex lattice. Once the splitting error left at stopping is small compared to that truncation error, tightening $\varepsilon_{\text{rel}}$ further buys only iterations. Let $\tau_{\text{trunc}}$ denote the relative local truncation error of a step: the error the discretization makes over one step divided by the step's increment, both in the $\mathbf{D}$-weighted norm. The splitting error left at stopping is the residual divided by the contraction margin of the iteration, which is of order one under Aitken relaxation, so

$$\varepsilon_{\text{rel}} \le \frac{\tau_{\text{trunc}}}{\eta}$$

with $\eta$ a safety margin of at least $10^2$. The aerodynamic side of the step is first order in time (the circulation time derivative is a backward difference and the wake advects by a forward Euler step) while MuJoCo's RK4 is fourth order, so the aerodynamic error dominates $\tau_{\text{trunc}}$. For a first-order method $\tau_{\text{trunc}}$ is on the order of $\Delta t / T$, with $T$ the timescale over which the circulation history changes and a method constant of at most order one. Substituting $\Delta t = h c_{\text{ref}} / V$ puts $\tau_{\text{trunc}}$ on the order of $h c_{\text{ref}} / \mathopen{}\mathclose{\left(V T\right)}$. For periodic motion $T$ is the period over $2\pi$, so $c_{\text{ref}} / \mathopen{}\mathclose{\left(V T\right)}$ is twice the reduced frequency based on the reference chord. It runs from about $10^{-2}$ for the rigid-body flight-dynamics modes of a glider to about $1$ for a flapping microrobot, and $h$ runs from about $2 \times 10^{-2}$ on a fine mesh to about $2 \times 10^{-1}$ on a coarse one. The bound is tightest for the slowest dynamics on the finest mesh, where $\tau_{\text{trunc}}$ is on the order of $4 \times 10^{-4}$ and the bound at $\eta = 10^2$ is $4 \times 10^{-6}$. The chosen $10^{-6}$ sits below that bound by a factor of $4$, so its margin against the truncation error itself is $400$ there and thousands elsewhere in the class. The admissible interval therefore never closes, from on the order of $10^{-14}$ to no lower than $4 \times 10^{-6}$, and $10^{-6}$ lies inside it for every vehicle in the class.

The vehicle dependence of the iteration is its rate: the contraction factor of the relaxed fixed-point iteration grows with the added-mass ratio $m_a/m$, so a very light vehicle, where $m_a/m$ can approach unity, may press against $k_{\text{max}}$. The lever is to raise $k_{\text{max}}$, and the per-step $k_{\text{max}}$ warning reports when the cap is the reason a step stopped.

### Absolute tolerance

$\varepsilon_{\text{abs}}$ exists for near-trim steps. As the motion over a step approaches trim, $\Delta\mathbf{x} \to 0$ and the relative term vanishes, so without a floor the test would demand an unreachable zero residual and exhaust $k_{\text{max}}$ every step. $\varepsilon_{\text{abs}}$ is the $\mathbf{D}$-weighted residual accepted once the increment is negligible.

The physical threshold fixes its value: take the smallest per-step motion still worth resolving relatively, $\delta_{\text{min}}$, to be on the order of $10^{-4}$, one part in $10^4$ of a reference chord per step. Below it the step is effectively trimmed, and $\varepsilon_{\text{abs}}$ is on the order of $\varepsilon_{\text{rel}}\delta_{\text{min}} = 10^{-6} \cdot 10^{-4} = 10^{-10}$. The floating-point floor only has to sit below this, and the well-conditioned solve above places it several orders lower, so the physical threshold governs and $\varepsilon_{\text{abs}} = 10^{-10}$ is comfortably reachable.

The floor must not interfere on ordinary steps. There the relative term $\varepsilon_{\text{rel}}\mathopen{}\mathclose{\left\lVert\mathbf{D}\Delta\mathbf{x}\right\rVert}$ is on the order of $10^{-6} \cdot 10^{-1} = 10^{-7}$, roughly three orders above $\varepsilon_{\text{abs}}$, so the absolute floor is dormant during normal flight and governs only genuinely near-trim steps, as intended.

### Relaxation-factor denominator guard

The Aitken update (Section 5, step 5) divides by $\mathopen{}\mathclose{\left\lVert \mathbf{D}\mathopen{}\mathclose{\left(\mathbf{r}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} - \mathbf{r}^{\mathopen{}\mathclose{\left(k-1\right)}}_{n+1}\right)} \right\rVert}^2$, the squared norm of the change in residual between consecutive sub-iterations. This is a distinct quantity from the residual magnitude the convergence test bounds, and it can collapse while the residual is still well above tolerance. A mode of $\operatorname{H}_n$ with eigenvalue near unity leaves $\mathbf{r}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} \approx \mathbf{r}^{\mathopen{}\mathclose{\left(k-1\right)}}_{n+1}$, and an update $\lambda^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}\mathbf{r}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$ smaller than the state's last representable bit reproduces the previous iterate exactly and drives the denominator to zero. The guard reverts the factor to $\lambda_{\text{init}}$ in these cases, where the extrapolation is untrustworthy.

The criterion is relative to the residual:

$$\mathopen{}\mathclose{\left\lVert \mathbf{D}\mathopen{}\mathclose{\left(\mathbf{r}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} - \mathbf{r}^{\mathopen{}\mathclose{\left(k-1\right)}}_{n+1}\right)} \right\rVert}^2 \le \varepsilon_{\text{div}} \mathopen{}\mathclose{\left\lVert \mathbf{D} \mathbf{r}^{\mathopen{}\mathclose{\left(k-1\right)}}_{n+1} \right\rVert}^2$$

A relative test is scale-invariant as the residual shrinks along the iteration, so it reads the same at every sub-iteration: fall back when consecutive residuals agree to within a fraction $\sqrt{\varepsilon_{\text{div}}}$ of their size. The recommended value is $\varepsilon_{\text{div}} = 10^{-20}$, which triggers when they agree to about one part in $10^{10}$, comfortably above the relative residual-noise floor $\operatorname{\kappa}\mathopen{}\mathclose{\left(\mathbf{A}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}\right)}\varepsilon_{\text{mach}}$, on the order of $10^{-15}$, established above. The denominator is bimodal in practice, either of order $\mathopen{}\mathclose{\left\lVert \mathbf{D} \mathbf{r}^{\mathopen{}\mathclose{\left(k-1\right)}}_{n+1} \right\rVert}^2$ under normal progress or collapsed to the noise floor at a stall, with a wide gap between, so any value from $10^{-16}$ to $10^{-24}$ behaves identically and $\varepsilon_{\text{div}}$ is a safety net rather than a tuned parameter. It is always well defined when reached, because a failed convergence test guarantees $\mathopen{}\mathclose{\left\lVert \mathbf{D} \mathbf{r}^{\mathopen{}\mathclose{\left(k-1\right)}}_{n+1} \right\rVert} > \varepsilon_{\text{abs}}$, comfortably above noise.

## 9. Initial relaxation factor

$\lambda_{\text{init}}$ is the relaxation factor applied to the first sub-iteration update only. From $k = 1$ onward the Aitken $\Delta^2$ formula (Section 5, step 5) recomputes the factor dynamically, so $\lambda_{\text{init}}$ sets just the $k = 0$ step, $\mathbf{x}^{\mathopen{}\mathclose{\left(1\right)}}_{n+1} = \mathbf{x}^{\mathopen{}\mathclose{\left(0\right)}}_{n+1} + \lambda_{\text{init}}\mathbf{r}^{\mathopen{}\mathclose{\left(0\right)}}_{n+1}$. Like the seed, it affects the path of the iteration, not the fixed point it converges to, and because Aitken adapts immediately afterward, it need not be the optimal factor, only a safe warm start. The choice is correspondingly forgiving.

Under-relaxation is the stabilizing direction. The un-relaxed step $\lambda = 1$ is plain fixed-point iteration on $\operatorname{H}_n$, which diverges once the added-mass feedback gain (of order $m_a/m$) exceeds unity (Section 10, item 1). Taking $\lambda < 1$ damps the update and shrinks the effective amplification. A value of $0.5$ sits firmly in the damped regime. The precise stable upper bound tightens below $1$ as $m_a/m$ grows, so $\lambda < 1$ is not an unconditional guarantee for an arbitrarily light vehicle, but since $\lambda_{\text{init}}$ seeds only the first step and Aitken drives the factor from there, the warm start need only be conservative, not globally stable.

The chosen $\lambda_{\text{init}} = 0.5$ is deliberately more conservative than comparable tools. SHARPy fixes its relaxation factor at $0.8$. We start lower because the target regime here is lighter vehicles with higher $m_a/m$ (Section 8), whereas SHARPy's very-flexible, large-aircraft regime carries a lower added-mass ratio and tolerates a more aggressive factor. The comparison is directional rather than exact: SHARPy relaxes a load-based interface residual, whereas the factor here relaxes the body-state 12-vector residual, so equal numerical values do not represent equal damping of the same quantity.

## 10. Practical notes

1. Velocities stay in the iterate because the added-mass feedback loop runs through the velocity block ($\delta \boldsymbol{\ell} \to \delta \mathbf{v}_{n+1} \to \delta \boldsymbol{\Gamma}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} \to \delta \boldsymbol{\ell}$, with step-size-independent gain on the order of $m_a/m$), so relaxing displacements while passing velocities through unrelaxed reintroduces the loose-coupling divergence. Slaving velocities to displacement increments via finite differences is exact only for a single semi-implicit Euler step and breaks for MuJoCo's default multi-stage RK4 integrator, whose velocity update is not a simple function of the displacement increment. Rejected. The state-dependent bound geometry and wake open a second feedback path, $\delta \mathbf{x} \to \delta \mathbf{B}, \delta W \to \delta \boldsymbol{\ell}$, but it is $O\mathopen{}\mathclose{\left(\Delta t\right)}$ geometric and subdominant to the velocity channel, whose gain is on the order of $m_a/m$, so it does not change this picture.
2. Velocity-only relaxation was considered and rejected. The converse splitting is sound: every closed feedback cycle through the displacement block of $\operatorname{H}_n$ carries the $O\mathopen{}\mathclose{\left(\Delta t^2\right)}$ configuration response to a load perturbation and is $O\mathopen{}\mathclose{\left(\Delta t\right)}$ overall, so relaxing only $\mathopen{}\mathclose{\left(\mathbf{v}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, \boldsymbol{\omega}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}\right)}$ while passing the configuration block through unrelaxed would converge, and would avoid relaxing the singular $\boldsymbol{\theta}$ block and the $\Delta t$-sensitive block structure of $\mathbf{D}$. Rejected in favor of relaxing the full 12-vector: full-state relaxation is the standard presentation in the partitioned-FSI literature, whereas the partial splitting would demand its own stability argument in peer review for a payoff that is only bookkeeping. The underlying observation (the stiff subspace is the velocity block) still informs the construction of $\mathbf{D}$ (Section 7).
3. The state's $\boldsymbol{\theta}$ block is intrinsic zyx Euler angles, converted to degrees only at the `OperatingPoint` seam, where each angle is validated to lie in $\mathopen{}\mathclose{\left(-180, 180\right]}$. We do not implement wrap-aware arithmetic on the $\boldsymbol{\theta}$ block. Residuals and updates use plain subtraction and addition, which is accurate as long as the orientation stays away from the representation singularities. If the relaxation drives any angle out of $\mathopen{}\mathclose{\left(-180, 180\right]}$, building the next `OperatingPoint` hard-errors, so a divergent orientation aborts the solve rather than corrupting it silently. Two in-range singularities are not caught by that validation and are accepted as out of scope for now: the gimbal-lock pole at pitch $= \pm 90^\circ$ (it passes validation, and the angle extraction resolves it with a non-erroring fallback, so the solve degrades silently rather than aborting), and the $\pm 180^\circ$ wrap seam (two candidates straddling it both pass validation but differ by nearly $360^\circ$ under plain arithmetic). These are acceptable because the intended flight regimes stay well away from both. Extreme or aerobatic orientations are not supported.
4. AIC reuse across sub-iterations is only partial, because the bound geometry is not fully iterate-independent. The AIC is assembled from $\mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1} = \operatorname{G}\mathopen{}\mathclose{\left(\mathbf{x}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}, t_{n+1}\right)}$, whose trailing-edge back legs carry the freestream-dependent offset (the geometry function $\operatorname{G}$), so they move with the trial state. The rest of $\mathbf{B}^{\mathopen{}\mathclose{\left(k\right)}}_{n+1}$ (front legs, interior back legs) is fixed at $t_{n+1}$. In free air this trailing-edge dependence is the only trial dependence, so most of the matrix could be built once and only its trailing-edge contributions and the right-hand sides refreshed per trial. With a surface active, the image system's geometry in GP1 additionally depends on the trial attitude and position, spreading the dependence across the whole matrix. We do not exploit any of this for now: each sub-iteration rebuilds the AIC matrix in full. The partial reuse (freezing all but the trailing-edge rings in free air) is left as a future optimization, to be considered only after the initial implementation and testing.

## Inspirations

- Causin, P., Gerbeau, J. F. & Nobile, F. Added-mass effect in the design of partitioned algorithms for fluid-structure problems. *Computer Methods in Applied Mechanics and Engineering* **194**, 4506-4527 (2005).
- Kuettler, U. & Wall, W. A. Fixed-point fluid-structure interaction solvers with dynamic relaxation. *Comput Mech* **43**, 61-72 (2008).
- Carre, A. del, Munoz-Simon, A., Goizueta, N. & Palacios, R. SHARPy: A dynamic aeroelastic simulation toolbox for very flexible aircraft and wind turbines. *Journal of Open Source Software* **4**, 1885 (2019).
