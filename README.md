# 2D XFEM Fracture Lab

**[Open the interactive simulator →](https://aiden-azarnoush.github.io/xfem-fracture-lab/)**

Explore crack initiation and propagation in a rectangular elastic specimen.
The browser runs the Python extended finite element method (XFEM) solver
locally, with a fixed triangular mesh, selectable materials, and repeated
hammer strikes. This version models **Mode III, out-of-plane shear**, with
energy-selected, piecewise-linear crack paths that can turn through the material.
A second page, the **impact experiment**, uses bond-based peridynamics to show
dynamic cracking with branching and fragmentation from the same materials and
hammer.

## Explore the web app

1. Choose **Glass**, **Ceramic**, or **Polymer**, or select **Custom material**
   to enter Young’s modulus, Poisson’s ratio, and fracture energy. The presets
   are illustrative elastic properties, not certified material data.
2. Select **Coarse**, **Medium**, or **Fine** mesh. The specimen stays
   **300 × 200 mm**, with a thickness of **5 mm**. Toggle **Show mesh** to
   display or hide the element edges.
3. Choose the hammer mass and drop height. Move the pointer up or down over
   the specimen to aim, then **left-click** to strike. The **Impact position**
   slider and **Strike specimen** button provide the same controls.
4. Strike to start a crack or advance existing tips. Each strike can produce
   up to **three growth steps**, with a new direction selected at each step.
   Aim farther away to make a new crack eligible for initiation. The default
   glass specimen, **0.5 kg hammer and 12 cm drop**, initiates a crack on
   all three meshes. Existing cracks can also respond to a new impact.
5. Watch the crack count, longest crack, and **energy release / resistance**
   ratio. A ratio below 1 means the trial extension is subcritical; increase
   the drop height or hammer mass to encourage growth. The color map shows
   the last strike’s element-averaged shear-stress magnitude in MPa.
6. Continue until a crack reaches the numerical fracture limit, or use
   **Reset specimen** to start again. Open **About the model** for the
   assumptions and the meaning of the simplified hammer load.

Choose **Varied toughness** for fixed stronger and weaker regions, or
**Uniform toughness** to remove this variation. Switch to **Toughness map**
to see the regions: light is weaker, dark is stronger. The variation enters
the fracture calculation; the paths are not random drawing effects.

Changing the material, material structure, custom properties, or mesh starts a fresh specimen.
Changing the hammer mass, drop height, or impact position preserves cracks.
The page fits the viewport; on small screens **Settings** opens a panel
that can scroll independently. With the canvas focused, **Space/Enter**
strikes and **↑/↓** adjusts the impact position.

The Python solver runs in a background worker through Pyodide, keeping the
interface responsive during calculations. The first visit downloads the
Python runtime and NumPy; calculations then run in the browser without a
Python server or account.


## The impact experiment (peridynamics)

**[Open the impact page →](https://aiden-azarnoush.github.io/xfem-fracture-lab/impact.html)**

The XFEM page grows one crack at a time from an energy criterion, which is the
right lesson for stable crack growth. A hammer blow on glass is a different
event: stress waves race outward, running cracks become unstable and branch,
and small variations in the material decide where. The second page models
that with **bond-based peridynamics**, running in plain JavaScript:

- The plate is a cloud of nodes, each bonded to every neighbour within a
  horizon of three spacings. A bond carries a force proportional to its
  stretch and breaks permanently once the stretch exceeds
  s₀ = √(4πG<sub>c</sub>/(9Eδ)). Newton’s second law is integrated
  explicitly. Broken bonds *are* the cracks; nothing about their path is
  prescribed, and branching, merging, and fragmentation emerge on their own.
- Click anywhere on the plate to strike. Nodes under the hammer head are
  pushed outward as the dent forms, with a push that grows with the impact
  energy m·g·h. Strike again to add damage; **New plate** starts over.
- Same presets (Glass, Ceramic, Polymer, Custom), same three mesh densities
  (5, 3.5, 2.5 mm node spacing), same hammer controls plus a head radius.
- **Material structure** scatters the node strengths randomly, as in real
  glass; that scatter is what makes the cracks wander and fork.

> [!NOTE]
> The model is two-dimensional plane stress with Poisson’s ratio fixed at
> 1/3 by the bond-based formulation; the in-plane push stands in for the
> bending and contact stresses of a face-on blow and reproduces radial
> cracks and branching, not the concentric ring cracks of a true face-on
> impact. A tapered no-fail zone near the free edges suppresses spurious
> edge damage from reflected waves, and node positions are jittered so
> cracks do not follow the grid. It is a teaching model, not a validated
> shatter simulation.

> [!TIP]
> Hit the glass preset twice in different places and watch the second set of
> cracks find the first. Then switch to the polymer preset: the same blow
> does nothing, because its fracture energy is thirty times higher.

## Physical model and normalization

The specimen occupies 0 ≤ x ≤ W and 0 ≤ y ≤ H. The left edge is fixed.
The scalar unknown w(x,y) is displacement perpendicular to the displayed
plane; it is not an in-plane displacement or a plate-bending solution.
For Young’s modulus E and Poisson’s ratio ν:

```math
\mu=\frac{E}{2(1+\nu)},\qquad
\boldsymbol{\tau}=\mu\nabla w,\qquad
\nabla\cdot(\mu\nabla w)=0.
```

The enriched approximation uses linear triangle shape functions Nᵢ,
shifted Heaviside functions across crack faces, and the Mode III
square-root function near each tip:

```math
w_h(\mathbf{x})=
\sum_i N_i(\mathbf{x})w_i
+\sum_c\sum_{j\in J_c}N_j(\mathbf{x})
\left[H_c(\mathbf{x})-H_c(\mathbf{x}_j)\right]a_{jc}
+\sum_c\sum_{k\in T_c}N_k(\mathbf{x})
\left[F_c(\mathbf{x})-F_c(\mathbf{x}_k)\right]b_{kc},
```

```math
F_c(r,\theta)=\sqrt{\frac{r}{h_e}}\sin\frac{\theta}{2},\qquad
K_{ij}=t\int_\Omega\mu\nabla\Phi_i\cdot\nabla\Phi_j\,dA.
```

Here t is thickness, hₑ is the horizontal element width, and Φᵢ denotes
an ordinary or enriched basis function. The tip’s local axis follows its
most recent segment. Tip enrichment covers nodes within 2hₑ of the tip.
The tip function has the displayed asymptotic form near the tip; away from
it, its sign is adjusted to route the displacement jump along the actual
polyline rather than an infinite extension of the last segment.
Heaviside enrichment also follows the polyline. These additional unknowns enter the stiffness matrix and
the elastic solution used to decide growth.

### Mesh choices

| Mesh | Grid divisions | Triangles | Ordinary nodes | Growth increment |
|---|---:|---:|---:|---:|
| Coarse | 12 × 8 | 192 | 117 | 16.25 mm |
| Medium | 18 × 12 | 432 | 247 | 10.83 mm |
| Fine | 24 × 16 | 768 | 425 | 8.13 mm |

The background mesh stays fixed as cracks grow. Cut elements receive
integration subcells, without remeshing the specimen. Mesh selection also
changes the finite trial-notch length and growth increment, so the results
are mesh dependent; refinement alone is not a validation claim.

## Hammer energy and crack growth

The drop height h_d and hammer mass m give an input energy. An assumed
stopping distance δ = **0.5 mm** converts it into a peak-equivalent load:

```math
E_{\mathrm{drop}}=mgh_d,\qquad
F=\frac{E_{\mathrm{drop}}}{\delta},\qquad g=9.81\ \mathrm{m/s^2}.
```

Two smooth, opposing sets of nodal loads on the right edge represent a
tearing impulse around the selected impact height. Each set has resultant
magnitude F. This is an illustrative load calibration; the hammer’s motion,
contact force history, and rebound are not resolved.

At each growth step the solver compares trial extensions at every active
crack tip. Candidate turns are **0°, ±25°, and ±50°** relative to the current
heading. Directions more than 65° away from the global leftward axis are
excluded so each path remains a single-valued graph without looping.
Each admissible trial requires a new elastic solve under the same applied
load. The load remains at the selected impact height; it is not snapped
to an existing crack.

```math
U=\frac12\mathbf{f}^{\mathsf T}\mathbf{u},\qquad
G\approx\max\left(0,\frac{U_{\mathrm{trial}}-U_{\mathrm{current}}}
{t\,\Delta a}\right).
```

The candidate with the largest G/Gc is selected. It is accepted only when
**G/Gc ≥ 1** and the remaining input-energy budget covers the new crack
area. A strike permits at most three accepted extensions, each costing
Gc × t × Δa from that budget. An existing-tip increment is **0.65hₑ**;
a new trial notch is **1.3hₑ** long. New nucleation is considered only on
the first step, when the impact is at least **1.5hₑ** from every existing
crack mouth. Existing tips compete with the new trial notch, so a remote
strike does not guarantee a new crack.

The selected base fracture energy is multiplied by a fixed smooth field:

```math
G_c(x,y)=G_{c0}\exp\left[v\left(0.6\sin(80x+23y)
+0.4\sin(44x-95y+1.7)\right)\right].
```

Coordinates are in meters. **Uniform toughness** uses v = 0;
**Varied toughness** uses v = 0.45. The latter is an illustrative spatial
variation, not a measured microstructure. Resistance along a candidate
segment is averaged at its quarter, midpoint, and three-quarter positions.

The displayed ratio is the last selected directional trial's G/Gc,
including a rejected trial if growth arrests after earlier accepted steps.
The longest-crack readout reports arc length along the polyline. The stress
map uses the final accepted geometry under the current strike load, with
its color scale normalized separately for each strike. Before the first
strike, the specimen displays the toughness field.

### Scope of the fracture model

Cracks start at the right edge, can turn upward or downward, and progress
leftward. Up to **eight separate cracks** interact through the global
elastic solve. A single tip does not split into branches. There is no
interior-impact nucleation, dynamic shattering, or intersection/coalescence
solver.

Growth stops at a numerical cutoff of **0.15hₑ** from an outer boundary
or before a proposed intersection with another crack. The interface calls
this **Fracture limit reached**, not a computed final separation. No crack
is artificially completed across the remaining ligament. Crack width is
exaggerated, and the growth animation is not physical time integration.

This model assumes linear elastic, quasi-static antiplane deformation.
It does not include plasticity, viscoelasticity, fatigue accumulation,
crack-face contact, or dynamic impact. Repeating a subcritical strike does
not accumulate damage. The polymer preset therefore illustrates an
elastic response rather than real polymer fracture behavior.

## Numerical method and verification

The Python module assembles the ordinary and enriched degrees of freedom
into a dense stiffness matrix. Integration polygons are split along crack
planes and through tip coordinates, triangulated, and integrated with
three-point triangle quadrature. Triangles near tips receive an additional
subdivision.

The system is diagonally scaled and stabilized with **10⁻¹¹** on the
scaled diagonal. A solution is accepted only when it is finite and its
relative free-degree-of-freedom equilibrium residual is below **10⁻⁴**.

The automated checks cover:

- conservation of area during subdivision around kinked cracks;
- rotated tip-function gradients against finite differences;
- displacement jumps along the actual polyline, without a false jump along
  the extrapolated tip tangent;
- linear elastic scaling with shear modulus and equilibrium residuals;
- approximate reflection symmetry of the intact specimen;
- fixed and uniform toughness fields;
- subcritical strikes, turning paths, and multiple cracks on **all three meshes**;
- sensitivity to impact position, boundary/intersection cutoffs, and arc length;
- rejection of invalid material and loading inputs.

```bash
python3 -m pip install -r python/requirements.txt
python3 -m unittest discover -s python -v
```

The earlier straight-path interface was checked in a browser. Browser
verification of this curved-path revision is incomplete; its Python tests
and JavaScript syntax are checked separately. These are consistency and
behavior checks, **not experimental validation or an independent crack-tip
benchmark**. Directional selection may be sensitive to discretization and
enrichment changes between trial geometries. Unresolved trial systems are
excluded and reported in the status message.

## Update the repository

Upload the contents of the complete ZIP, keeping the `python/` folder.
This revision requires **index.html, style.css, app.js, worker.js, and
python/xfem.py together**. Replacing only the HTML and README leaves the
old crack solver in place. The `?v=2` asset references refresh cached files.

The link at the top points to GitHub Pages. It shows the new version only
after you upload the files and GitHub finishes rebuilding the site.

## Run locally

```bash
python3 -m http.server 8000
```

Open **http://localhost:8000**. Serve the repository over HTTP rather than
opening `index.html` directly, because workers and Python-source fetching
require an HTTP origin. The static website also works on GitHub Pages.

The web worker downloads the pinned **Pyodide 0.27.7** runtime and NumPy
from jsDelivr. The local Python requirements are only needed when running
the solver or its checks outside the browser.

## Repository layout

```
index.html, style.css, app.js   XFEM page: interface, canvas, and hammer controls
impact.html, impact.js         impact page: bond-based peridynamics, all in JavaScript
worker.js                      Python runtime and background computation (XFEM page)
python/xfem.py                 XFEM solver and crack-growth calculation
python/test_xfem.py             numerical consistency and regression checks
python/requirements.txt        dependency for local Python use
LICENSE                        MIT license
```

All Python source stays in `python/`. The interface and static assets stay
at the repository root, matching the layout of the author's other web
calculators. GitHub Pages serves the root of the `main` branch.

## References

- [GetFEM: level sets and XFEM](https://getfem.org/userdoc/xfem.html) —
  enrichment and cut-cell integration concepts.
- [University of Liège: energy release rate in Mode III](https://www.ltas-cm3.ulg.ac.be/FractureMechanics/print.php?p=Lecture3_P4) —
  antiplane crack-tip fields and fracture energy.
- [Pyodide: using Python in a web worker](https://pyodide.org/en/stable/usage/webworker.html) —
  running the numerical calculation outside the interface thread.

## Author and license

Aiden Azarnoush. MIT — see [LICENSE](LICENSE).
