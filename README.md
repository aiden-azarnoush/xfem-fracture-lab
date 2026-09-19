# 2D Fracture Lab

**[Open the interactive simulator →](https://aiden-azarnoush.github.io/xfem-fracture-lab/)**

Strike a clamped plate with a hammer and watch cracks run, turn, and branch.
The specimen is a 300 × 200 mm plate, fixed along its left edge and struck on
its right edge at the height you aim at. Everything runs in the browser in
plain JavaScript: choose the material, the mesh, the hammer mass and drop
height, and click.

The cracks are computed with **bond-based peridynamics**, a dynamic model in
which the plate is a cloud of nodes bonded to their neighbours and a bond
breaks permanently once its stretch exceeds the critical value set by the
fracture energy. Broken bonds *are* the cracks; nothing about their path,
turning, or branching is prescribed.

## Explore the web app

1. Choose **Glass**, **Ceramic**, or **Polymer**, or **Custom material** to
   enter Young’s modulus, density, and fracture energy. The presets are
   illustrative, not certified material data. The bond-based model fixes
   Poisson’s ratio at 1/3.
2. Choose **Varied strength** (node strengths scattered randomly, as in real
   glass; this is what makes cracks wander and fork) or **Nearly uniform
   strength** to see how tidy the same blow becomes without it.
3. Select **Coarse**, **Medium**, or **Fine** mesh (5, 3.5, or 2.5 mm node
   spacing). Toggle **Show mesh** to see the nodes.
4. Choose the hammer mass and drop height. Move the pointer up or down over
   the specimen to aim, then **click** to strike (or use the Impact position
   slider and the Strike specimen button). The dent the hammer makes grows
   with the impact energy m·g·h.
5. Watch the cracks run over a few seconds. Strike again to grow them; the
   readouts give the number of cracks, the longest crack, and the fraction of
   broken bonds. When a crack reaches the fixed edge the specimen is broken.
6. Switch the color map between **Strain** (bond stretch during and after
   the strike) and **Strength** (the scattered strength field).

> [!TIP]
> Try the glass preset with the standard 0.5 kg hammer from 12 cm: one crack
> runs leftward from the hammer with small side branches. Raise it to 1 kg
> from 40 cm and the crack forks into a Y near the fixed edge. Then switch to
> the polymer preset: the same blow does nothing, because its fracture energy
> is thirty times higher.

## The model

- Nodes on a jittered grid; every pair closer than the horizon δ = 3 spacings
  is a bond. Jittering keeps cracks from following the grid.
- Bond force ∝ stretch with the plane-stress bond constant c = 9E/(πtδ³);
  a bond breaks for good when its stretch exceeds
  s₀ = √(4πG<sub>c</sub>/(9Eδ)). Node strengths are lognormally scattered
  (12 % for the varied option).
- Explicit velocity-Verlet time stepping at 0.4 of the wave-transit time per
  spacing, with light damping.
- The hammer head dents the right edge inward over a short ramp and holds;
  the dent depth scales with √(m·g·h). The left edge is clamped. The free
  top and bottom edges carry a tapered no-fail zone so waves reflecting from
  them do not produce spurious edge damage.

> [!NOTE]
> This is a teaching model, not a validated fracture simulation: two
> dimensional, plane stress, Poisson’s ratio fixed by the formulation, broken
> pieces stay in place. It was checked qualitatively against the known
> behaviour of the model (crack speed, branching under stronger loading,
> the polymer preset staying intact) rather than against experiments.

## The earlier XFEM version

The lab started as a quasi-static **Mode III XFEM** solver written in Python
and run in the browser through Pyodide (`worker.js`, `python/xfem.py`, tests
in `python/test_xfem.py`). That solver is kept in the repository for
reference and its tests still run locally, but the page no longer uses it:
a quasi-static, single-tip model cannot produce the dynamic branching that a
hammer blow on glass shows. The sections below document that version.

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
index.html, style.css, app.js   the web page: interface, canvas, and the peridynamics solver
worker.js, python/             the earlier XFEM version (not used by the page)
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
