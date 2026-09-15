# 2D XFEM Fracture Lab

**[Open the interactive simulator →](https://aiden-azarnoush.github.io/xfem-fracture-lab/)**

Explore crack initiation and propagation in a rectangular elastic specimen.
The browser runs the Python extended finite element method (XFEM) solver
locally, with a fixed triangular mesh, selectable materials, and repeated
hammer strikes. This version models **Mode III, out-of-plane shear**, with
prescribed straight crack paths.

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
4. Strike near an existing crack to extend it. Aim farther away to create
   another crack. The default glass specimen, **0.5 kg hammer and 12 cm
   drop**, initiates a crack on all three meshes.
5. Watch the crack count, longest crack, and **energy release / resistance**
   ratio. A ratio below 1 means the trial extension is subcritical; increase
   the drop height or hammer mass to encourage growth. The color map shows
   the last strike’s element-averaged shear-stress magnitude in MPa.
6. Continue until the specimen is classified as separated, or use
   **Reset specimen** to start again. Open **About the model** for the
   assumptions and the meaning of the simplified hammer load.

Changing the material, custom properties, or mesh starts a fresh specimen.
Changing the hammer mass, drop height, or impact position preserves cracks.
The page fits the viewport; on small screens **Settings** opens a panel
that can scroll independently. With the canvas focused, **Space/Enter**
strikes and **↑/↓** adjusts the impact position.

The Python solver runs in a background worker through Pyodide, keeping the
interface responsive during calculations. The first visit downloads the
Python runtime and NumPy; calculations then run in the browser without a
Python server or account.

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
an ordinary or enriched basis function. The tip’s local axis points left,
along the prescribed growth direction. Tip enrichment covers nodes within
2hₑ of the tip. Heaviside enrichment represents the displacement jump
behind the tip. These additional unknowns enter the stiffness matrix and
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

Each strike solves the current geometry and a virtual crack extension
under the same force. For this load-controlled linear elastic calculation:

```math
U=\frac12\mathbf{f}^{\mathsf T}\mathbf{u},\qquad
G\approx\max\left(0,\frac{U_{\mathrm{trial}}-U_{\mathrm{current}}}
{t\,\Delta a}\right).
```

An extension is accepted when **G ≥ Gc**, where Gc is the selected fracture
energy, and the input energy covers the new crack area:

```math
E_{\mathrm{drop}}\ge G_c\,t\,\Delta a.
```

The displayed ratio is G/Gc for the attempted extension, including rejected
attempts. Each strike advances an existing crack by at most **0.65hₑ**;
a new trial notch is **1.3hₑ** long. A strike within **1.5hₑ** of an existing
crack targets it. Strikes at other heights can initiate additional cracks.
The stress map uses the current strike load on the accepted geometry and
normalizes its color scale separately for each strike.

### Scope of the fracture model

Cracks start at the right edge and remain horizontal. Multiple cracks
interact through the global elastic solution, with at most 12 separate
cracks. There is no curved-path selection, branching, crack intersection,
or special coalescence treatment.

A remaining ligament of **one element width or less** is classified as
separated, and the displayed terminal crack is completed to the left edge.
The final stress map remains that of the last solved, nearly separated
specimen. This last visual completion is not a resolved ligament-failure
or fragment-motion calculation. Crack width is exaggerated for visibility.

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

The six automated checks cover:

- conservation of area during integration subdivision;
- the analytic tip-function gradient against finite differences;
- linear elastic scaling with shear modulus and equilibrium residuals;
- approximate reflection symmetry of the intact specimen;
- subcritical strikes, initiation, repeated growth, multiple cracks, and
  terminal separation on **all three meshes**;
- rejection of invalid material and loading inputs.

```bash
python3 -m pip install -r python/requirements.txt
python3 -m unittest discover -s python -v
```

Browser checks also cover custom-material editing, invalid inputs,
repeated strikes, multiple impact positions, and desktop/mobile layouts.
These are consistency and behavior checks, **not experimental validation
or an independent crack-tip benchmark**.

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
index.html, style.css, app.js   web interface, canvas, and hammer controls
worker.js                      Python runtime and background computation
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
