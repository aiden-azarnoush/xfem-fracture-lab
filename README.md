# 2D Fracture Lab

Interactive, single-screen **Mode III XFEM teaching experiment** by Aiden Azarnoush. Choose a material, select a coarse/medium/fine triangular mesh, aim the hammer vertically, and click to create or extend cracks. Multiple impact locations produce multiple straight cracks. The footer follows the author's other mechanics websites.

## Run

Serve the repository as a static website (GitHub Pages works):

```sh
python3 -m http.server 8000
```

Open `http://localhost:8000`. Do not open `index.html` directly with a `file:` URL, because web workers and Python-source fetching require an HTTP origin. An internet connection is needed to download the pinned Pyodide 0.27.7 runtime and NumPy from jsDelivr. Computation runs locally in a web worker; there is no Python server or account requirement.

The default glass specimen, 0.5 kg hammer, and 12 cm drop initiate a crack on all three meshes. Strike again to extend it. Aim sufficiently far from an existing crack to start another. If a strike does not cause growth, increase the drop height or hammer mass. Space/Enter strikes and up/down arrows aim when the specimen canvas has keyboard focus. A position slider and Strike button provide alternatives to pointer interaction.

Changing material, custom properties, or mesh resets the specimen. Changing hammer mass, height, or position preserves cracks. Reset starts a new experiment. The page itself does not scroll; on small screens the Settings button opens a panel that can scroll independently. Strike and Reset remain outside that panel.

## Files

- `index.html`, `style.css`: interface and footer.
- `app.js`: canvas, hammer, interactions, and worker messages.
- `worker.js`: loads Python/NumPy and runs calculations away from the interface thread.
- `python/xfem.py`: numerical model and crack-growth decision.
- `python/test_xfem.py`: numerical consistency and behavior tests.
- `python/requirements.txt`: dependency for running the solver/tests in local Python.

## Numerical formulation

The 300 × 200 × 5 mm rectangle has a fixed left edge. The unknown is scalar out-of-plane displacement `w(x,y)`, with `mu = E / [2(1+nu)]`. The antiplane stiffness is:

`K_ij = thickness × integral(mu × grad(Phi_i) · grad(Phi_j) dA)`.

The approximation combines linear triangle shape functions with **shifted Heaviside** enrichment behind each crack tip and the **Mode III tip function** `sqrt(r/h) sin(theta/2)` within two element widths of a tip. The local tip axis points left, along the prescribed direction of growth. Extra degrees of freedom are assembled into the actual stiffness matrix and solved; cracks are not merely drawn on an unchanged FEM solution.

Integration polygons are split along crack planes and through tip coordinates, triangulated, and evaluated with three-point triangle quadrature. Triangles near a tip receive an additional subdivision. These are integration subcells; the underlying mesh is fixed. The meshes contain 192, 432, or 768 triangles. The dense system is diagonally scaled, stabilized with `1e-11` on the scaled diagonal, and checked for a relative free-DOF residual below `1e-4`.

Two smooth, opposing sets of nodal loads on the right edge idealize a tearing impulse. Each set has unit resultant before scaling. For a drop energy `E_drop = m g h_drop`, the peak-equivalent load magnitude for each set is `F = E_drop / 0.0005 m`. This assumed stopping distance is an illustrative load calibration, not a contact calculation.

Each strike solves the existing geometry and a trial crack extension under the same force. For load-controlled linear elasticity, released potential energy is the increase in strain energy:

`G ≈ (U_trial − U_current) / (thickness × delta_a)`.

Growth is accepted when `G >= Gc` and the drop energy is sufficient to create the incremental crack area. A strike advances at most one increment (`0.65 h`); a newly nucleated trial notch is `1.3 h` long. The UI shows `G/Gc` from that trial, including unsuccessful trials. The heat map shows element-averaged shear-stress magnitude from the last accepted geometry under the current strike load, in MPa. Its color scale is normalized per strike.

## Scope and limitations

This is an educational, quasi-static antiplane model, **not a validated dynamic hammer/contact or general plane-stress fracture solver**.

- Cracks begin at the right edge and follow prescribed straight horizontal paths. There is no branching or curved path selection.
- Finite-notch nucleation, enrichment radius, quadrature, and extension increments depend on the mesh. No convergence claim is made.
- Impacts near an existing crack (within `1.5 h`) target that crack; there is no separate closely spaced-crack nucleation law.
- At most 12 separate cracks are supported. Interaction enters through the global elastic solve, without a special intersection/coalescence model.
- A remaining ligament of one element width or less is classified as separated; the displayed terminal crack is completed to the left boundary. The stress field remains that of the last solved, nearly separated specimen. This final visual completion is not a resolved last-ligament failure calculation.
- There is no plasticity, fatigue accumulation, crack-face contact, time integration, rebound, or moving-fragment calculation. Repeating a subcritical strike does not accumulate damage.
- The glass, ceramic, and polymer properties are illustrative elastic presets. The polymer preset does not model real viscoelastic/plastic polymer fracture.
- Crack width is exaggerated for visibility. Dimensions of the undeformed specimen remain fixed.

## Tests

With NumPy installed:

```sh
python3 -m unittest discover -s python -v
```

Checks cover integration-area conservation, analytic tip gradients against finite differences, elastic material scaling, equilibrium residuals, intact reflection symmetry, rejection of invalid inputs, subcritical strikes, initiation, repeated growth, multiple cracks, and terminal separation on all meshes. These checks establish numerical consistency and intended behavior; they are not experimental validation.

## References

- [GetFEM: level sets and XFEM](https://getfem.org/userdoc/xfem.html) — enrichment and cut-cell integration concepts.
- [University of Liège: Mode III energy release rate](https://www.ltas-cm3.ulg.ac.be/FractureMechanics/print.php?p=Lecture3_P4) — antiplane crack-tip field and energy concepts.
- [Pyodide web workers](https://pyodide.org/en/stable/usage/webworker.html) — running Python outside the interface thread.
