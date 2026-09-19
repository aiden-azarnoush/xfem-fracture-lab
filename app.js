"use strict";
// 2D Fracture Lab. The specimen is a 300 x 200 mm plate, clamped on its left
// edge and struck by a hammer on its right edge at the aimed height.
// The cracks come from a dynamic bond-based peridynamics model: the plate is a
// cloud of nodes bonded to every neighbour within a horizon of three spacings;
// a bond breaks permanently once its stretch exceeds the critical value
// s0 = sqrt(4 pi Gc / (9 E delta)) set by the fracture energy. Broken bonds are
// the cracks; nothing about their path, turning, or branching is prescribed.
const $ = (id) => document.getElementById(id);
const canvas = $("canvas"), ctx = canvas.getContext("2d");
const LX = 0.3, LY = 0.2, HORIZON = 3, JITTER = 0.15;
const meshes = { coarse: 5e-3, medium: 3.5e-3, fine: 2.5e-3 };          // node spacing, m
const presets = { glass: [70, 2500, 10], ceramic: [200, 3900, 35], polymer: [3, 1200, 300] }; // E GPa, rho kg/m3, Gc J/m2
let mesh = "medium", aim = 0.5, P = null, busy = false, strikes = 0, fractured = false, flash = 0, solveStart = 0;

// ------------------------------------------------------------------ plate
function buildPlate() {
  const dx = meshes[mesh], E = +$("young").value * 1e9, rho = +$("density").value, Gc = +$("toughness").value, scatter = +$("structure").value;
  const nx = Math.round(LX / dx), ny = Math.round(LY / dx), N = nx * ny;
  const x = new Float64Array(N), y = new Float64Array(N);
  let seed = 20260918;
  const rnd = () => { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 4294967296; };
  const gauss = () => Math.sqrt(-2 * Math.log(Math.max(rnd(), 1e-12))) * Math.cos(2 * Math.PI * rnd());
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const a = j * nx + i; x[a] = (i + 0.5) * dx + JITTER * dx * (2 * rnd() - 1); y[a] = (j + 0.5) * dx + JITTER * dx * (2 * rnd() - 1); }
  const delta = HORIZON * dx;
  const bx = Math.ceil(LX / delta) + 1, by = Math.ceil(LY / delta) + 1, head = new Int32Array(bx * by).fill(-1), next = new Int32Array(N);
  const cx = (a) => Math.min(bx - 1, Math.max(0, Math.floor(x[a] / delta))), cy = (a) => Math.min(by - 1, Math.max(0, Math.floor(y[a] / delta)));
  for (let a = 0; a < N; a++) { const c = cy(a) * bx + cx(a); next[a] = head[c]; head[c] = a; }
  const I = [], J = [], L0 = [];
  for (let a = 0; a < N; a++) for (let ox = -1; ox <= 1; ox++) for (let oy = -1; oy <= 1; oy++) {
    const gx = cx(a) + ox, gy = cy(a) + oy; if (gx < 0 || gy < 0 || gx >= bx || gy >= by) continue;
    for (let b = head[gy * bx + gx]; b !== -1; b = next[b]) {
      if (b === a) continue;
      const d = Math.hypot(x[b] - x[a], y[b] - y[a]);
      if (d < delta && d > 0.5 * dx) { I.push(a); J.push(b); L0.push(d); }
    }
  }
  const M = I.length, thick = 1;
  const c = (9 * E) / (Math.PI * thick * delta ** 3), s0 = Math.sqrt((4 * Math.PI * Gc) / (9 * E * delta));
  const strength = new Float64Array(N), fixed = new Uint8Array(N);
  for (let a = 0; a < N; a++) {
    strength[a] = s0 * Math.exp(scatter * gauss());
    const edge = Math.min(y[a], LY - y[a]);
    if (edge < 2 * delta) strength[a] *= 1 + 9 * (1 - edge / (2 * delta));   // no-fail zone at the free top/bottom edges
    if (x[a] < 2 * delta) fixed[a] = 1;                                        // clamped left edge
  }
  const Ib = Int32Array.from(I), Jb = Int32Array.from(J), L0b = Float64Array.from(L0), s0b = new Float64Array(M), alive = new Uint8Array(M).fill(1);
  for (let k = 0; k < M; k++) s0b[k] = 0.5 * (strength[Ib[k]] + strength[Jb[k]]);
  const cw = Math.sqrt(E / rho), dt = (0.4 * dx) / cw;
  P = { N, M, dx, delta, x, y, I: Ib, J: Jb, L0: L0b, s0b, alive, c, vol: dx * dx * thick, rho, thick, dt, s0, strength, fixed,
        ux: new Float64Array(N), uy: new Float64Array(N), vx: new Float64Array(N), vy: new Float64Array(N),
        fx: new Float64Array(N), fy: new Float64Array(N), cnt: new Float64Array(N), dmg: new Float64Array(N), stretch: new Float64Array(N), hits: [] };
  for (let k = 0; k < M; k++) P.cnt[Ib[k]]++;
  $("meshCount").textContent = `${N.toLocaleString()} nodes · ${(dx * 1e3).toFixed(1)} mm`;
}

// ------------------------------------------------------------------ time stepping
function step() {
  const { N, M, I, J, L0, s0b, alive, c, vol, rho, thick, dt, x, y, ux, uy, vx, vy, fx, fy, fixed } = P;
  for (const H of P.hits) {
    if (H.step >= H.total) continue;
    const amp = H.push * Math.min(1, H.step / H.ramp);            // the hammer head dents the edge inward
    for (let a = 0; a < N; a++) { if (Math.hypot(x[a] - LX, y[a] - H.hy) < H.r) { ux[a] = -amp; uy[a] = 0; vx[a] = 0; vy[a] = 0; } }
    H.step++;
  }
  fx.fill(0); fy.fill(0);
  for (let k = 0; k < M; k++) {
    if (!alive[k]) continue;
    const a = I[k], b = J[k];
    const rx = x[b] + ux[b] - x[a] - ux[a], ry = y[b] + uy[b] - y[a] - uy[a];
    const Ln = Math.hypot(rx, ry), s = (Ln - L0[k]) / L0[k];
    if (s >= s0b[k]) { alive[k] = 0; P.dmg[a]++; continue; }
    const f = (c * s * vol) / Ln;
    fx[a] += f * rx; fy[a] += f * ry;
  }
  const inv = dt / (rho * thick);
  for (let a = 0; a < N; a++) {
    if (fixed[a]) { ux[a] = uy[a] = vx[a] = vy[a] = 0; continue; }
    vx[a] = (vx[a] + fx[a] * inv) * 0.998; vy[a] = (vy[a] + fy[a] * inv) * 0.998; ux[a] += dt * vx[a]; uy[a] += dt * vy[a];
  }
}
function stretchField() {           // mean |stretch| of live bonds per node, for the strain map
  const { N, M, I, J, L0, alive, x, y, ux, uy, stretch, cnt } = P; stretch.fill(0);
  for (let k = 0; k < M; k++) { if (!alive[k]) continue; const a = I[k], b = J[k];
    const Ln = Math.hypot(x[b] + ux[b] - x[a] - ux[a], y[b] + uy[b] - y[a] - uy[a]); stretch[a] += Math.abs(Ln - L0[k]) / L0[k]; }
  for (let a = 0; a < N; a++) stretch[a] /= cnt[a];
}
function crackStats() {             // clusters of broken bonds -> number of cracks, longest crack, breakthrough
  const { N, M, I, J, alive, x, delta } = P;
  const parent = Int32Array.from({ length: N }, (_, i) => i);
  const find = (a) => { while (parent[a] !== a) { parent[a] = parent[parent[a]]; a = parent[a]; } return a; };
  const nodeBroken = new Int32Array(N);
  for (let k = 0; k < M; k++) if (!alive[k]) nodeBroken[I[k]]++;
  for (let k = 0; k < M; k++) if (!alive[k]) { const ra = find(I[k]), rb = find(J[k]); if (ra !== rb) parent[ra] = rb; }
  const size = new Map(), minx = new Map(), maxx = new Map();
  for (let a = 0; a < N; a++) if (nodeBroken[a] > 2) { const r = find(a); size.set(r, (size.get(r) || 0) + 1);
    minx.set(r, Math.min(minx.has(r) ? minx.get(r) : 1, x[a])); maxx.set(r, Math.max(maxx.has(r) ? maxx.get(r) : 0, x[a])); }
  let cracks = 0, longest = 0, through = false;
  for (const [r, n] of size) if (n >= 6) { cracks++; longest = Math.max(longest, maxx.get(r) - minx.get(r)); if (minx.get(r) < 2.5 * delta) through = true; }
  let broken = 0; for (let k = 0; k < M; k++) if (!alive[k]) broken++;
  return { cracks, longest, through, brokenFrac: broken / M };
}

// ------------------------------------------------------------------ drawing
function frame() {
  const b = canvas.getBoundingClientRect(), scale = Math.min((b.width - 115) / 300, (b.height - 78) / 200);
  return { w: b.width, h: b.height, x: (b.width - 300 * scale) / 2 - 16, y: (b.height - 200 * scale) / 2 - 4, s: Math.max(0.1, scale) };
}
function draw() {
  const b = canvas.getBoundingClientRect(), dpr = devicePixelRatio || 1;
  if (canvas.width !== Math.round(b.width * dpr) || canvas.height !== Math.round(b.height * dpr)) { canvas.width = Math.round(b.width * dpr); canvas.height = Math.round(b.height * dpr); }
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, b.width, b.height);
  const f = frame(), w = 300 * f.s, h = 200 * f.s, x = f.x, y = f.y, sx = w / LX, sy = h / LY;
  ctx.fillStyle = "#dbe5f3"; ctx.fillRect(x, y, w, h);
  if (P) {
    const view = $("fieldView").value, px = Math.max(1.5, P.dx * sx * 0.95);
    if (view === "stress") stretchField();
    let vmax = 1e-9; if (view === "stress") for (let a = 0; a < P.N; a++) vmax = Math.max(vmax, P.stretch[a]);
    for (let a = 0; a < P.N; a++) {
      let t;
      if (view === "stress") { t = Math.min(1, P.stretch[a] / vmax); if (t < 0.03) continue; }
      else t = Math.min(1, Math.max(0, (P.strength[a] / P.s0 - 0.7) / 0.9));
      ctx.fillStyle = view === "stress" ? `hsl(${214 - t * 22} ${40 + t * 30}% ${89 - t * 25}%)` : `hsl(214 30% ${92 - t * 18}%)`;
      ctx.fillRect(x + (P.x[a] + P.ux[a]) * sx - px / 2, y + h - (P.y[a] + P.uy[a]) * sy - px / 2, px, px);
    }
    if ($("showMesh").checked) {
      ctx.fillStyle = "#9eb1c8";
      for (let a = 0; a < P.N; a++) ctx.fillRect(x + P.x[a] * sx - 0.6, y + h - P.y[a] * sy - 0.6, 1.2, 1.2);
    }
    // cracks: broken nearest-neighbour bonds drawn as short segments, white halo then red
    const { I, J, alive, L0, dx } = P;
    for (const [lw, col] of [[3.2, "#fff"], [1.6, "#bc3a2f"]]) {
      ctx.lineWidth = lw; ctx.strokeStyle = col; ctx.lineCap = "round"; ctx.beginPath();
      for (let k = 0; k < P.M; k++) {
        if (alive[k] || I[k] > J[k] || L0[k] > 1.6 * dx) continue;
        const a = I[k], bb = J[k];
        ctx.moveTo(x + (P.x[a] + P.ux[a]) * sx, y + h - (P.y[a] + P.uy[a]) * sy);
        ctx.lineTo(x + (P.x[bb] + P.ux[bb]) * sx, y + h - (P.y[bb] + P.uy[bb]) * sy);
      }
      ctx.stroke();
    }
  }
  ctx.strokeStyle = "#7192b5"; ctx.lineWidth = 1; ctx.strokeRect(x, y, w, h);
  ctx.fillStyle = "#526b88"; ctx.fillRect(x - 8, y, 8, h);
  ctx.strokeStyle = "#9eb1c8"; ctx.beginPath();
  for (let j = 0; j < h; j += 10) { ctx.moveTo(x - 17, y + j + 7); ctx.lineTo(x - 8, y + j); }
  ctx.stroke();
  ctx.fillStyle = "#8193a8"; ctx.font = "11px system-ui"; ctx.textAlign = "center";
  ctx.fillText("300 mm", x + w / 2, y - 16);
  ctx.save(); ctx.translate(x - 26, y + h / 2); ctx.rotate(-Math.PI / 2); ctx.fillText("200 mm · fixed edge", 0, 0); ctx.restore();
  const hy = y + h * (1 - aim);
  ctx.setLineDash([3, 4]); ctx.strokeStyle = "#3a6ac476"; ctx.beginPath(); ctx.moveTo(x + w - 28, hy); ctx.lineTo(x + w + 30, hy); ctx.stroke(); ctx.setLineDash([]);
  ctx.save();
  ctx.translate(x + w + 36 - Math.sin(flash * Math.PI) * 16, hy); ctx.rotate(-0.22);
  ctx.fillStyle = "#b68b57"; ctx.fillRect(1, -3, 8, 48);
  ctx.fillStyle = "#284b72"; ctx.beginPath(); ctx.roundRect(-12, -12, 34, 22, 4); ctx.fill();
  ctx.fillStyle = "#8da9c5"; ctx.fillRect(-12, -8, 5, 14);
  ctx.restore();
  if (fractured) {
    ctx.fillStyle = "#fffef1ec"; ctx.fillRect(x + 8, y + h / 2 - 20, w - 16, 40);
    ctx.fillStyle = "#a45c23"; ctx.textAlign = "center"; ctx.font = "600 14px system-ui";
    ctx.fillText("Cracked through to the fixed edge", x + w / 2, y + h / 2 + 5);
  }
  if (flash > 0) { flash = Math.max(0, flash - 0.06); requestAnimationFrame(draw); }
}
new ResizeObserver(draw).observe($("stage"));

// ------------------------------------------------------------------ interface
function setStatus(m) { $("status").textContent = m; }
function controls() {
  document.querySelectorAll("[data-mesh],#material,#structure,#hammer,#height,#position,#reset").forEach((el) => (el.disabled = busy));
  ["young", "density", "toughness"].forEach((id) => ($(id).disabled = busy || $("material").value !== "custom"));
  $("strike").disabled = busy || fractured;
  $("strike").textContent = busy ? "Cracking…" : fractured ? "Specimen broken" : "Strike specimen";
}
function readouts() {
  const s = crackStats();
  $("strikes").textContent = strikes;
  $("cracks").textContent = s.cracks;
  $("length").innerHTML = `${Math.round(s.longest * 1000)} <small>mm</small>`;
  $("ratio").innerHTML = strikes ? `${(100 * s.brokenFrac).toFixed(2)} <small>%</small>` : "—";
  return s;
}
function reset() {
  if (busy) return;
  buildPlate(); strikes = 0; fractured = false;
  $("solveTime").textContent = "";
  $("growthSummary").textContent = "Aim the hammer at any height and strike. Cracks run, turn, and branch on their own.";
  readouts(); controls(); updateLegend(); draw();
  setStatus("Fresh specimen. Aim the hammer and strike.");
}
function hammerValues() {
  aim = +$("position").value / 100;
  $("positionValue").value = `${Math.round(aim * 100)}%`;
  $("heightValue").value = `${$("height").value} cm`;
  $("energy").textContent = ((+$("hammer").value * 9.81 * +$("height").value) / 100).toFixed(3) + " J";
  draw();
}
function strike() {
  if (busy || fractured) return;
  for (const id of ["young", "density", "toughness"]) if (!$(id).checkValidity()) { $(id).reportValidity(); return; }
  const energy = (+$("hammer").value * 9.81 * +$("height").value) / 100;
  const push = Math.min(110e-6, Math.max(8e-6, 40e-6 * Math.sqrt(energy / 0.589)));   // dent depth grows with impact energy
  P.hits.push({ hy: aim * LY, r: 0.010, push, step: 0, ramp: 400, total: 2400 });
  strikes++; busy = true; flash = 1; solveStart = performance.now();
  controls(); draw();
  setStatus(`Strike ${strikes}: ${energy.toFixed(2)} J at ${Math.round(aim * 100)} % height. Cracks running…`);
  requestAnimationFrame(loop);
}
function loop() {
  const active = P.hits.some((H) => H.step < H.total);
  const start = performance.now();
  while (active && performance.now() - start < 40) for (let i = 0; i < 5; i++) step();
  const s = readouts(); draw();
  if (active) { requestAnimationFrame(loop); return; }
  busy = false;
  $("solveTime").textContent = `${((performance.now() - solveStart) / 1000).toFixed(1)} s`;
  if (s.through) { fractured = true; setStatus("The crack reached the fixed edge. Reset for a fresh specimen."); }
  else if (s.cracks === 0) setStatus("No crack: the blow was too light for this material. Raise the drop height or the mass.");
  else setStatus(`${s.cracks} crack${s.cracks > 1 ? "s" : ""}, longest ${Math.round(s.longest * 1000)} mm. Strike again to grow them.`);
  controls();
}
function updateLegend() {
  $("stressLegend").hidden = false;
  $("stressLegend").textContent = $("fieldView").value === "stress" ? "Bond stretch: light = relaxed · dark = highly strained"
    : +$("structure").value > 0.05 ? "Strength: light = weaker · dark = stronger" : "Nearly uniform strength";
}
$("strike").onclick = strike;
$("reset").onclick = reset;
$("material").onchange = () => { const p = presets[$("material").value]; if (p) ["young", "density", "toughness"].forEach((id, i) => ($(id).value = p[i])); reset(); };
["young", "density", "toughness"].forEach((id) => $(id).addEventListener("change", reset));
document.querySelectorAll("[data-mesh]").forEach((button) => {
  button.setAttribute("aria-pressed", button.dataset.mesh === mesh);
  button.onclick = () => { if (busy) return; mesh = button.dataset.mesh;
    document.querySelectorAll("[data-mesh]").forEach((b) => { b.classList.toggle("active", b === button); b.setAttribute("aria-pressed", b === button); });
    reset(); };
});
$("showMesh").onchange = draw;
["hammer", "height", "position"].forEach((id) => ($(id).oninput = hammerValues));
function aimAt(e) {
  const f = frame(), b = canvas.getBoundingClientRect();
  aim = Math.max(0.08, Math.min(0.92, 1 - (e.clientY - b.top - f.y) / (200 * f.s)));
  $("position").value = aim * 100; $("positionValue").value = `${Math.round(aim * 100)}%`; draw();
}
canvas.onpointermove = (e) => { if (!busy && !fractured) aimAt(e); };
canvas.onpointerdown = (e) => { if (e.button !== 0) return; e.preventDefault(); canvas.focus(); if (!busy && !fractured) { aimAt(e); strike(); } };
canvas.onkeydown = (e) => {
  if (e.code === "Space" || e.code === "Enter") { e.preventDefault(); strike(); }
  else if (["ArrowUp", "ArrowDown"].includes(e.code)) { e.preventDefault(); $("position").value = Math.max(8, Math.min(92, +$("position").value + (e.code === "ArrowUp" ? 2 : -2))); hammerValues(); }
};
$("settingsButton").onclick = () => { const open = document.body.classList.toggle("settings-open"); $("settingsButton").setAttribute("aria-expanded", open); $("settingsButton").textContent = open ? "Close settings" : "Settings"; };
document.addEventListener("keydown", (e) => { if (e.key === "Escape" && document.body.classList.contains("settings-open")) $("settingsButton").click(); });
$("methodButton").onclick = () => $("method").showModal();
$("closeMethod").onclick = () => $("method").close();
$("structure").onchange = reset;
$("fieldView").onchange = () => { updateLegend(); draw(); };
$("loading").hidden = true;
hammerValues();
reset();
