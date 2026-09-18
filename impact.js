"use strict";
// 2D Fracture Lab · Impact: bond-based peridynamics, plane stress.
// The plate is 300 x 200 mm. Nodes sit on a jittered grid; every pair closer
// than the horizon (3 spacings) is a bond. A bond breaks permanently when its
// stretch exceeds s0 = sqrt(4 pi Gc / (9 E delta)). Explicit time stepping.
const $ = (id) => document.getElementById(id);
const canvas = $("canvas"), ctx = canvas.getContext("2d");
const LX = 0.3, LY = 0.2, HORIZON = 3, JITTER = 0.15;
const presets = { glass: [70, 2500, 10], ceramic: [200, 3900, 35], polymer: [3, 1200, 300] }; // GPa, kg/m3, J/m2
const spacing = { coarse: 5e-3, medium: 3.5e-3, fine: 2.5e-3 };
let mesh = "medium", P = null, running = false, lastHit = null, strikes = 0, t0 = 0;

// ---------------------------------------------------------------- plate
function buildPlate() {
  const dx = spacing[mesh], E = +$("young").value * 1e9, rho = +$("density").value, Gc = +$("toughness").value, scatter = +$("structure").value;
  const nx = Math.round(LX / dx), ny = Math.round(LY / dx), N = nx * ny;
  const x = new Float64Array(N), y = new Float64Array(N);
  let seed = 20260918;
  const rnd = () => { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 4294967296; };
  const gauss = () => Math.sqrt(-2 * Math.log(Math.max(rnd(), 1e-12))) * Math.cos(2 * Math.PI * rnd());
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const a = j * nx + i; x[a] = (i + 0.5) * dx + JITTER * dx * (2 * rnd() - 1); y[a] = (j + 0.5) * dx + JITTER * dx * (2 * rnd() - 1); }
  const delta = HORIZON * dx;
  const bx = Math.ceil(LX / delta) + 1, by = Math.ceil(LY / delta) + 1, head = new Int32Array(bx * by).fill(-1), next = new Int32Array(N);
  const cell = (a) => Math.min(by - 1, Math.max(0, Math.floor(y[a] / delta))) * bx + Math.min(bx - 1, Math.max(0, Math.floor(x[a] / delta)));
  for (let a = 0; a < N; a++) { const c = cell(a); next[a] = head[c]; head[c] = a; }
  const I = [], J = [], L0 = [];
  for (let a = 0; a < N; a++) {
    const cx = Math.min(bx - 1, Math.max(0, Math.floor(x[a] / delta))), cy = Math.min(by - 1, Math.max(0, Math.floor(y[a] / delta)));
    for (let ox = -1; ox <= 1; ox++) for (let oy = -1; oy <= 1; oy++) {
      const gx = cx + ox, gy = cy + oy; if (gx < 0 || gy < 0 || gx >= bx || gy >= by) continue;
      for (let b = head[gy * bx + gx]; b !== -1; b = next[b]) {
        if (b === a) continue;
        const d = Math.hypot(x[b] - x[a], y[b] - y[a]);
        if (d < delta && d > 0.5 * dx) { I.push(a); J.push(b); L0.push(d); }
      }
    }
  }
  const M = I.length, thick = 1;
  const c = (9 * E) / (Math.PI * thick * delta ** 3), s0 = Math.sqrt((4 * Math.PI * Gc) / (9 * E * delta));
  const strength = new Float64Array(N);
  for (let a = 0; a < N; a++) {
    strength[a] = s0 * Math.exp(scatter * gauss());
    const edge = Math.min(x[a], LX - x[a], y[a], LY - y[a]);
    if (edge < 2 * delta) strength[a] *= 1 + 9 * (1 - edge / (2 * delta)); // tapered no-fail zone at free edges
  }
  const Ib = Int32Array.from(I), Jb = Int32Array.from(J), L0b = Float64Array.from(L0), s0b = new Float64Array(M), alive = new Uint8Array(M).fill(1);
  for (let k = 0; k < M; k++) s0b[k] = 0.5 * (strength[Ib[k]] + strength[Jb[k]]);
  const cw = Math.sqrt(E / rho), dt = (0.4 * dx) / cw;
  P = { N, M, dx, x, y, I: Ib, J: Jb, L0: L0b, s0b, alive, c, vol: dx * dx * thick, rho, thick, dt, s0, cw,
        ux: new Float64Array(N), uy: new Float64Array(N), vx: new Float64Array(N), vy: new Float64Array(N),
        fx: new Float64Array(N), fy: new Float64Array(N), cnt: new Float64Array(N), dmg: new Float64Array(N), hits: [] };
  for (let k = 0; k < M; k++) P.cnt[Ib[k]]++;
  strikes = 0; lastHit = null; running = false;
  $("meshCount").textContent = `${N.toLocaleString()} nodes · ${(M / 2).toLocaleString()} bonds · ${(dx * 1e3).toFixed(1)} mm spacing`;
  $("s0").innerHTML = `${(s0 * 1e6).toFixed(0)} <small>µε</small>`;
  $("strikes").textContent = "0"; $("broken").innerHTML = "0 <small>%</small>"; $("area").innerHTML = "0 <small>%</small>";
  setStatus("Fresh plate. Click on it to strike.");
  draw();
}

// ---------------------------------------------------------------- hammer and time stepping
function hammerValues() {
  $("heightValue").value = `${$("height").value} cm`;
  $("radiusValue").value = `${$("radius").value} mm`;
  $("energy").textContent = ((+$("hammer").value * 9.81 * +$("height").value) / 100).toFixed(3) + " J";
}
function strike(hx, hy) {
  const energy = (+$("hammer").value * 9.81 * +$("height").value) / 100;
  const push = Math.min(60e-6, Math.max(4e-6, 20e-6 * Math.sqrt(energy / 1.0))); // dent push-out, grows with energy
  P.hits.push({ hx, hy, r: +$("radius").value * 1e-3, push, step: 0, ramp: 400, total: 2600 });
  lastHit = { hx, hy }; strikes++; $("strikes").textContent = strikes;
  setStatus(`Strike ${strikes}: ${energy.toFixed(2)} J at (${(hx * 1e3).toFixed(0)}, ${(hy * 1e3).toFixed(0)}) mm. Cracking…`);
  t0 = performance.now();
  if (!running) { running = true; requestAnimationFrame(loop); }
}
function step() {
  const { N, M, I, J, L0, s0b, alive, c, vol, rho, thick, dt, x, y, ux, uy, vx, vy, fx, fy } = P;
  for (const H of P.hits) {
    if (H.step >= H.total) continue;
    const amp = H.push * Math.min(1, H.step / H.ramp);
    for (let a = 0; a < N; a++) {
      const dxh = x[a] - H.hx, dyh = y[a] - H.hy, rr = Math.hypot(dxh, dyh);
      if (rr < H.r) { const s = amp / Math.max(rr, 1e-9); ux[a] = s * dxh; uy[a] = s * dyh; vx[a] = 0; vy[a] = 0; }
    }
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
  for (let a = 0; a < N; a++) { vx[a] = (vx[a] + fx[a] * inv) * 0.998; vy[a] = (vy[a] + fy[a] * inv) * 0.998; ux[a] += dt * vx[a]; uy[a] += dt * vy[a]; }
}
function loop() {
  if (!running) return;
  if (!P.hits.some((H) => H.step < H.total)) {
    running = false; draw();
    $("solveTime").textContent = `${((performance.now() - t0) / 1000).toFixed(1)} s`;
    setStatus(`Done. ${strikes} strike${strikes > 1 ? "s" : ""} so far; strike again or start a new plate.`);
    return;
  }
  const start = performance.now();
  while (performance.now() - start < 40) for (let i = 0; i < 5; i++) step();
  draw(); requestAnimationFrame(loop);
}

// ---------------------------------------------------------------- drawing (same frame conventions as the XFEM page)
function frame() {
  const b = canvas.getBoundingClientRect(), scale = Math.min((b.width - 115) / 300, (b.height - 78) / 200);
  return { w: b.width, h: b.height, x: (b.width - 300 * scale) / 2 - 16, y: (b.height - 200 * scale) / 2 - 4, s: Math.max(0.1, scale) };
}
function draw() {
  const b = canvas.getBoundingClientRect(), dpr = devicePixelRatio || 1;
  if (canvas.width !== Math.round(b.width * dpr) || canvas.height !== Math.round(b.height * dpr)) { canvas.width = Math.round(b.width * dpr); canvas.height = Math.round(b.height * dpr); }
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, b.width, b.height);
  const f = frame(), w = 300 * f.s, h = 200 * f.s, x0 = f.x, y0 = f.y;
  ctx.fillStyle = "#dbe5f3"; ctx.fillRect(x0, y0, w, h);
  if (P) {
    const { N, x, y, ux, uy, dmg, cnt, dx } = P, sx = w / LX, sy = h / LY, px = Math.max(1.5, dx * sx * 0.75);
    let damaged = 0;
    for (let a = 0; a < N; a++) {
      const fr = dmg[a] / cnt[a];
      if (fr < 0.12) continue;
      damaged++;
      const shade = Math.max(0, 1 - fr * 1.6);
      ctx.fillStyle = `rgb(${Math.round(40 + 120 * shade)}, ${Math.round(20 + 90 * shade)}, ${Math.round(60 + 80 * shade)})`;
      ctx.fillRect(x0 + (x[a] + ux[a]) * sx - px / 2, y0 + h - (y[a] + uy[a]) * sy - px / 2, px, px);
    }
    ctx.strokeStyle = "#8ea3bd"; ctx.lineWidth = 1; ctx.strokeRect(x0, y0, w, h);
    if (lastHit) { ctx.strokeStyle = "#d3541f"; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(x0 + lastHit.hx * sx, y0 + h - lastHit.hy * sy, +$("radius").value * 1e-3 * sx, 0, 2 * Math.PI); ctx.stroke(); }
    let broken = 0; for (let k = 0; k < P.M; k++) if (!P.alive[k]) broken++;
    $("broken").innerHTML = `${((100 * broken) / P.M).toFixed(2)} <small>%</small>`;
    $("area").innerHTML = `${((100 * damaged) / N).toFixed(1)} <small>%</small>`;
  }
  // axes labels as on the XFEM page
  ctx.fillStyle = "#5e748c"; ctx.font = "12px Inter, sans-serif"; ctx.textAlign = "center";
  ctx.fillText("300 mm", x0 + w / 2, y0 + h + 22);
  ctx.save(); ctx.translate(x0 - 22, y0 + h / 2); ctx.rotate(-Math.PI / 2); ctx.fillText("200 mm", 0, 0); ctx.restore();
}
function setStatus(t) { $("status").textContent = t; }

// ---------------------------------------------------------------- wiring
canvas.addEventListener("click", (e) => {
  const b = canvas.getBoundingClientRect(), f = frame(), w = 300 * f.s, h = 200 * f.s;
  const hx = ((e.clientX - b.left - f.x) / w) * LX, hy = (1 - (e.clientY - b.top - f.y) / h) * LY;
  if (hx < 0 || hx > LX || hy < 0 || hy > LY) return;
  strike(hx, hy);
});
$("strike").onclick = () => strike(lastHit ? lastHit.hx : LX / 2, lastHit ? lastHit.hy : LY / 2);
$("reset").onclick = buildPlate;
document.querySelectorAll("[data-mesh]").forEach((btn) => (btn.onclick = () => {
  document.querySelectorAll("[data-mesh]").forEach((b) => { b.classList.toggle("active", b === btn); b.setAttribute("aria-pressed", b === btn); });
  mesh = btn.dataset.mesh; buildPlate();
}));
$("material").onchange = () => {
  const p = presets[$("material").value], custom = !p;
  ["young", "density", "toughness"].forEach((id, i) => { $(id).disabled = !custom; if (p) $(id).value = p[i]; });
  buildPlate();
};
["young", "density", "toughness", "structure"].forEach((id) => ($(id).onchange = buildPlate));
["hammer", "height", "radius"].forEach((id) => ($(id).oninput = hammerValues));
$("settingsButton").onclick = () => {
  const open = document.body.classList.toggle("settings-open");
  $("settingsButton").setAttribute("aria-expanded", open);
  $("settingsButton").textContent = open ? "Close settings" : "Settings";
};
document.addEventListener("keydown", (e) => {
  if (e.key === "Escape" && document.body.classList.contains("settings-open")) $("settingsButton").click();
});
$("methodButton").onclick = () => $("method").showModal();
$("closeMethod").onclick = () => $("method").close();
window.addEventListener("resize", draw);
hammerValues();
buildPlate();
