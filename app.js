"use strict";
const $ = (id) => document.getElementById(id);
const canvas = $("canvas"),
  ctx = canvas.getContext("2d");
let mesh = "medium",
  aim = 0.5,
  crackData = [],
  field = null,
  busy = false,
  ready = false,
  strikes = 0,
  fractured = false,
  flash = 0,
  growthStart = 0,
  previousLengths = [];
const meshes = { coarse: [12, 8], medium: [18, 12], fine: [24, 16] };
const presets = {
  glass: [70, 0.22, 10],
  ceramic: [200, 0.22, 35],
  polymer: [3, 0.35, 300],
};
function frame() {
  const b = canvas.getBoundingClientRect(),
    scale = Math.min((b.width - 115) / 300, (b.height - 78) / 200);
  return {
    w: b.width,
    h: b.height,
    x: (b.width - 300 * scale) / 2 - 16,
    y: (b.height - 200 * scale) / 2 - 4,
    s: Math.max(0.1, scale),
  };
}
function draw() {
  const b = canvas.getBoundingClientRect(),
    dpr = devicePixelRatio || 1;
  if (
    canvas.width !== Math.round(b.width * dpr) ||
    canvas.height !== Math.round(b.height * dpr)
  ) {
    canvas.width = Math.round(b.width * dpr);
    canvas.height = Math.round(b.height * dpr);
  }
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, b.width, b.height);
  const f = frame(),
    w = 300 * f.s,
    h = 200 * f.s,
    x = f.x,
    y = f.y;
  ctx.fillStyle = "#dbe5f3";
  ctx.fillRect(x, y, w, h);
  if (field) {
    const vmax = Math.max(...field.values.map(Math.abs), 1e-10);
    field.triangles.forEach((tri, i) => {
      const t = Math.min(1, Math.abs(field.values[i]) / vmax);
      ctx.fillStyle = `hsl(${214 - t * 22} ${40 + t * 30}% ${89 - t * 25}%)`;
      ctx.beginPath();
      tri.forEach((n, j) => {
        const p = field.nodes[n];
        j
          ? ctx.lineTo(x + p[0] * 1000 * f.s, y + (200 - p[1] * 1000) * f.s)
          : ctx.moveTo(x + p[0] * 1000 * f.s, y + (200 - p[1] * 1000) * f.s);
      });
      ctx.closePath();
      ctx.fill();
    });
  }
  const [nx, ny] = meshes[mesh];
  if ($("showMesh").checked) {
    ctx.strokeStyle = "#587da535";
    ctx.lineWidth = 0.6;
    ctx.beginPath();
    for (let i = 0; i <= nx; i++) {
      ctx.moveTo(x + (i * w) / nx, y);
      ctx.lineTo(x + (i * w) / nx, y + h);
    }
    for (let j = 0; j <= ny; j++) {
      ctx.moveTo(x, y + (j * h) / ny);
      ctx.lineTo(x + w, y + (j * h) / ny);
    }
    for (let i = 0; i < nx; i++)
      for (let j = 0; j < ny; j++) {
        ctx.moveTo(x + (i * w) / nx, y + ((j + 1) * h) / ny);
        ctx.lineTo(x + ((i + 1) * w) / nx, y + (j * h) / ny);
      }
    ctx.stroke();
  }
  ctx.strokeStyle = "#7192b5";
  ctx.lineWidth = 1;
  ctx.strokeRect(x, y, w, h);
  ctx.fillStyle = "#526b88";
  ctx.fillRect(x - 8, y, 8, h);
  ctx.strokeStyle = "#9eb1c8";
  ctx.beginPath();
  for (let j = 0; j < h; j += 10) {
    ctx.moveTo(x - 17, y + j + 7);
    ctx.lineTo(x - 8, y + j);
  }
  ctx.stroke();
  ctx.fillStyle = "#8193a8";
  ctx.font = "11px system-ui";
  ctx.textAlign = "center";
  ctx.fillText("300 mm", x + w / 2, y - 16);
  ctx.save();
  ctx.translate(x - 26, y + h / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText("200 mm · fixed edge", 0, 0);
  ctx.restore();
  const growth = Math.min(1, (performance.now() - growthStart) / 280);
  crackData.forEach((c, i) => {
    const shownLength =
      (previousLengths[i] || 0) +
      (c.length - (previousLengths[i] || 0)) * growth;
    const yy = y + (200 - c.y * 1000) * f.s,
      tip = x + (300 - shownLength * 1000) * f.s;
    ctx.lineWidth = 4;
    ctx.strokeStyle = "#f8fafc";
    ctx.beginPath();
    ctx.moveTo(x + w, yy);
    ctx.lineTo(tip, yy);
    ctx.stroke();
    ctx.lineWidth = 1.5;
    ctx.strokeStyle = "#d64a36";
    ctx.stroke();
    ctx.fillStyle = "#d64a36";
    ctx.beginPath();
    ctx.arc(tip, yy, 3, 0, Math.PI * 2);
    ctx.fill();
    ctx.font = "10px system-ui";
    ctx.textAlign = "left";
    ctx.fillText(`C${i + 1}`, x + w + 5, yy + 4);
  });
  const hy = y + h * (1 - aim);
  ctx.setLineDash([3, 4]);
  ctx.strokeStyle = "#3a6ac476";
  ctx.beginPath();
  ctx.moveTo(x + w - 28, hy);
  ctx.lineTo(x + w + 30, hy);
  ctx.stroke();
  ctx.setLineDash([]);
  ctx.save();
  ctx.translate(x + w + 36 - Math.sin(flash * Math.PI) * 16, hy);
  ctx.rotate(-0.22);
  ctx.fillStyle = "#b68b57";
  ctx.fillRect(1, -3, 8, 48);
  ctx.fillStyle = "#284b72";
  ctx.beginPath();
  ctx.roundRect(-12, -12, 34, 22, 4);
  ctx.fill();
  ctx.fillStyle = "#8da9c5";
  ctx.fillRect(-12, -8, 5, 14);
  ctx.restore();
  if (fractured) {
    ctx.fillStyle = "#fffef1ec";
    ctx.fillRect(x + 8, y + h / 2 - 20, w - 16, 40);
    ctx.fillStyle = "#a45c23";
    ctx.textAlign = "center";
    ctx.font = "600 14px system-ui";
    ctx.fillText("Specimen separated", x + w / 2, y + h / 2 + 5);
  }
  if (flash > 0 || growth < 1) {
    flash = Math.max(0, flash - 0.06);
    requestAnimationFrame(draw);
  }
}
new ResizeObserver(draw).observe($("stage"));
$("methodButton").onclick = () => $("method").showModal();
$("closeMethod").onclick = () => $("method").close();
let worker,
  requestId = 0,
  timer;
function setStatus(message) {
  $("status").textContent = message;
}
function controls() {
  document
    .querySelectorAll("[data-mesh],#material,#hammer,#height,#reset")
    .forEach((el) => (el.disabled = busy));
  ["young", "poisson", "toughness"].forEach(
    (id) => ($(id).disabled = busy || $("material").value !== "custom"),
  );
  $("strike").disabled = busy || !ready || fractured;
  $("strike").textContent = busy
    ? "Calculating…"
    : fractured
      ? "Specimen separated"
      : ready
        ? "Strike specimen"
        : "Loading solver…";
}
function readouts() {
  $("strikes").textContent = strikes;
  $("cracks").textContent = crackData.length;
  $("length").innerHTML =
    `${Math.round(Math.max(0, ...crackData.map((c) => c.length)) * 1000)} <small>mm</small>`;
}
function reset() {
  if (busy) return;
  requestId++;
  crackData = [];
  previousLengths = [];
  field = null;
  $("stressLegend").hidden = true;
  strikes = 0;
  fractured = false;
  $("ratio").textContent = "—";
  $("solveTime").textContent = "";
  readouts();
  controls();
  draw();
  setStatus("Fresh specimen. Aim the hammer and strike.");
}
function hammerValues() {
  aim = +$("position").value / 100;
  $("positionValue").value = `${Math.round(aim * 100)}%`;
  $("heightValue").value = `${$("height").value} cm`;
  $("energy").textContent =
    ((+$("hammer").value * 9.81 * +$("height").value) / 100).toFixed(3) + " J";
  draw();
}
function strike() {
  if (!ready || busy || fractured) return;
  for (const id of ["young", "poisson", "toughness"]) {
    if (!$(id).checkValidity()) {
      $(id).reportValidity();
      return;
    }
  }
  busy = true;
  controls();
  flash = 1;
  draw();
  setStatus("Solving the current specimen and a virtual crack extension…");
  const id = ++requestId;
  worker.postMessage({
    id,
    payload: {
      mesh,
      young: +$("young").value,
      poisson: +$("poisson").value,
      toughness: +$("toughness").value,
      mass: +$("hammer").value,
      height: +$("height").value / 100,
      aim,
      cracks: crackData,
    },
  });
  timer = setTimeout(() => {
    if (id === requestId && busy) {
      worker.terminate();
      busy = false;
      ready = false;
      controls();
      $("loading").hidden = false;
      $("loading").textContent =
        "Calculation took too long. Reload to restart the solver.";
      setStatus("The previous crack state has been preserved.");
    }
  }, 120000);
}
function startWorker() {
  worker = new Worker("worker.js");
  worker.onmessage = ({ data }) => {
    if (data.type === "loading") {
      $("loading").textContent = data.message;
      return;
    }
    if (data.type === "ready") {
      ready = true;
      $("loading").hidden = true;
      controls();
      setStatus("Ready. Aim the hammer and click to strike.");
      return;
    }
    if (data.type === "fatal") {
      ready = false;
      $("loading").textContent =
        "Solver unavailable. Check your connection and reload.";
      controls();
      setStatus(
        "The numerical runtime could not load. No simulated cracks have been substituted.",
      );
      return;
    }
    if (data.id !== requestId) return;
    clearTimeout(timer);
    busy = false;
    if (data.type === "error") {
      setStatus(
        "Calculation failed. Reset the specimen or try a coarser mesh.",
      );
      console.error(data.message);
      controls();
      return;
    }
    const r = data.result;
    strikes++;
    previousLengths = crackData.map((c) => c.length);
    growthStart = performance.now();
    crackData = r.cracks;
    field = r.field;
    fractured = r.fractured;
    $("ratio").textContent = r.ratio.toFixed(2) + "×";
    $("solveTime").textContent = `Solved in ${r.seconds.toFixed(1)} s`;
    $("stressLegend").hidden = false;
    $("stressLegend").textContent =
      `Shear stress · 0–${Math.max(...field.values).toFixed(1)} MPa`;
    setStatus(r.message);
    readouts();
    controls();
    draw();
  };
  worker.onerror = () => {
    clearTimeout(timer);
    busy = false;
    ready = false;
    controls();
    $("loading").hidden = false;
    $("loading").textContent = "Solver unavailable. Reload to retry.";
  };
}
$("strike").onclick = strike;
$("reset").onclick = reset;
$("material").onchange = () => {
  const p = presets[$("material").value];
  if (p)
    ["young", "poisson", "toughness"].forEach((id, i) => ($(id).value = p[i]));
  reset();
};
["young", "poisson", "toughness"].forEach((id) =>
  $(id).addEventListener("change", reset),
);
document.querySelectorAll("[data-mesh]").forEach((button) => {
  button.setAttribute("aria-pressed", button.dataset.mesh === mesh);
  button.onclick = () => {
    if (busy) return;
    mesh = button.dataset.mesh;
    document.querySelectorAll("[data-mesh]").forEach((b) => {
      b.classList.toggle("active", b === button);
      b.setAttribute("aria-pressed", b === button);
    });
    const [nx, ny] = meshes[mesh];
    $("meshCount").textContent = `${nx * ny * 2} triangles`;
    reset();
  };
});
$("showMesh").onchange = draw;
["hammer", "height", "position"].forEach(
  (id) => ($(id).oninput = hammerValues),
);
function aimAt(e) {
  const f = frame(),
    b = canvas.getBoundingClientRect();
  aim = Math.max(
    0.08,
    Math.min(0.92, 1 - (e.clientY - b.top - f.y) / (200 * f.s)),
  );
  $("position").value = aim * 100;
  $("positionValue").value = `${Math.round(aim * 100)}%`;
  draw();
}
canvas.onpointermove = (e) => {
  if (!busy && !fractured) aimAt(e);
};
canvas.onpointerdown = (e) => {
  if (e.button !== 0) return;
  e.preventDefault();
  canvas.focus();
  if (!busy && !fractured) {
    aimAt(e);
    strike();
  }
};
canvas.onkeydown = (e) => {
  if (e.code === "Space" || e.code === "Enter") {
    e.preventDefault();
    strike();
  } else if (["ArrowUp", "ArrowDown"].includes(e.code)) {
    e.preventDefault();
    $("position").value = Math.max(
      8,
      Math.min(92, +$("position").value + (e.code === "ArrowUp" ? 2 : -2)),
    );
    hammerValues();
  }
};
hammerValues();
startWorker();

$("settingsButton").onclick = () => {
  const open = document.body.classList.toggle("settings-open");
  $("settingsButton").setAttribute("aria-expanded", open);
  $("settingsButton").textContent = open ? "Close settings" : "Settings";
};
document.addEventListener("keydown", (e) => {
  if (e.key === "Escape" && document.body.classList.contains("settings-open"))
    $("settingsButton").click();
});
