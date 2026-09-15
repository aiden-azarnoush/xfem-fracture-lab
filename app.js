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
  previousLengths = [],
  growthDuration = 900;
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
  if (field && $("fieldView").value === "stress") {
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
  if ($("fieldView").value === "toughness" || !field) {
    const variation = +$("structure").value;
    const dx = w / 60,
      dy = h / 40;
    for (let ix = 0; ix < 60; ix++)
      for (let iy = 0; iy < 40; iy++) {
        const px = ((ix + 0.5) * 0.3) / 60,
          py = ((iy + 0.5) * 0.2) / 40;
        const v = Math.exp(
          variation *
            (0.6 * Math.sin(80 * px + 23 * py) +
              0.4 * Math.sin(44 * px - 95 * py + 1.7)),
        );
        const t = Math.max(0, Math.min(1, (v - 0.6) / 1.0));
        ctx.fillStyle = `hsl(${210 - t * 10} ${35 + t * 10}% ${92 - t * 20}%)`;
        ctx.fillRect(x + ix * dx, y + h - (iy + 1) * dy, dx + 0.5, dy + 0.5);
      }
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
  const growth = Math.min(
    1,
    (performance.now() - growthStart) / growthDuration,
  );
  crackData.forEach((c, i) => {
    const shownLength =
      (previousLengths[i] || 0) +
      (c.length - (previousLengths[i] || 0)) * growth;
    let remaining = shownLength;
    const origin = c.points[0];
    let tip = origin;
    ctx.beginPath();
    ctx.moveTo(x + origin[0] * 1000 * f.s, y + (200 - origin[1] * 1000) * f.s);
    for (let j = 1; j < c.points.length; j++) {
      const a = c.points[j - 1],
        b = c.points[j],
        distance = Math.hypot(b[0] - a[0], b[1] - a[1]);
      const amount = Math.min(1, remaining / distance);
      tip = [a[0] + (b[0] - a[0]) * amount, a[1] + (b[1] - a[1]) * amount];
      ctx.lineTo(x + tip[0] * 1000 * f.s, y + (200 - tip[1] * 1000) * f.s);
      remaining -= distance;
      if (remaining <= 0) break;
    }
    ctx.lineJoin = "round";
    ctx.lineCap = "round";
    ctx.lineWidth = 4;
    ctx.strokeStyle = "#fff";
    ctx.stroke();
    ctx.lineWidth = 1.8;
    ctx.strokeStyle = "#bc3a2f";
    ctx.stroke();
    ctx.fillStyle = c.stopped ? "#854831" : "#d64a36";
    ctx.beginPath();
    ctx.arc(
      x + tip[0] * 1000 * f.s,
      y + (200 - tip[1] * 1000) * f.s,
      3.4,
      0,
      Math.PI * 2,
    );
    ctx.fill();
    ctx.font = "11px system-ui";
    ctx.textAlign = "left";
    ctx.fillText(
      `C${i + 1}`,
      x + w + 5,
      y + (200 - origin[1] * 1000) * f.s + 4,
    );
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
    ctx.fillText("Fracture limit reached", x + w / 2, y + h / 2 + 5);
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
    .querySelectorAll(
      "[data-mesh],#material,#structure,#hammer,#height,#position,#reset",
    )
    .forEach((el) => (el.disabled = busy));
  ["young", "poisson", "toughness"].forEach(
    (id) => ($(id).disabled = busy || $("material").value !== "custom"),
  );
  $("strike").disabled = busy || !ready || fractured;
  $("strike").textContent = busy
    ? "Calculating…"
    : fractured
      ? "Fracture limit reached"
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
  updateLegend();
  strikes = 0;
  fractured = false;
  $("ratio").textContent = "—";
  $("growthSummary").textContent =
    "Cracks choose a direction at each step. Move the hammer to change the loading.";
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
  setStatus("Comparing growth directions at the crack tips…");
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
      variation: +$("structure").value,
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
  worker = new Worker("worker.js?v=2");
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
    growthDuration = Math.max(600, (r.events?.length || 1) * 400);
    crackData = r.cracks;
    field = r.field;
    fractured = r.fractured;
    $("ratio").textContent = r.ratio.toFixed(2) + "×";
    $("solveTime").textContent = `Solved in ${r.seconds.toFixed(1)} s`;
    updateLegend();
    $("growthSummary").textContent = r.events?.length
      ? r.events
          .map(
            (e) =>
              `C${e.crack}: ${e.new ? "started " : ""}${e.turn === 0 ? "continued" : `turned ${Math.abs(e.turn)}°`}`,
          )
          .join(" · ")
      : "No direction exceeded the fracture resistance.";
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
updateLegend();
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

function updateLegend() {
  $("stressLegend").hidden = false;
  if (field && $("fieldView").value === "stress")
    $("stressLegend").textContent =
      `Shear stress · 0–${Math.max(...field.values).toFixed(1)} MPa`;
  else
    $("stressLegend").textContent = +$("structure").value
      ? "Toughness: light = weaker · dark = stronger"
      : "Uniform fracture toughness";
}
$("structure").onchange = reset;
$("fieldView").onchange = () => {
  updateLegend();
  draw();
};
