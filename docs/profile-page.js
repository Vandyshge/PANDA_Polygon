import { computeProfile, listMorphologies } from "./model-runtime.js?v=20261005";

const geometry = document.body.dataset.geometry;
const DEFAULTS = { theta: 180, delta: 0.017, length: 1, ly: 1, phi: 0.25 };
const COLORS = { triangle: "#4e8c78", square: "#557aa4" };
const LABELS = {
  filament_w0c0: "Free filament", filament_w0c0_straight: "Straight free filament",
  filament_w1c0: "One-wall filament", filament_w2c0: "Two-wall filament",
  filament_w2c0_opposite: "Opposite-wall filament", filament_w2c1: "Two-wall corner filament",
  filament_w3c0: "Three-wall filament", filament_w3c1: "Three-wall one-corner filament",
  filament_w3c2: "Three-wall two-corner filament", filament_w4c0: "Four-wall filament",
  filament_w4c1: "Four-wall one-corner filament", filament_w4c2: "Four-wall two-corner filament",
  filament_w4c3: "Four-wall three-corner filament", droplet_w1c0: "Wall droplet",
  droplet_w2c0: "Two-wall droplet", droplet_w2c1: "Corner droplet",
  capsule_w3c0: "Three-wall capsule", capsule_w4c0: "Four-wall capsule",
  capsule_w4c4: "Four-corner capsule", perforation: "Perforation",
};

const state = { ...DEFAULTS, morphology: "", result: null, requestId: 0 };
const els = Object.fromEntries([...document.querySelectorAll("[id]")].map((node) => [node.id.replace(/-([a-z])/g, (_, c) => c.toUpperCase()), node]));
let plot = null;
let inputTimer = 0;
let resizeFrame = 0;

function labelFor(name) { return LABELS[name] || name.replaceAll("_", " "); }

function readInputs() {
  state.morphology = els.morphology.value;
  state.theta = Number(els.theta.value);
  state.delta = Number(els.delta.value);
  state.length = Number(els.length.value);
  state.ly = Number(els.ly.value);
  state.phi = Number(els.phi.value);
  els.thetaOutput.value = `${state.theta.toFixed(0)}°`;
  els.deltaOutput.value = state.delta.toFixed(3);
  els.phiOutput.value = state.phi.toFixed(3);
  els.lengthOutput.value = state.length.toFixed(2);
  els.profileTitle.textContent = labelFor(state.morphology);
  els.geometryLabel.textContent = geometry === "triangle" ? "triangular pore" : "square pore";
  els.parameterLine.textContent = `morphology=${state.morphology} · θ=${state.theta.toFixed(1)}° · δ=${state.delta.toFixed(3)} · l=${state.length.toFixed(3)} · lᵧ=${state.ly.toFixed(3)} · φ=${state.phi.toFixed(3)}`;
}

async function calculate() {
  readInputs();
  const id = ++state.requestId;
  els.loading.hidden = false;
  els.loading.textContent = "Calculating density profile…";
  try {
    const result = await computeProfile({ geometry, morphology: state.morphology, thetaDeg: state.theta, delta: state.delta, length: state.length, ly: state.ly, phi: state.phi });
    if (id !== state.requestId) return;
    state.result = result;
    els.existenceStatus.classList.toggle("invalid", !result.exists);
    els.existenceStatus.textContent = result.exists
      ? "The selected morphology exists for these parameters."
      : "This morphology does not satisfy its geometrical existence condition for these parameters.";
    renderExistenceRanges(result);
    draw();
    els.loading.hidden = true;
  } catch (error) {
    els.loading.textContent = "The analytical profile could not be evaluated.";
    console.error(error);
  }
}

function renderExistenceRanges(result) {
  const controls = {
    theta: { input: els.theta, note: els.thetaRegion, digits: 0, suffix: "°" },
    phi: { input: els.phi, note: els.phiRegion, digits: 3, suffix: "" },
    delta: { input: els.delta, note: els.deltaRegion, digits: 3, suffix: "" },
    length: { input: els.length, note: els.lengthRegion, digits: 2, suffix: "" },
  };
  for (const [axis, control] of Object.entries(controls)) {
    const domain = result.domains[axis];
    const ranges = result.ranges[axis] || [];
    control.input.style.setProperty("--existence-track", gradientForRanges(ranges, domain));
    control.note.textContent = ranges.length
      ? `Existence region: ${ranges.map(([start, end]) => `${start.toFixed(control.digits)}–${end.toFixed(control.digits)}${control.suffix}`).join(", ")}`
      : "No existence region for the current parameters.";
    control.note.classList.toggle("empty", !ranges.length);
  }
}

function gradientForRanges(ranges, [minimum, maximum]) {
  const inactive = "#dce1e5";
  const active = "#79a985";
  if (!ranges.length) return inactive;
  const stops = [`${inactive} 0%`];
  for (const [start, end] of ranges) {
    const from = Math.max(0, Math.min(100, (start - minimum) / (maximum - minimum) * 100));
    const to = Math.max(0, Math.min(100, (end - minimum) / (maximum - minimum) * 100));
    stops.push(`${inactive} ${from}%`, `${active} ${from}%`, `${active} ${to}%`, `${inactive} ${to}%`);
  }
  stops.push(`${inactive} 100%`);
  return `linear-gradient(to right, ${stops.join(", ")})`;
}

function scheduleCompute() { clearTimeout(inputTimer); inputTimer = setTimeout(calculate, 140); }

function canvasContext() {
  const rect = els.plotWrap.getBoundingClientRect();
  const dpr = Math.min(devicePixelRatio || 1, 2);
  const width = Math.max(320, Math.round(rect.width));
  const height = Math.max(340, Math.round(rect.height));
  els.profileCanvas.width = width * dpr;
  els.profileCanvas.height = height * dpr;
  const ctx = els.profileCanvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, width, height };
}

function draw() {
  if (!state.result) return;
  const { ctx, width, height } = canvasContext();
  const margin = { left: 74, right: 25, top: 30, bottom: 62 };
  plot = { x: margin.left, y: margin.top, w: width - margin.left - margin.right, h: height - margin.top - margin.bottom };
  ctx.clearRect(0, 0, width, height);
  ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, width, height);
  const xs = state.result.z;
  const ys = state.result.rho;
  // Physical pore walls are fixed at delta = 0. The wetting-layer thickness
  // changes the analytical profile, not the displayed pore geometry.
  const wallPositions = geometry === "square"
    ? [-0.5, 0.5]
    : [-Math.sqrt(3) / 6, Math.sqrt(3) / 3];
  const wallSpan = wallPositions[1] - wallPositions[0];
  const xMin = wallPositions[0] - 0.06 * wallSpan, xMax = wallPositions[1] + 0.06 * wallSpan;
  const finite = ys.filter(Number.isFinite);
  const yMax = Math.max(1, ...finite) * 1.08;
  const xToPx = (x) => plot.x + (x - xMin) / (xMax - xMin) * plot.w;
  const yToPx = (y) => plot.y + plot.h - Math.max(0, y) / yMax * plot.h;

  ctx.strokeStyle = "rgba(25,39,45,.14)"; ctx.lineWidth = 1;
  ctx.font = "14px Arial, Helvetica, sans-serif"; ctx.fillStyle = "#18272d";
  for (let i = 0; i <= 5; i++) {
    const f = i / 5, x = plot.x + f * plot.w, value = xMin + f * (xMax - xMin);
    ctx.beginPath(); ctx.moveTo(x, plot.y); ctx.lineTo(x, plot.y + plot.h); ctx.stroke();
    ctx.textAlign = "center"; ctx.textBaseline = "top"; ctx.fillText(value.toFixed(2), x, plot.y + plot.h + 10);
  }
  for (let i = 0; i <= 5; i++) {
    const f = i / 5, y = plot.y + plot.h - f * plot.h;
    ctx.beginPath(); ctx.moveTo(plot.x, y); ctx.lineTo(plot.x + plot.w, y); ctx.stroke();
    ctx.textAlign = "right"; ctx.textBaseline = "middle"; ctx.fillText((f * yMax).toFixed(2), plot.x - 10, y);
  }

  ctx.save(); ctx.setLineDash([5, 5]); ctx.strokeStyle = "#7d837c";
  wallPositions.forEach((value) => { const x = xToPx(value); ctx.beginPath(); ctx.moveTo(x, plot.y); ctx.lineTo(x, plot.y + plot.h); ctx.stroke(); });
  ctx.restore();

  ctx.beginPath();
  let drawing = false;
  xs.forEach((x, index) => {
    const insidePhysicalPore = x >= wallPositions[0] && x <= wallPositions[1];
    if (!insidePhysicalPore || !Number.isFinite(ys[index])) { drawing = false; return; }
    const px = xToPx(x), py = yToPx(ys[index]);
    drawing ? ctx.lineTo(px, py) : ctx.moveTo(px, py);
    drawing = true;
  });
  ctx.strokeStyle = COLORS[geometry]; ctx.lineWidth = 3; ctx.lineJoin = "round"; ctx.lineCap = "round"; ctx.stroke();
  ctx.strokeStyle = "#26332d"; ctx.lineWidth = 1; ctx.strokeRect(plot.x, plot.y, plot.w, plot.h);
  ctx.font = "italic 16px Arial, Helvetica, sans-serif"; ctx.fillStyle = "#18272d"; ctx.textAlign = "center"; ctx.textBaseline = "bottom";
  ctx.fillText("z/a", plot.x + plot.w / 2, height - 8);
  drawDensityAxisLabel(ctx, 20, plot.y + plot.h / 2);
}

function drawDensityAxisLabel(ctx, x, y) {
  ctx.save();
  ctx.translate(x, y);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = "left";
  ctx.fillStyle = "#18272d";
  ctx.font = "italic 16px Arial, Helvetica, sans-serif";
  const main = "ρ/ρ";
  const mainWidth = ctx.measureText(main).width;
  ctx.font = "italic 11px Arial, Helvetica, sans-serif";
  const subWidth = ctx.measureText("bulk").width;
  const start = -(mainWidth + subWidth) / 2;
  ctx.font = "italic 16px Arial, Helvetica, sans-serif";
  ctx.textBaseline = "alphabetic";
  ctx.fillText(main, start, 0);
  ctx.font = "italic 11px Arial, Helvetica, sans-serif";
  ctx.textBaseline = "top";
  ctx.fillText("bulk", start + mainWidth, 1);
  ctx.restore();
}

function reset() {
  els.theta.value = DEFAULTS.theta; els.delta.value = DEFAULTS.delta; els.length.value = DEFAULTS.length;
  els.ly.value = DEFAULTS.ly; els.phi.value = DEFAULTS.phi;
  els.morphology.value = els.morphology.options[0].value;
  calculate();
}

function exportPng() {
  const link = document.createElement("a");
  link.download = `${geometry}_${state.morphology}_density_profile.png`;
  link.href = els.profileCanvas.toDataURL("image/png"); link.click();
}

async function init() {
  try {
    const names = await listMorphologies(geometry, "RHO");
    names.forEach((name) => els.morphology.add(new Option(`${labelFor(name)}  ·  ${name}`, name)));
    els.morphology.value = names[0];
    await calculate();
  } catch (error) {
    els.loading.textContent = "The analytical model could not be loaded. Check the internet connection and reload the page.";
    console.error(error);
  }
}

[els.theta, els.delta, els.phi, els.length].forEach((input) => input.addEventListener("input", scheduleCompute));
[els.morphology, els.ly].forEach((input) => input.addEventListener("change", scheduleCompute));
els.reset.addEventListener("click", reset);
els.export.addEventListener("click", exportPng);
new ResizeObserver(() => { cancelAnimationFrame(resizeFrame); resizeFrame = requestAnimationFrame(draw); }).observe(els.plotWrap);
init();
