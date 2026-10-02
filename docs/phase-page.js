import { computePhase, listMorphologies } from "./model-runtime.js?v=20261003";

const geometry = document.body.dataset.geometry;
const DEFAULTS = { theta: 180, delta: 0.017, lMax: 5, resolution: 180 };
const PALETTE = ["#3f8f7a", "#c69134", "#846a99", "#557aa4", "#789caf", "#aa6868", "#8e708c", "#8e7569", "#aa7694", "#aea047", "#75996a", "#95565a", "#667f6c", "#96764e", "#6e7990", "#8f8f8f"];
const LABELS = {
  filament_w0c0: "Free filament",
  filament_w0c0_straight: "Straight free filament",
  filament_w1c0: "One-wall filament",
  filament_w2c0: "Two-wall filament",
  filament_w2c0_opposite: "Opposite-wall filament",
  filament_w2c1: "Two-wall corner filament",
  filament_w3c0: "Three-wall filament",
  filament_w3c1: "Three-wall one-corner filament",
  filament_w3c2: "Three-wall two-corner filament",
  filament_w4c0: "Four-wall filament",
  filament_w4c1: "Four-wall one-corner filament",
  filament_w4c2: "Four-wall two-corner filament",
  filament_w4c3: "Four-wall three-corner filament",
  droplet_w1c0: "Wall droplet",
  droplet_w2c0: "Two-wall droplet",
  droplet_w2c1: "Corner droplet",
  droplet_w3c0: "Three-wall droplet",
  capsule_w3c0: "Three-wall capsule",
  capsule_w3c3: "Three-corner capsule",
  capsule_w4c0: "Four-wall capsule",
  capsule_w4c4: "Four-corner capsule",
  perforation: "Perforation",
};

const state = { ...DEFAULTS, morphologies: [], enabled: new Set(), colors: new Map(), result: null, selected: null, requestId: 0 };
const els = Object.fromEntries([...document.querySelectorAll("[id]")].map((node) => [node.id.replace(/-([a-z])/g, (_, c) => c.toUpperCase()), node]));
let plot = null;
let resizeFrame = 0;
let inputTimer = 0;

function labelFor(name) { return LABELS[name] || name.replaceAll("_", " "); }

function renderToggles() {
  els.morphologyToggles.replaceChildren(...state.morphologies.map((name) => {
    const label = document.createElement("label");
    label.className = "morphology-toggle";
    label.innerHTML = `<input type="checkbox" value="${name}" checked><span class="swatch" style="background:${state.colors.get(name)}"></span><span>${labelFor(name)}</span><span class="code">${name}</span>`;
    label.querySelector("input").addEventListener("change", (event) => {
      event.target.checked ? state.enabled.add(name) : state.enabled.delete(name);
      if (!state.enabled.size) { event.target.checked = true; state.enabled.add(name); }
      scheduleCompute();
    });
    return label;
  }));
  renderLegend();
}

function renderLegend() {
  const names = state.morphologies.filter((name) => state.enabled.has(name));
  els.legend.replaceChildren(...names.map((name) => {
    const item = document.createElement("div");
    item.className = "legend-item";
    item.innerHTML = `<span style="background:${state.colors.get(name)}"></span>${labelFor(name)}`;
    return item;
  }));
  els.candidateCount.textContent = `${names.length} analytical candidate${names.length === 1 ? "" : "s"}`;
}

function readInputs() {
  state.theta = Number(els.theta.value);
  state.delta = Number(els.delta.value);
  state.lMax = Math.max(1, Math.min(10, Number(els.lMax.value) || DEFAULTS.lMax));
  state.resolution = Number(els.resolution.value);
  els.thetaOutput.value = `${state.theta.toFixed(0)}°`;
  els.deltaOutput.value = state.delta.toFixed(3);
  els.diagramTheta.textContent = `θ = ${state.theta.toFixed(0)}°`;
  els.diagramDelta.textContent = `δ = ${state.delta.toFixed(3)}`;
}

async function calculate() {
  readInputs();
  renderLegend();
  const id = ++state.requestId;
  els.loading.hidden = false;
  els.loading.textContent = "Calculating stability map…";
  try {
    const result = await computePhase({ geometry, thetaDeg: state.theta, delta: state.delta, lMax: state.lMax, resolution: state.resolution, enabled: [...state.enabled] });
    if (id !== state.requestId) return;
    state.result = result;
    draw();
  } catch (error) {
    els.loading.textContent = "The analytical model could not be evaluated.";
    console.error(error);
    return;
  }
  els.loading.hidden = true;
}

function scheduleCompute() {
  clearTimeout(inputTimer);
  inputTimer = setTimeout(calculate, 140);
}

function canvasContext() {
  const rect = els.plotWrap.getBoundingClientRect();
  const dpr = Math.min(devicePixelRatio || 1, 2);
  const width = Math.max(320, Math.round(rect.width));
  const height = Math.max(340, Math.round(rect.height));
  els.phaseCanvas.width = width * dpr;
  els.phaseCanvas.height = height * dpr;
  const ctx = els.phaseCanvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, width, height };
}

function draw() {
  if (!state.result) return;
  const { ctx, width, height } = canvasContext();
  const margin = { left: 66, right: 24, top: 24, bottom: 58 };
  plot = { x: margin.left, y: margin.top, w: width - margin.left - margin.right, h: height - margin.top - margin.bottom };
  ctx.clearRect(0, 0, width, height);
  ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, width, height);

  const n = state.result.phase.length;
  const offscreen = document.createElement("canvas");
  offscreen.width = n; offscreen.height = n;
  const off = offscreen.getContext("2d");
  const image = off.createImageData(n, n);
  for (let iy = 0; iy < n; iy++) {
    const sourceRow = state.result.phase[n - 1 - iy];
    for (let ix = 0; ix < n; ix++) {
      const phaseIndex = sourceRow[ix];
      const name = phaseIndex >= 0 ? state.result.names[phaseIndex] : null;
      const rgb = name ? hexToRgb(state.colors.get(name)) : [237, 240, 237];
      const offset = (iy * n + ix) * 4;
      image.data[offset] = rgb[0]; image.data[offset + 1] = rgb[1]; image.data[offset + 2] = rgb[2]; image.data[offset + 3] = 255;
    }
  }
  off.putImageData(image, 0, 0);
  ctx.imageSmoothingEnabled = false;
  ctx.drawImage(offscreen, plot.x, plot.y, plot.w, plot.h);
  if (state.result.phiAccessible < 1) drawHatching(ctx, plot.x + state.result.phiAccessible * plot.w, plot.y, (1 - state.result.phiAccessible) * plot.w, plot.h);
  drawAxes(ctx);
  if (state.selected) drawMarker(ctx, state.selected.phi, state.selected.l);
}

function drawAxes(ctx) {
  ctx.fillStyle = "#18272d"; ctx.strokeStyle = "rgba(25,39,45,.2)"; ctx.lineWidth = 1;
  ctx.font = "13px Arial, sans-serif"; ctx.textAlign = "center"; ctx.textBaseline = "top";
  for (let i = 0; i <= 5; i++) {
    const value = i / 5, x = plot.x + value * plot.w;
    ctx.beginPath(); ctx.moveTo(x, plot.y); ctx.lineTo(x, plot.y + plot.h); ctx.stroke();
    ctx.fillText(value.toFixed(1), x, plot.y + plot.h + 10);
  }
  ctx.textAlign = "right"; ctx.textBaseline = "middle";
  for (let i = 0; i <= 5; i++) {
    const value = state.lMax * i / 5, y = plot.y + plot.h - i / 5 * plot.h;
    ctx.beginPath(); ctx.moveTo(plot.x, y); ctx.lineTo(plot.x + plot.w, y); ctx.stroke();
    ctx.fillText(Number.isInteger(value) ? value.toFixed(0) : value.toFixed(1), plot.x - 10, y);
  }
  ctx.font = "italic 17px Georgia, serif"; ctx.textAlign = "center"; ctx.textBaseline = "bottom";
  ctx.fillText("φ", plot.x + plot.w / 2, plot.y + plot.h + 49);
  ctx.save(); ctx.translate(18, plot.y + plot.h / 2); ctx.rotate(-Math.PI / 2); ctx.fillText("l", 0, 0); ctx.restore();
  ctx.strokeStyle = "#26332d"; ctx.strokeRect(plot.x, plot.y, plot.w, plot.h);
}

function drawHatching(ctx, x, y, width, height) {
  ctx.save(); ctx.fillStyle = "rgba(255,255,255,.78)"; ctx.fillRect(x, y, width, height); ctx.beginPath(); ctx.rect(x, y, width, height); ctx.clip();
  ctx.strokeStyle = "rgba(85,95,89,.58)";
  for (let d = -height; d < width + height; d += 9) { ctx.beginPath(); ctx.moveTo(x + d, y + height); ctx.lineTo(x + d + height, y); ctx.stroke(); }
  ctx.restore();
}

function pointFromEvent(event) {
  if (!plot || !state.result) return null;
  const rect = els.phaseCanvas.getBoundingClientRect();
  const x = event.clientX - rect.left, y = event.clientY - rect.top;
  if (x < plot.x || x > plot.x + plot.w || y < plot.y || y > plot.y + plot.h) return null;
  const phi = (x - plot.x) / plot.w;
  const l = (plot.y + plot.h - y) / plot.h * state.lMax;
  const n = state.result.phase.length;
  const ix = Math.max(0, Math.min(n - 1, Math.round(phi * (n - 1))));
  const iy = Math.max(0, Math.min(n - 1, Math.round(l / state.lMax * (n - 1))));
  const phaseIndex = state.result.phase[iy][ix];
  return { phi, l, x, y, name: phaseIndex >= 0 ? state.result.names[phaseIndex] : null };
}

function showTooltip(event, pin = false) {
  const point = pointFromEvent(event);
  if (!point) { if (!pin) els.tooltip.hidden = true; return; }
  els.tooltip.innerHTML = `<b>${point.name ? labelFor(point.name) : "No valid candidate"}</b><br>φ = ${point.phi.toFixed(3)} · l = ${point.l.toFixed(3)}`;
  els.tooltip.hidden = false;
  const tw = els.tooltip.offsetWidth, th = els.tooltip.offsetHeight;
  els.tooltip.style.left = `${Math.min(point.x + 14, els.plotWrap.clientWidth - tw - 8)}px`;
  els.tooltip.style.top = `${Math.max(8, point.y - th - 12)}px`;
  if (pin) { state.selected = point; draw(); }
}

function drawMarker(ctx, phi, l) {
  const x = plot.x + phi * plot.w, y = plot.y + plot.h - l / state.lMax * plot.h;
  ctx.beginPath(); ctx.arc(x, y, 5.5, 0, 2 * Math.PI); ctx.fillStyle = "#fff"; ctx.fill(); ctx.lineWidth = 2; ctx.strokeStyle = "#17201b"; ctx.stroke();
}

function reset() {
  els.theta.value = DEFAULTS.theta; els.delta.value = DEFAULTS.delta; els.lMax.value = DEFAULTS.lMax; els.resolution.value = DEFAULTS.resolution;
  state.enabled = new Set(state.morphologies); state.selected = null;
  document.querySelectorAll(".morphology-toggle input").forEach((input) => { input.checked = true; });
  calculate();
}

function exportPng() {
  const link = document.createElement("a");
  link.download = `${geometry}_stability_theta${state.theta}_delta${state.delta.toFixed(3)}.png`;
  link.href = els.phaseCanvas.toDataURL("image/png"); link.click();
}

function registerWebMcp() {
  const context = document.modelContext;
  if (!context?.registerTool) return;
  void Promise.resolve(context.registerTool({
    name: `read_${geometry}_stability_state`, title: "Read stability diagram state", description: `Read parameters and enabled morphologies for the ${geometry} stability diagram.`,
    inputSchema: { type: "object", properties: {}, additionalProperties: false }, annotations: { readOnlyHint: true, untrustedContentHint: false },
    execute: () => ({ geometry, thetaDeg: state.theta, delta: state.delta, lMax: state.lMax, morphologies: [...state.enabled] }),
  })).catch(() => {});
}

function hexToRgb(hex) { const value = parseInt(hex.slice(1), 16); return [(value >> 16) & 255, (value >> 8) & 255, value & 255]; }

async function init() {
  try {
    state.morphologies = await listMorphologies(geometry, "AREA");
    state.morphologies.forEach((name, index) => state.colors.set(name, PALETTE[index % PALETTE.length]));
    state.enabled = new Set(state.morphologies);
    renderToggles();
    await calculate();
  } catch (error) {
    els.loading.textContent = "The analytical model could not be loaded. Check the internet connection and reload the page.";
    console.error(error);
  }
}

[els.theta, els.delta].forEach((input) => input.addEventListener("input", scheduleCompute));
[els.lMax, els.resolution].forEach((input) => input.addEventListener("change", scheduleCompute));
els.selectAll.addEventListener("click", () => { state.enabled = new Set(state.morphologies); document.querySelectorAll(".morphology-toggle input").forEach((input) => { input.checked = true; }); calculate(); });
els.reset.addEventListener("click", reset); els.export.addEventListener("click", exportPng);
els.phaseCanvas.addEventListener("pointermove", (event) => showTooltip(event));
els.phaseCanvas.addEventListener("pointerleave", () => { els.tooltip.hidden = true; });
els.phaseCanvas.addEventListener("pointerdown", (event) => showTooltip(event, true));
new ResizeObserver(() => { cancelAnimationFrame(resizeFrame); resizeFrame = requestAnimationFrame(draw); }).observe(els.plotWrap);
registerWebMcp();
init();
