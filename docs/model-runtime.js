const PYODIDE_VERSION = "v314.0.7";
const PYODIDE_BASE = `https://cdn.jsdelivr.net/pyodide/${PYODIDE_VERSION}/full/`;

let runtimePromise = null;
const loadedModels = new Set();

function setRuntimeState(text, kind = "") {
  const node = document.querySelector("#runtime-state");
  if (!node) return;
  node.textContent = text;
  node.className = `runtime-state ${kind}`.trim();
}

export async function getRuntime(geometry) {
  if (!runtimePromise) {
    setRuntimeState("loading analytical model");
    runtimePromise = (async () => {
      const pyodide = await loadPyodide({ indexURL: PYODIDE_BASE });
      await pyodide.loadPackage("numpy");
      return pyodide;
    })().catch((error) => {
      setRuntimeState("model failed to load", "error");
      throw error;
    });
  }

  const pyodide = await runtimePromise;
  if (!loadedModels.has(geometry)) {
    const source = await fetch(`./model/${geometry}.py`).then((response) => {
      if (!response.ok) throw new Error(`Cannot load ${geometry}.py`);
      return response.text();
    });
    pyodide.FS.writeFile(`${geometry}_model.py`, source);
    await pyodide.runPythonAsync(`import ${geometry}_model`);
    loadedModels.add(geometry);
  }
  setRuntimeState("analytical model ready", "ready");
  return pyodide;
}

function setGlobals(pyodide, values) {
  for (const [key, value] of Object.entries(values)) pyodide.globals.set(key, value);
}

export async function listMorphologies(geometry, registry = "AREA") {
  const pyodide = await getRuntime(geometry);
  const registryName = `${registry}_${geometry.toUpperCase()}`;
  const result = await pyodide.runPythonAsync(`list(${geometry}_model.${registryName}.keys())`);
  const values = result.toJs();
  result.destroy();
  return [...values];
}

export async function computePhase({ geometry, thetaDeg, delta, lMax, resolution, enabled }) {
  const pyodide = await getRuntime(geometry);
  setGlobals(pyodide, { __theta_deg: thetaDeg, __delta: delta, __l_max: lMax, __resolution: resolution, __enabled: enabled });
  const upper = geometry.toUpperCase();
  const json = await pyodide.runPythonAsync(`
import json, numpy as np
_module = ${geometry}_model
_registry = _module.AREA_${upper}
_names = [str(name) for name in __enabled if str(name) in _registry]
_n = int(__resolution)
_phi = np.linspace(0.001, 1.0, _n)
_length = np.linspace(0.001, float(__l_max), _n)
_PHI, _L = np.meshgrid(_phi, _length)
_theta = np.deg2rad(float(__theta_deg))
_areas = []
for _name in _names:
    _area = np.asarray(_registry[_name](_L, _PHI, _theta, float(__delta), False), dtype=float)
    _area[~np.isfinite(_area)] = np.inf
    _areas.append(_area)
if _areas:
    _all = np.stack(_areas)
    _valid = np.isfinite(_all).any(axis=0)
    _phase = np.where(_valid, np.argmin(_all, axis=0), -1).astype(np.int16)
else:
    _phase = np.full((_n, _n), -1, dtype=np.int16)
_phi_accessible = max(0.0, min(1.0, (1.0 - (2.0*np.sqrt(3.0) if '${geometry}' == 'triangle' else 2.0)*float(__delta))**2))
json.dumps({'names': _names, 'phase': _phase.tolist(), 'phiAccessible': _phi_accessible})
  `);
  return JSON.parse(json);
}

export async function computeProfile({ geometry, morphology, thetaDeg, delta, length, ly, phi, samples = 520 }) {
  const pyodide = await getRuntime(geometry);
  setGlobals(pyodide, { __morphology: morphology, __theta_deg: thetaDeg, __delta: delta, __length: length, __ly: ly, __phi: phi, __samples: samples });
  const upper = geometry.toUpperCase();
  const json = await pyodide.runPythonAsync(`
import json, numpy as np
_module = ${geometry}_model
_rho_registry = _module.RHO_${upper}
_cond_registry = _module.COND_${upper}
_name = str(__morphology)
_theta = np.deg2rad(float(__theta_deg))
_z_min, _z_max = (-0.65, 0.65) if '${geometry}' == 'square' else (-0.58, 0.72)
_z = np.linspace(_z_min, _z_max, int(__samples))
_exists = bool(np.asarray(_cond_registry[_name](float(__length), float(__phi), _theta, float(__delta))).item())
_rho = np.asarray(_rho_registry[_name](_z, float(__length), float(__ly), float(__phi), _theta, float(__delta)), dtype=float)
_rho[~np.isfinite(_rho)] = 0.0
json.dumps({'z': _z.tolist(), 'rho': _rho.tolist(), 'exists': _exists})
  `);
  return JSON.parse(json);
}
