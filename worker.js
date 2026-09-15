"use strict";
let runtime;
async function initialize() {
  postMessage({ type: "loading", message: "Downloading the Python runtime…" });
  importScripts("https://cdn.jsdelivr.net/pyodide/v0.27.7/full/pyodide.js");
  runtime = await loadPyodide({
    indexURL: "https://cdn.jsdelivr.net/pyodide/v0.27.7/full/",
  });
  postMessage({ type: "loading", message: "Loading numerical libraries…" });
  await runtime.loadPackage("numpy");
  const response = await fetch("./python/xfem.py?v=2");
  if (!response.ok)
    throw new Error("The Python solver could not be downloaded.");
  await runtime.runPythonAsync(await response.text());
  postMessage({ type: "ready" });
}
const initialized = initialize().catch((error) => {
  postMessage({ type: "fatal", message: String(error.message || error) });
  throw error;
});
onmessage = async ({ data }) => {
  try {
    await initialized;
    runtime.globals.set("input_json", JSON.stringify(data.payload));
    const result = await runtime.runPythonAsync("run_json(input_json)");
    postMessage({ type: "result", id: data.id, result: JSON.parse(result) });
  } catch (error) {
    postMessage({
      type: "error",
      id: data.id,
      message: String(error.message || error),
    });
  }
};
