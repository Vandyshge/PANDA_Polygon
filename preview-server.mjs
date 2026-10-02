import { createServer } from "node:http";
import { readFile, stat } from "node:fs/promises";
import { extname, join, normalize } from "node:path";

const root = new URL("./docs/", import.meta.url).pathname.replace(/^\/(.:\/)/, "$1");
const mime = { ".html": "text/html; charset=utf-8", ".js": "text/javascript; charset=utf-8", ".css": "text/css; charset=utf-8", ".py": "text/x-python; charset=utf-8", ".png": "image/png" };
createServer(async (request, response) => {
  try {
    const relative = decodeURIComponent(new URL(request.url, "http://localhost").pathname).replace(/^\/+/, "") || "index.html";
    const path = normalize(join(root, relative));
    if (!path.startsWith(normalize(root)) || !(await stat(path)).isFile()) throw new Error("not found");
    response.writeHead(200, { "Content-Type": mime[extname(path)] || "application/octet-stream", "Cache-Control": "no-store" });
    response.end(await readFile(path));
  } catch { response.writeHead(404); response.end("Not found"); }
}).listen(4174, "127.0.0.1", () => console.log("PANDA Polygon preview: http://127.0.0.1:4174/"));
