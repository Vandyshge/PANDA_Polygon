# PANDA Polygon

Interactive analytical tools for confined liquid morphologies in polygonal nanopores.

The site contains four linked explorers:

- stability diagrams for a triangular pore;
- analytical density profiles for a triangular pore;
- stability diagrams for a square pore;
- analytical density profiles for a square pore.

The diagrams compare interfacial areas of all morphologies implemented in the model. Density-profile pages expose every morphology for which a `rho` function is implemented. The original NumPy equations are evaluated directly in the browser through Pyodide, so no server-side calculation or data upload is required.

## Website

GitHub Pages serves the contents of [`docs/`](docs/). Open `docs/index.html` through a local HTTP server for development; opening it directly as a `file://` URL will not work because the model modules are fetched at runtime.

## Structure

```text
docs/
  index.html                 triangular stability diagram
  triangle-profiles.html     triangular density profiles
  square-stability.html      square stability diagram
  square-profiles.html       square density profiles
  model/                     analytical NumPy implementations
```

## License

Scientific use and redistribution terms will be added before a tagged release.
