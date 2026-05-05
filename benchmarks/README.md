# Benchmarks

Benchmarks should make numerical claims reproducible. Each benchmark reports the computed integral, absolute error against an analytic reference value, runtime, and the quadrature configuration used.

List available benchmark cases:

```bash
python -m benchmarks.run_benchmarks --list
```

Run the default quick suite:

```bash
python -m benchmarks.run_benchmarks
```

Run the baseline suite and emit JSON:

```bash
python -m benchmarks.run_benchmarks --suite baseline --format json
```

The legacy sphere-area entry point still works:

```bash
python benchmarks/benchmark_sphere_area.py
```
