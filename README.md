# Crash when slicing a microscopic object with two raft layers

## Reproduce

1. Open `empty-layer-crash.3mf` as a project and retain its settings.
2. Click **Slice plate**.
3. Orca exits instead of reporting that the object has no printable layers.

The project has one synthetic tetrahedron, four vertices and four triangles.
It measures 0.0002 × 0.0002 × 0.002 mm. It is closed, consistently oriented,
and has positive volume, but is intentionally too small to print.
Layer height and first-layer height are 0.2 mm; raft layers are 2.
Supports and brim are disabled. No Voronoi or adaptive infill is involved.
No geometry from the original project is included.

CLI equivalent:

```sh
orca-slicer --slice 1 --outputdir output empty-layer-crash.3mf
```

## Expected result

Reject the object with a readable "No layers were detected" error, without crashing.

## Controls

- `no-raft-control.3mf`: identical microscopic geometry, but zero raft layers.
  Produces the expected no-layers error.
- `normal-tetrahedron-control.3mf`: a 10 × 10 × 1 mm tetrahedron with the same
  two-raft-layer configuration. Slices successfully.

`make_reproducer.py` regenerates all three files using only Python's standard library:

```sh
python make_reproducer.py
```

## Observed failure

Official Linux x86_64 AppImages: OrcaSlicer 2.4.2 and 2.5.0 beta.
The failing case exits with SIGSEGV (signal 11). See `verification.json` for
exit codes; 156 is the CLI's handled no-layers error (-100 modulo 256).
GUI reproduction of these synthetic files has not been independently tested.

The original, larger project also crashes on Windows ARM64 according to its
reporter; the synthetic files have not yet been tested there.

Backtrace from the failing synthetic project:

```text
Slic3r::groupingVolumesForBrim(...)
Slic3r::PrintObject::slice()
Slic3r::Print::process(...)
Slic3r::CLI::run(...)
main
```

`groupingVolumesForBrim()` uses `layers.front()->lslices` without checking
whether the layer vector is empty. The empty-layer guard in `PrintObject::slice()`
comes after this call. Brim grouping runs even with brim disabled.
