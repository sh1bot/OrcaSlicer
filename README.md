# OrcaSlicer KD-tree rebuild demonstrations

Open either .3mf **as a project**, retain its embedded settings, and Slice plate.
Both projects select Bambu Lab A1 / 0.4 mm, PLA, 0.20 mm layers and 35% infill.
They have zero perimeters, zero top/bottom shell layers, no anchors, no brim,
no supports, no prime tower and disabled adaptive density. These are slicing
workloads intended to expose ordering cost, rather than finished printable parts.

## Stock-pattern case (recommended for the standalone KD-tree PR)

`kdtree-stock-perforated-plate.3mf`

- Ordinary rectilinear infill; no Voronoi feature required.
- One connected 190 x 190 x 2 mm plate containing 2,500 square through-holes.
- Holes have deterministic position jitter to reduce regular-grid ties.
- Fit is inside the A1 build area; object centred at X=128, Y=128 mm.
- The mesh is manifold: every edge has exactly two incident triangles, with
  opposite directions; enclosed mesh volume is 52,200 mm³.

Verification uses the mesh read back from the packaged 3MF. Orca's native mesh
slicer produces one island with 2,500 holes at three tested heights. Native
rectilinear infill at 35%, 0.45 mm width and 0.20 mm height produces 6,436 or
6,454 paths depending on layer direction.

The original before/after Release ordering executables were run five times,
alternating execution order. Each invocation reports the best of three calls;
the numbers below are medians of those five reports:

| Direction | Paths | Before | With rebuild | Ordering time reduction |
|---|---:|---:|---:|---:|
| First | 6,436 | 140.9 ms | 90.4 ms | 35.9% |
| Alternating | 6,454 | 160.5 ms | 88.4 ms | 45.0% |

Both tested directions produce identical routes and travel between the two
versions. Other data can have different equal-distance tie choices.
These measurements exercise the isolated greedy endpoint-ordering routine using
paths generated from this mesh, not complete project import, slicing or G-code
export. Full slice-time improvement will be smaller and depend on the machine.
The exact ordering output/timing is not promised across versions or platforms.

## Original geometry-backed comparison

`kdtree-voronoi-50mm-cube.3mf`

- 50 x 50 x 50 mm cube, uniform Random Voronoi at 35%, angle 45 degrees.
- Needs a build containing the Voronoi feature (such as whole-stack).
- Corner smoothing, relaxation, wall snapping, shell tracing and void cutoff
  are disabled.
- Matches the geometry and nominal extrusion settings of the earlier
  geometry-backed benchmark. That benchmark sampled ten heights and tested
  generated path endpoints; it was not a full Orca GUI slicing benchmark.
  Subsequent Voronoi changes mean current paths need not match the old dumps.

## Comparing application builds

1. Use equivalent Release builds with and without only the KD-tree change.
2. Open the same project and preserve its settings.
3. Time slicing through completion, including G-code generation. Repeat with
   freshly loaded projects to avoid cached work; alternate the builds.
4. Compare measured elapsed time. The preview's estimated PRINT time is not
   the performance measurement for this CPU-time optimization.

STLs are supplied as geometry-only fallbacks. Importing an STL does NOT carry
these settings, so use the .3mf project for reproducibility.

Both 3MF archives, their relationships/XML/settings and their closed mesh
geometry were checked. Native slicing and infill generation were checked locally;
a full Orca GUI import or complete Print/G-code run was not available here.
