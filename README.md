# Two independent OrcaSlicer crash reproducers

Both failures are verified in the official Linux x86_64 nightly built from
upstream main commit `d81a212f53d15ad841c88df9e71325e37f0efe14`.
The application log identifies the build as `2.5.0-dev d81a212f`.

## 1. Slicing crash: empty layer vector with a raft

Open **slice-crash.3mf** as a project, retaining its settings, then click
**Slice plate**. The file loads successfully; slicing crashes in
`groupingVolumesForBrim()` when it reads `layers.front()->lslices` on an empty
layer vector. Expected: a readable no-layers error.

Geometry is one synthetic, closed tetrahedron: four vertices, four triangles,
0.0002 × 0.0002 × 0.002 mm. The project retains complete profile settings so
GUI import preserves the same slicing setup. Layer height and first-layer
height are 0.2 mm; raft layers are 2; brim and supports are disabled.
No original user geometry is included.

Official 2.4.2 and 2.5.0 beta also crash when this project is sliced through the CLI.

Controls:
- **no-raft-control.3mf**: same geometry, zero raft layers. CLI reports no layers.
- **normal-tetrahedron-control.3mf**: 10 × 10 × 1 mm tetrahedron, otherwise the
  same two-raft-layer setup. CLI slicing succeeds.

CLI: `orca-slicer --slice 1 --outputdir output slice-crash.3mf`

See `slice-backtrace.txt` and `slice-cli-verification.json`.

## 2. Loading crash: missing filament_colour setting

Open **load-crash.3mf** as a project. Orca crashes before slicing, in
`Slic3r::GUI::Plater::priv::load_files()`.

The project settings deliberately omit `filament_colour`. The loader assumes
that option exists and dereferences it without checking:

```cpp
project_filament_count = config_loaded.option<ConfigOptionStrings>("filament_colour")->size();
```

Source at the tested commit:
https://github.com/OrcaSlicer/OrcaSlicer/blob/d81a212f53d15ad841c88df9e71325e37f0efe14/src/slic3r/GUI/Plater.cpp#L9449

Expected: default the missing option or reject the project with a readable error.

**load-ok-control.3mf** adds only `filament_colour: ["#F2754E"]` to the settings.
That control loads successfully in the same GUI.

This load failure was independently verified on the current nightly GUI.
No claim is made here about reproducing the load failure on older Linux releases.
The reporter also observed the load failure on Windows ARM64.

See `load-backtrace.txt` and `gui-verification.json`.

## Regeneration

`python make_reproducers.py` regenerates all five projects using only Python's
standard library. Keep `slice-settings.json` alongside the script.

These are synthetic test cases, not printable models. No fork source changes
are needed to reproduce either crash.
