# Adaptive infill distances

`AdaptiveInfillDistanceField` retains the voxel distances used by
`AdaptiveInfillField`. It exposes physical distances and interpolation weights
for consumers needing a sampled distance field. Signed sampling is negative
inside the object and positive outside.

Slices provide Orca's processed outlines, including floors, roofs and holes.
The existing raster and separable distance transform compute interior costs.
An optional exterior transform supplies signed samples without changing the
interior costs or the existing adaptive depth, body and lobe profiles.

The grid starts at 0.5 mm spacing and coarsens to stay within about one million
nodes. Narrow features smaller than a voxel may disappear. Physical samples are
capped at 1000 mm. Construction reports progress and supports cancellation.
The default interior-only calculation avoids the exterior transform.
