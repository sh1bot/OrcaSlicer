# Shared infill material density

Voronoi and Smooth blend Gyroid, Schwarz D and Fischer-Koch consume an
object-wide material-density field independent of point sampling and routing.
Other TPMS adaptive modes use their own density handling.

Surface density is the sparse infill setting; interior density is an absolute
percentage. Linear, Quadratic and Exponential grade density using the selected adaptive
spatial profile. Distance modes use wall depth normalised to the deepest point
of each connected body; Lobes and Normal X/Y/Z use lobe or section centers. Equal surface
and interior values give a uniform base. Adaptive density defaults to Disabled.
Voronoi supports every adaptive mode. Distance warp, Smooth blend and Stepped
shells give the same continuous wall-depth point density. Lobes and Normal modes
blend density contributions from neighbouring lobes or sections. None warp the
Voronoi cloud coordinates. Corner reinforcement remains independent of
this switch.

Corner reinforcement compares a Gaussian-smoothed
physical signed distance field with its original values. The difference boosts both convex and concave corners and
may exceed surface density, capped at 100%. Sigma is in millimetres, defaults to zero (disabled), and is limited to 64 mm.
The distance grid coarsens to limit memory on large objects. When its spacing
approaches sigma, the sampled Gaussian and corner reinforcement become less
accurate; a small sigma may no longer resolve useful detail.

Density queries evaluate the selected upstream profile and gradient directly.
With nonzero sigma, OpenCV filters a float copy of physical signed distances;
queries subtract the interpolated smoothed distance from the raw distance to
calculate reinforcement. Zero sigma allocates no smoothing buffer.
Immutable geometry and density evaluators are shared across layer workers.
Construction and OpenCV filtering report progress and support cancellation.

Voronoi uses density to select points. Smooth blend TPMS converts density through
each pattern's frequency calibration and blends fixed lattice levels, extending
beyond the surface/interior range for reinforcement. Spatial variation does not
directly rescale lattice coordinates. Zero-density regions produce no infill.

With zero sigma and no fused modifiers, Smooth blend uses TPMS depth directly.

## Fused modifiers

Compositing material density uses a float grid only when a participating modifier
is present. Modifier targets are evaluated from their own settings and blended
without density encoding or quantization. Queries interpolate the resulting
material fractions.

Material density blend width is available on modifiers with a Voronoi or
Smooth blend TPMS/Gyroid parent. Negative values retain independent print
regions. Zero replaces density at field resolution; positive values blend
entirely inside the modifier with smoothstep, reaching its target at the
specified inward distance.

Surface density, adaptive mode, interior density, gradient and corner smoothing remain editable;
other settings inherit the parent while fusion is enabled. Modifiers use the
parent's physical distance and adaptive spatial fields. Their own sides, floor and
ceiling determine blend distance. Fields blend successively in Orca's modifier
order without splitting the parent's infill into separate regions.

## Distance calculation

Adaptive infill, corner reinforcement and fused material density share the
object's adaptive voxel geometry. The existing raster and distance transform
supply interior distances. Consumers needing signed distances additionally
transform the exterior, preserving the interior costs used by adaptive grading.
Modifier blend masks use the same voxel calculation.

The grid starts at 0.5 mm spacing and coarsens to stay within about one million
nodes. Narrow features and gaps smaller than a voxel may disappear; snapping,
blend boundaries and corner reinforcement have this resolution limit too.
When grid spacing approaches Gaussian sigma, smoothing becomes less accurate.
Normal X/Y/Z measure within their sections; physical corner reinforcement uses
three-dimensional distances. Density settings reuse the object's immutable grid.
