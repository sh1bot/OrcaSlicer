# Three-dimensional Voronoi sparse infill

The `voronoi` pattern prints layer sections of faces between nearest 3D control
points. Point generation, cell geometry, wall routing and extrusion are separate
stages. Density is statistical: clipping and floor-row rounding affect the
material deposited.

## Point cloud and geometry

`PointCloudProvider` supplies reproducible XYZ points through half-open box
queries. The random provider uses spatially seeded Poisson sampling, calibrated
to infill density and bead count. Sparse infill direction rotates the cloud
about Z through the object's bounding-box centre; model alignment adds the
object's XY rotation. Per-layer rotation templates do not rotate the 3D cloud.

Exterior points participate so that the model boundary does not distort cells.
Covering bounds determine how far the known cloud must extend. Resulting walls
and floors are clipped to the fill surface. Object-wide XY coverage lets a
layer section serve multiple fill surfaces, at extra cost on sparse islands.

Each engine retains the points that can affect its working slab. Covering
radii bound ingestion above and retirement below it. A persistent CGAL 3D
Delaunay triangulation inserts new sites and batch-removes retired clusters,
retaining unaffected connectivity. Its shallow faces supply floor polygons.
A weighted 2D power diagram supplies wall sections and covering bounds;
it is rebuilt at each height because its weights depend on vertical distance.

Workers process ascending layer batches with independent mutable engines.
Temporary fillers reuse the batch's engines; copies own independent networks.

CGAL supplies the Delaunay and power-diagram geometry and is required at build
time. Slicing cost and memory grow with point count; large, dense objects are
more expensive. Fill Multiline thickens walls without adding more cells.

## Extrusion and routing

Steep faces produce wall lines. Shallow faces produce floor bands rounded
outward to native flow spacing while keeping the supported boundary fixed.
At constant layer height, continuing bands overlap by at least half a bead
before clipping; face ends, holes and changing layer heights limit that guarantee.

Native Rectilinear fill prints alternating floor rows along the face section,
advancing from the supported side. Reachable row ends extend into previous-layer
perimeter footprints. This does not guarantee support at both ends or select
bridge flow. Orca's bead width, spacing and flow are retained.

Wall edges form trails through paired junctions. Pairing favours straight
crossings; three-way junctions cycle turns across layers using their generating
points. Closed-trail starts are hashed by layer to scatter seams. Each wall
edge is traversed once.

Orca orders wall trails and complete floor sweeps as reversible endpoint chords.
Wall trails can reverse; floors retain their supported-side-first
order, so their actual entry and exit may differ from the scheduled chord.

Fill Multiline uses Orca's native path offsets: odd counts retain the central
line, while even counts use half-spacing offsets. Boundary connectors respect
infill anchor length and maximum anchor length. Wider floor bands are filled
once at solid row spacing, bypassing wall offsets and connectors.

Ordering is local to each fill surface. Sparse-infill anchoring queries and extrusion
use the same nominal region flow and cloud settings, so first-layer widths and adaptive
layer heights do not move points.
Layer completion reports throttled progress and checks cancellation.

## Spatial density

The [shared material-density field](infill-density.md) grades density by wall
depth or lobe/section profiles when Adaptive density is enabled, and independently reinforces
corners when smoothing is enabled. The defaults (Disabled and zero smoothing)
give a uniform cloud without preparing a distance field.

Fused modifiers change this field without splitting the point cloud. Candidate
bounds include every participating density, interior target and corner boost.
The parent retains its point distribution, extrusion and geometry settings.

## Structured points

FCC and Cubic lattices have their [111] body diagonal vertical and rotate with
sparse infill direction. Ordered dithering selects sites according to the cube
of local density. FCC uses tetrahedral parity groups and coarser checkerboards;
Cubic uses a three-dimensional Bayer order. Density calibration is approximate.

Complete FCC lattices have vertical faces and faces inclined about 35 degrees
above horizontal, below the usual 45-degree guideline. Complete Cubic faces
incline about 55 degrees, making Cubic preferable for steeper walls. Thinning can
introduce shallower faces. Neither distribution guarantees self-supporting infill.

## Point adjustment and shell bonding

Relaxation iterations defaults to zero and is limited to 50. Each iteration balances Delaunay edge
lengths against local preferred spacing. Simultaneous damped steps are limited
to 10% of nominal spacing; convex-hull sites remain fixed. The saved object-wide
point cloud costs additional memory and computation and cannot eliminate
every acute cell section. Preparation collects the point cloud through local
slabs; relaxation itself operates on the complete saved point cloud.

Snap points to walls defaults off. After relaxation, points within half the
nominal spacing on either side of a wall take a bounded step toward the raw
distance field's zero contour. Grid resolution, curves and corners limit
accuracy. Relaxation and snapping preserve reproducible queries and leave the
density field unchanged; either can alter face angles.

Trace model shell defaults off. It adds one bead along the physical shell,
including holes, and excludes modifier boundaries. Touching foam walls join
the trace through cyclic junction routing, allowing trails to turn along the
shell. Traces are clipped to the fill surface and add material.

## Experimental void refinement

Voronoi hull density threshold is an absolute material-density percentage;
zero disables it. After optional relaxation and snapping, a saved point cloud
is filtered at its control points. Interior points below the cutoff are removed;
exterior neighbours remain. Fused modifiers affect density before filtering and
inherit the parent's cutoff.

Refinement examines Delaunay edges whose dual Voronoi faces incline less than
35 degrees above horizontal and span more than 10 mm in XY. Faces are clipped
to the model bounding box; triangle-centre signed-distance samples approximate
their area inside the model. The score is that area multiplied by excess span.
For the worst few faces, trial points at the face centre, lateral offsets and
below the face are evaluated against the resulting diagram. One point is kept
only if it reduces the total score. Repair points may lie below the density
cutoff.

This greedy search stops when no tested insertion improves the score or after
256 insertions. It does not find a globally minimum point set or guarantee
vaulting, support or a closed cavity. Preparation can be expensive and holds the
complete cloud in memory. The saved result is shared across layers; ordinary
face extraction, shallow floors and native extrusion remain unchanged.

Point ingestion advances in complete generation cubes. Sites remain unchanged between batches; obsolete lower cube rows are retired immediately before the next batch is generated. The persistent Delaunay network removes those vertices before inserting new ones.

Adaptive point selection reuses the object’s prepared spatial profile. A provider owns its density evaluator; direct profile evaluation has no separate result cache. Gaussian smoothing and modifier compositing cache their prepared buffers.

The object and density evaluators share the same prepared adaptive profile arrays. The object borrows its selected array; clearing the density fields also clears those entries.
