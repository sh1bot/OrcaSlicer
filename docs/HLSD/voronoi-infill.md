# Three-dimensional Voronoi sparse infill

The `voronoi` pattern prints layer sections of faces between nearest 3D control
points. Point generation, cell geometry, wall routing and extrusion are separate
stages, so other point distributions can reuse the same infill implementation.

## Point cloud and geometry

`PointCloudProvider` supplies reproducible XYZ points through half-open box
queries and defines the density-to-spacing mapping. The default provider uses
spatially seeded Poisson sampling. Other providers may be structured, spatially
varying or finite; queries must reproduce the same points regardless of order.
An unbounded provider must eventually supply points.

Sparse infill direction rotates the random cloud about Z through the object's
bounding-box centre; model alignment adds the object's XY rotation. Density is
sampled at rotated positions while candidate identities remain unchanged. The
angle stays fixed across layers; per-layer rotation templates do not rotate the
3D constellation.

`ActivePointCloud` expands its known region until nearest-site covering bounds
prove that unknown points cannot affect the working slab. Points outside the
model participate; resulting walls and floors are clipped to the fill surface.
Points sufficiently far below the slab are retired. Each infill engine owns its
active constellation and a persistent CGAL 3D Delaunay triangulation. Point-set
changes insert new vertices and batch-remove retired clusters, retaining unaffected
connectivity. Vertex handles are indexed in a sorted vector; sampled XYZ values
remain in a contiguous vector. CGAL stores vertices and tetrahedra in chunked
containers with reusable slots. Copies of a filler own independent networks.
During normal slicing, each worker processes an ascending batch of layers and
caches its engines by point provider; temporary surface fillers reuse those
engines. The cache is released at the end of the batch and is never shared between workers.

CGAL's 2D power diagram, weighted by vertical distance, supplies wall sections
and covering bounds. Its weights change with Z, so this planar diagram is rebuilt
for each section. Points are retired before updating the persistent 3D network,
which supplies shallow-face floor polygons. Geometry clipping removes overlapping
floors and wall sections within floors while retaining their supporting boundaries.

## Extrusion and routing

Faces become single wall lines when successive layers have at least half-bead
overlap. Shallower faces become floor bands, rounded outward to Orca's flow row
spacing while keeping the supported boundary fixed. At constant layer height,
continuing bands have at least 50% overlap before clipping; face ends, holes and
changing layer heights limit that guarantee.

Native Rectilinear fill prints floor rows along the face section in alternating
zigzags, advancing from the supported side. Orca's bead width, flow and spacing
are retained. Floor row ends extend into nearby previous-layer perimeter bead
footprints when reachable; sparse slice interiors are not treated as solid support.
This does not select bridge flow or guarantee support at both ends.

Fill Multiline uses Orca's native offset outlines around routed wall paths, while
floors are filled once at native solid row spacing. Nominal point spacing grows
with bead count to approximately retain the selected density. Native boundary
connectors respect infill anchor length and maximum anchor length; a zero maximum
disables them.

Wall graph vertices snap to fixed representatives. Branch pairing favours
straight crossings; three-way junctions cycle turns across layers using their
generating points. Trails traverse each wall edge once. Wall trails and complete
floor groups are ordered together by nearest permitted entry, using each group's
actual exit. Walls may reverse; floors retain their supported-side-first sweep.
Ordering is local to a fill surface, and disconnected paths still require travel.

Density is statistical: spacing calibration accounts for floor material, but
clipping and row rounding affect the result. Nominal extrusion dimensions keep
first-layer widths and adaptive layer heights from changing the point field.
Infill reports completed layers with throttled UI updates and cancellation.

## Wall proximity and smoothing

Wall-proximity decay defaults to 6 mm; zero gives a uniform base. Gaussian sigma
defaults to 8 mm; zero disables smoothing. With both zero, the original Poisson
provider is used without a distance table.

Along flat walls, local density stays at the selected setting within 2 mm,
then decays exponentially to a 1% tail after five decay lengths. The magnitude
of the difference between original and smoothed signed distances reduces the
effective depth near both convex and concave corners. This can boost density to
approximately 1.39 times the setting. Zero decay retains the corner boost with
a 6 mm response length. A reproducibly thinned Poisson candidate field follows
this density; thinning does not constrain face slopes or ensure vaulted roofs.

The distance field uses Orca's completed slice outlines, including holes,
floors and ceilings. `EdgeGrid` supplies the approximate planar signed distance
for each layer slab. Two Z sweeps accumulate distance to material and air using
lower envelopes of weighted squared endpoint distances. Identical consecutive
profiles merge; thin slabs between output planes and empty gaps still contribute.
Accumulation is linear in the number of profiles and output heights per XY
column, rather than scanning the whole volume for every profile.

The immutable cubic tables normally have 1 mm spacing, logarithmic distance
bytes covering 0–1,000 mm, and sign bits. A 32-million-node limit coarsens large
objects. Sampling interpolates decoded signed distances; density saturation
happens afterwards. The field is approximate, with no geometry fallback.

Gaussian smoothing uses separable CPU OpenCV filters with fixed-point working
values. Temporary filter borders extrapolate outward slopes; zero sigma reuses
the original field. Tables are shared by smoothing width, and point providers
retain their own decay setting. Reslicing clears the per-object cache.
Construction and smoothing report phase progress and check cancellation.
