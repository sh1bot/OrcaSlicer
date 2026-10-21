# Shared adaptive infill density

Voronoi uses Orca's common adaptive modes, interior density and gradient.
Disabled gives a uniform cloud. Distance warp, Smooth blend and Stepped shells
all select the same continuous wall-depth point density. Lobes and Normal X/Y/Z
blend density contributions from neighbouring lobes or sections. These profiles
change point selection without warping the cloud coordinates.

Surface density is the sparse infill setting; interior density is an absolute
percentage. Linear, Quadratic and Exponential use the existing adaptive gradient.
Distance modes normalise depth to the deepest point of each connected body;
Lobes and Normal modes use lobe and section centers.

Density queries evaluate the selected upstream profile and gradient directly.
There is no additional material-density grid. Adaptive geometry and point providers
are cached by settings and shared across layer workers. Construction reports
progress and supports cancellation.
The voxel grid starts at 0.5 mm spacing and coarsens to stay within about one
million nodes. Narrow features smaller than a voxel may disappear.
