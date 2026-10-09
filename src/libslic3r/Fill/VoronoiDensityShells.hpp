#pragma once

#include "VoronoiPointCloud.hpp"
#include <functional>
#include "../BoundingBox.hpp"
#include <optional>
#include <vector>

namespace Slic3r::Voronoi {

// Prepared once per object/density. Surface meshes are temporary; layer workers
// query the immutable, height-sorted sites using the normal provider interface.
class DensityShellPointCloud : public WallDistancePointCloud {
public:
    // angle is the object-wide XY lattice orientation, in radians.
    DensityShellPointCloud(WallDistanceSamples sample, const BoundingBoxf3 &bounds,
                          double site_spacing, double decay,
                          const std::function<void(double)> &progress = {}, double angle = 0.);
    std::vector<Vec3d> points_in(const BoundingBoxf3 &box, double spacing) const override;
    std::optional<BoundingBoxf3> extent() const override { return m_bounds; }
private:
    BoundingBoxf3 m_bounds;
    std::vector<Vec3d> m_sites;
};

} // namespace Slic3r::Voronoi
