#pragma once

#include "libslic3r/InfillDensityField.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "VoronoiPointCloud.hpp"
#include "../VoronoiCloudMethod.hpp"
#include <Eigen/Core>
#include <vector>

namespace Slic3r::Voronoi {

// Cubic or FCC sites with a vertical [111] axis. Ordered thinning follows density.
class LatticePointCloud : public WallDistancePointCloud {
public:
    explicit LatticePointCloud(VoronoiCloudMethod method, WallDistanceSamples distance = {}, double maximum_ratio = InfillDensityField::density_ratio(0., 2.),
                               double angle = 0., const Vec2d &center = Vec2d::Zero());
    std::vector<Vec3d> points_in(const BoundingBoxf3 &box, double spacing) const override;
    BoundingBoxf3 generation_bounds(const BoundingBoxf3 &box) const override { return box; }
    std::vector<Vec3d> points_in_cells(const Vec3i64 &first, const Vec3i64 &last, double spacing) const override
        { return PointCloudProvider::points_in_cells(first, last, spacing); }

private:
    VoronoiCloudMethod m_method;
    bool m_variable;
    Eigen::Matrix3d m_basis, m_inverse;
};

} // namespace Slic3r::Voronoi
