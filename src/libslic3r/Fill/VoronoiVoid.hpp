#pragma once

#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Point.hpp"
#include "VoronoiPointCloud.hpp"
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <vector>

namespace Slic3r::Voronoi {

struct VoidRefinement {
    double minimum_angle = 35.;
    double maximum_span = 10.;
    size_t maximum_insertions = 256;
};

struct VoidRefinementResult {
    size_t removed = 0, inserted = 0, remaining_faces = 0;
    double initial_cost = 0., final_cost = 0.;
};

// Materialise, filter, then refine a saved constellation before layer extraction.
class VoidPointCloud : public PointCloudProvider {
public:
    VoidPointCloud(std::shared_ptr<const PointCloudProvider> source, const BoundingBoxf3 &bounds,
                   double spacing, WallDistanceSamples samples, double relative_threshold,
                   const VoidRefinement &settings = {}, const std::function<void(double)> &progress = {});
    double spacing(double width, double height, double density) const override
        { return m_source->spacing(width, height, density); }
    std::vector<Vec3d> points_in(const BoundingBoxf3 &box, double spacing) const override;
    std::optional<BoundingBoxf3> extent() const override { return m_bounds; }
    std::vector<double> local_spacing(const std::vector<Vec3d> &points, double spacing) const override
        { return m_source->local_spacing(points, spacing); }
    const VoidRefinementResult &result() const { return m_result; }
private:
    std::shared_ptr<const PointCloudProvider> m_source;
    BoundingBoxf3 m_bounds;
    std::vector<Vec3d> m_sites;
    VoidRefinementResult m_result;
};

} // namespace Slic3r::Voronoi
