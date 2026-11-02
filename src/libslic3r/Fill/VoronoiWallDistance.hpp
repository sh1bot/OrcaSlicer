#pragma once

#include <vector>
#include "libslic3r/Fill/AdaptiveInfillField.hpp"
#include <tuple>
#include <mutex>
#include <memory>
#include <map>
#include <functional>
#include "../InfillDensityField.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "VoronoiPointCloud.hpp"
#include "../VoronoiCloudMethod.hpp"

namespace Slic3r::Voronoi {

struct PointCloudSettings {
    InfillDensitySettings density;
    double site_spacing = 0.;
    VoronoiCloudMethod method = VoronoiCloudMethod::Random;
    int region_id = -1;
};

// Owned by a PrintObject. Initialization is serialized; providers and distance
// data are immutable and can subsequently be used by parallel layer workers.
class WallDistanceCloudCache {
public:
    std::shared_ptr<const PointCloudProvider> get(const std::function<std::vector<WallSlice>()> &load_slices,
        const PointCloudSettings &settings = {}, const WallDistanceProgress &progress = {},
        const std::vector<DensityModifier> &modifiers = {});
    InfillDensityCache &density_fields() { return m_fields; }
    void clear();
private:
    std::mutex m_mutex;
    InfillDensityCache m_fields;
    std::map<std::tuple<AdaptiveInfillMode, double, double, int, VoronoiCloudMethod, double, AdaptiveInfillGradient>, std::shared_ptr<const PointCloudProvider>> m_providers;
};

} // namespace Slic3r::Voronoi
