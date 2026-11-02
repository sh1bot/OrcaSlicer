#include <vector>
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Fill/AdaptiveInfillField.hpp"
#include <mutex>
#include <memory>
#include <functional>
#include "VoronoiWallDistance.hpp"
#include "VoronoiPointCloud.hpp"
#include "libslic3r/InfillDensityField.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace Slic3r::Voronoi {

std::shared_ptr<const PointCloudProvider> WallDistanceCloudCache::get(const std::function<std::vector<WallSlice>()> &load_slices,
    const PointCloudSettings &settings, const WallDistanceProgress &progress)
{
    auto [density, site_spacing] = settings;
    const auto [mode, material_density, interior_density, gradient] = density;
    if (!std::isfinite(material_density) || material_density <= 0. || material_density > 1. ||
        !std::isfinite(interior_density) || interior_density <= 0. || interior_density > 1.)
        throw std::invalid_argument("Invalid Voronoi distance settings");
    std::lock_guard<std::mutex> lock(m_mutex);
    const auto cached = m_providers.find({mode, material_density, interior_density, gradient});
    if (cached != m_providers.end() && cached->second) return cached->second;
    std::shared_ptr<const InfillDensityField> table;
    WallDistanceSamples sample;
    if (mode != AdaptiveInfillMode::Disabled) {
        table = m_fields.get(load_slices, {mode, material_density, interior_density, gradient}, progress);
        sample = [table](const auto &points) { return table->sample_fields(points); };
    }
    auto &provider = m_providers[{mode, material_density, interior_density, gradient}];
    if (!provider) {
        const bool variable = mode != AdaptiveInfillMode::Disabled;
        if (variable)
            provider = std::make_shared<WallDistancePointCloud>(sample, table->maximum_density_ratio());
        else
            provider = std::make_shared<PoissonPointCloud>();
    }
    return provider;
}

void WallDistanceCloudCache::clear()
{
    std::lock_guard<std::mutex> lock(m_mutex);
    m_providers.clear();
    m_fields.clear();
}

} // namespace Slic3r::Voronoi
