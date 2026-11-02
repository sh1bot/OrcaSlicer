#include <vector>
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Fill/AdaptiveInfillField.hpp"
#include <mutex>
#include <memory>
#include <functional>
#include "VoronoiWallDistance.hpp"
#include "libslic3r/VoronoiCloudMethod.hpp"
#include "VoronoiLattice.hpp"
#include "VoronoiPointCloud.hpp"
#include "libslic3r/InfillDensityField.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace Slic3r::Voronoi {

std::shared_ptr<const PointCloudProvider> WallDistanceCloudCache::get(const std::function<std::vector<WallSlice>()> &load_slices,
    const PointCloudSettings &settings, const WallDistanceProgress &progress, const std::vector<DensityModifier> &modifiers)
{
    auto [density, site_spacing, method, region_id] = settings;
    const auto [mode, sigma, material_density, interior_density, gradient] = density;
    if (!std::isfinite(material_density) || material_density <= 0. || material_density > 1. ||
        !std::isfinite(interior_density) || interior_density <= 0. || interior_density > 1. ||
            !std::isfinite(sigma) || sigma < 0. || sigma > InfillDensityField::maximum_smoothing_sigma)
        throw std::invalid_argument("Invalid Voronoi distance settings");
    std::lock_guard<std::mutex> lock(m_mutex);
    const auto cached = m_providers.find({mode, sigma, material_density, region_id, method, interior_density, gradient});
    if (cached != m_providers.end() && cached->second) return cached->second;
    std::shared_ptr<const InfillDensityField> table;
    WallDistanceSamples sample;
    const bool fused = std::any_of(modifiers.begin(), modifiers.end(), [region_id = settings.region_id](const auto &m) { return m.region_id == region_id; });
    if (mode != AdaptiveInfillMode::Disabled || sigma != 0. || fused) {
        table = m_fields.get(load_slices, {mode, sigma, material_density, interior_density, gradient}, progress, modifiers, region_id);
        sample = [table](const auto &points) { return table->sample_fields(points); };
    }
    auto &provider = m_providers[{mode, sigma, material_density, region_id, method, interior_density, gradient}];
    if (!provider) {
        const bool variable = mode != AdaptiveInfillMode::Disabled || sigma != 0. || fused;
        if (method != VoronoiCloudMethod::Random)
            provider = std::make_shared<LatticePointCloud>(method, variable ? sample : WallDistanceSamples{}, table ? table->maximum_density_ratio() : 1.);
        else if (variable)
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
