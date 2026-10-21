#include <vector>
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Fill/AdaptiveInfillField.hpp"
#include <mutex>
#include <memory>
#include <functional>
#include "VoronoiWallDistance.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/VoronoiCloudMethod.hpp"
#include "VoronoiLattice.hpp"
#include "VoronoiWallSnap.hpp"
#include "VoronoiRelaxation.hpp"
#include "VoronoiVoid.hpp"
#include "VoronoiPointCloud.hpp"
#include "libslic3r/InfillDensityField.hpp"
#include "libslic3r/libslic3r.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace Slic3r::Voronoi {

std::shared_ptr<const PointCloudProvider> WallDistanceCloudCache::get(const std::function<std::vector<WallSlice>()> &load_slices,
    const PointCloudSettings &settings, const WallDistanceProgress &progress, const std::vector<DensityModifier> &modifiers)
{
    auto [density, site_spacing, method, relaxation_iterations, angle, snap, region_id, hull_threshold] = settings;
    const auto [mode, sigma, material_density, interior_density, gradient] = density;
    if (!std::isfinite(site_spacing) || site_spacing < 0. || relaxation_iterations < 0 || relaxation_iterations > RelaxedPointCloud::maximum_iterations ||
        ((snap || relaxation_iterations > 0 || hull_threshold > 0.) && site_spacing == 0.))
        throw std::invalid_argument("Invalid Voronoi point preparation settings");
    if (!std::isfinite(material_density) || material_density <= 0. || material_density > 1. ||
        !std::isfinite(interior_density) || interior_density <= 0. || interior_density > 1. ||
            !std::isfinite(sigma) || sigma < 0. || sigma > InfillDensityField::maximum_smoothing_sigma)
        throw std::invalid_argument("Invalid Voronoi distance settings");
    if (!std::isfinite(hull_threshold) || hull_threshold < 0. || hull_threshold > 1.)
        throw std::invalid_argument("Invalid Voronoi hull threshold");
    if (!std::isfinite(angle)) throw std::invalid_argument("Invalid Voronoi cloud angle");
    angle = std::fmod(angle, 2. * std::acos(-1.));
    if (angle < 0.) angle += 2. * std::acos(-1.);
    std::lock_guard<std::mutex> lock(m_mutex);
    const auto cached = m_providers.find({mode, sigma, material_density, region_id, hull_threshold > 0. || snap || relaxation_iterations > 0 ? site_spacing : 0., method, relaxation_iterations, angle, snap, interior_density, gradient, hull_threshold});
    if (cached != m_providers.end() && cached->second) return cached->second;
    const auto bounds = [&] {
        BoundingBoxf3 box;
        for (const auto &slice : load_slices())
            for (const auto &polygon : *slice.expolygons) {
                const auto xy = polygon.contour.bounding_box();
                if (!xy.defined) continue;
                box.merge(Vec3d(unscale<double>(xy.min.x()), unscale<double>(xy.min.y()), slice.bottom_z));
                box.merge(Vec3d(unscale<double>(xy.max.x()), unscale<double>(xy.max.y()), slice.top_z));
            }
        return box;
    };
    std::shared_ptr<const InfillDensityField> table;
    WallDistanceSamples sample;
    const bool fused = std::any_of(modifiers.begin(), modifiers.end(), [region_id = settings.region_id](const auto &m) { return m.region_id == region_id; });
    if (mode != AdaptiveInfillMode::Disabled || sigma != 0. || snap || fused || hull_threshold > 0.) {
        table = m_fields.get(load_slices, {mode, sigma, material_density, interior_density, gradient}, progress, modifiers, region_id);
        sample = [table](const auto &points) { return table->sample_fields(points); };
    }
    auto &provider = m_providers[{mode, sigma, material_density, region_id, 0., method, 0, angle, false, interior_density, gradient, 0.}];
    if (!provider) {
        Vec2d center = Vec2d::Zero();
        if (angle != 0.) {
            const BoundingBoxf3 box = bounds();
            if (box.defined) center = box.center().head<2>();
        }
        const bool variable = mode != AdaptiveInfillMode::Disabled || sigma != 0. || fused;
        if (method != VoronoiCloudMethod::Random)
            provider = std::make_shared<LatticePointCloud>(method, variable ? sample : WallDistanceSamples{}, table ? table->maximum_density_ratio() : 1., angle, center);
        else if (variable)
            provider = std::make_shared<WallDistancePointCloud>(sample, table->maximum_density_ratio(), angle, center);
        else
            provider = std::make_shared<PoissonPointCloud>(angle, center);
    }
    std::shared_ptr<const PointCloudProvider> result = provider;
    if (relaxation_iterations > 0) {
        auto &relaxed = m_providers[{mode, sigma, material_density, region_id, site_spacing, method, relaxation_iterations, angle, false, interior_density, gradient, 0.}];
        if (!relaxed)
            relaxed = std::make_shared<RelaxedPointCloud>(provider, bounds(), site_spacing, relaxation_iterations,
                [&](double fraction) { if (progress) progress(WallDistanceStage::Relaxation, fraction); });
        result = relaxed;
    }
    if (snap) {
        auto &snapped = m_providers[{mode, sigma, material_density, region_id, site_spacing, method, relaxation_iterations, angle, true, interior_density, gradient, 0.}];
        if (!snapped)
            snapped = std::make_shared<SnappedPointCloud>(result, sample, site_spacing, table->grid_spacing());
        result = snapped;
    }
    if (hull_threshold > 0.) {
        auto &refined = m_providers[{mode, sigma, material_density, region_id, site_spacing, method,
            relaxation_iterations, angle, snap, interior_density, gradient, hull_threshold}];
        if (!refined)
            refined = std::make_shared<VoidPointCloud>(result, bounds(), site_spacing, sample,
                hull_threshold / material_density, VoidRefinement{},
                [&](double fraction) { if (progress) progress(WallDistanceStage::VoidRefinement, fraction); });
        result = refined;
    }
    return result;
}

void WallDistanceCloudCache::clear()
{
    std::lock_guard<std::mutex> lock(m_mutex);
    m_providers.clear();
    m_fields.clear();
}

} // namespace Slic3r::Voronoi
