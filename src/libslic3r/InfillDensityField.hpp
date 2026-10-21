#pragma once

#include "Point.hpp"
#include "ExPolygon.hpp"
#include "Fill/AdaptiveInfillField.hpp"
#include "BoundingBox.hpp"
#include "PrintConfig.hpp"
#include <algorithm>
#include <functional>
#include <memory>
#include <map>
#include <mutex>
#include <optional>
#include <utility>
#include <tuple>
#include <vector>

namespace Slic3r {

// Raw signed distance in millimeters (negative inside) and prepared material
// density relative to the parent setting.
struct InfillDensitySample {
    double distance;
    double density;
};

using WallSlice = AdaptiveInfillDistanceField::Slice;

enum class WallDistanceStage { Distance, Encoding, Complete };
// Called synchronously during construction, including within long phases.
// The callback may throw to cancel construction; no partial table is cached.
using WallDistanceProgress = std::function<void(WallDistanceStage, double)>;

// Adaptive material density evaluated from shared object geometry.
class InfillDensityField {
public:
    explicit InfillDensityField(std::vector<WallSlice> slices, const WallDistanceProgress &progress = {}, AdaptiveInfillMode mode = AdaptiveInfillMode::SmoothBlend, double material_density = 0.25,
                               double interior_density = 0.05, AdaptiveInfillGradient gradient = AdaptiveInfillGradient::Linear, AdaptiveInfillFields geometry = {});
    // Absolute material fraction, from zero to one.
    double density(const Vec3d &point) const;
    double maximum_density() const;
    double grid_spacing() const;
    double maximum_density_ratio() const;
    std::vector<double> sample(const std::vector<Vec3d> &points) const;
    // Signed distances and densities relative to the parent setting for point sampling.
    std::vector<InfillDensitySample> sample_fields(const std::vector<Vec3d> &points) const;
private:
    friend class InfillDensityCache;
    struct Data;
    std::shared_ptr<const Data> m_data;
};

struct InfillDensitySettings {
    AdaptiveInfillMode mode = AdaptiveInfillMode::Disabled;
    double material_density = 0.25;
    double interior_density = 0.05;
    AdaptiveInfillGradient gradient = AdaptiveInfillGradient::Linear;

    static InfillDensitySettings from_config(const PrintRegionConfig &config)
    {
        return {config.adaptive_infill_mode.value,
            0.01 * config.sparse_infill_density.value, 0.01 * std::max(1., config.adaptive_infill_interior_density.value),
            config.adaptive_infill_gradient.value};
    }
};

class InfillDensityCache {
public:
    std::shared_ptr<const InfillDensityField> get(const std::function<std::vector<WallSlice>()> &load_slices,
        const InfillDensitySettings &settings, const WallDistanceProgress &progress = {});
    std::shared_ptr<const AdaptiveInfillField> geometry(const std::function<std::vector<WallSlice>()> &load_slices,
        AdaptiveInfillMode mode, bool signed_distances = false,
        const WallDistanceProgress &progress = {});
    void clear();
private:
    std::mutex m_mutex;
    std::map<std::pair<AdaptiveInfillMode, bool>, std::shared_ptr<const AdaptiveInfillField>> m_geometry;
    std::map<std::pair<int, bool>, std::shared_ptr<const AdaptiveInfillDistanceField>> m_distances;
    std::map<std::tuple<AdaptiveInfillMode, double, double, AdaptiveInfillGradient>, std::shared_ptr<const InfillDensityField>> m_fields;
};

} // namespace Slic3r
