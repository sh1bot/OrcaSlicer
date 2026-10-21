#pragma once

#include "Point.hpp"
#include "ExPolygon.hpp"
#include "Fill/AdaptiveInfillField.hpp"
#include "BoundingBox.hpp"
#include "PrintConfig.hpp"
#include <algorithm>
#include <array>
#include <cstddef>
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

enum class WallDistanceStage { Distance, Encoding, SmoothXY, SmoothZ, Complete };
// Called synchronously during construction, including within long phases.
// The callback may throw to cancel construction; no partial table is cached.
using WallDistanceProgress = std::function<void(WallDistanceStage, double)>;

// Adaptive material density evaluated from shared object geometry.
class InfillDensityField {
public:
    static constexpr double maximum_smoothing_sigma = 64.;
    explicit InfillDensityField(std::vector<WallSlice> slices, double smoothing_sigma = 0.,
                               const WallDistanceProgress &progress = {}, AdaptiveInfillMode mode = AdaptiveInfillMode::SmoothBlend, double material_density = 0.25,
                               double interior_density = 0.05, AdaptiveInfillGradient gradient = AdaptiveInfillGradient::Linear, AdaptiveInfillFields geometry = {});
    // Absolute material fraction, from zero to one.
    double density(const Vec3d &point) const;
    double maximum_density() const;
    static constexpr double near_wall_distance = 2.;
    // Relative material density; zero maximum_depth selects a uniform base.
    static double density_ratio(double distance, double sampled_mean, double maximum_depth = 0., double interior_ratio = 0.01,
                                AdaptiveInfillGradient gradient = AdaptiveInfillGradient::Exponential);
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
    double sigma = 0.;
    double material_density = 0.25;
    double interior_density = 0.05;
    AdaptiveInfillGradient gradient = AdaptiveInfillGradient::Linear;

    static InfillDensitySettings from_config(const PrintRegionConfig &config)
    {
        return {config.adaptive_infill_mode.value, config.infill_smoothing_sigma.value,
            0.01 * config.sparse_infill_density.value, 0.01 * std::max(1., config.adaptive_infill_interior_density.value),
            config.adaptive_infill_gradient.value};
    }
};

class InfillDensityCache {
public:
    std::shared_ptr<const InfillDensityField> get(const std::function<std::vector<WallSlice>()> &load_slices,
        const InfillDensitySettings &settings, const WallDistanceProgress &progress = {});
    std::shared_ptr<const AdaptiveInfillField> geometry(const std::function<std::vector<WallSlice>()> &load_slices,
        AdaptiveInfillMode mode, bool signed_distances = true,
        const WallDistanceProgress &progress = {});
    const AdaptiveInfillFields &profiles(bool signed_distances) const { return m_geometry[size_t(signed_distances)]; }
    void clear();
private:
    std::mutex m_mutex;
    // The object and density evaluators use these same prepared profile arrays.
    std::array<AdaptiveInfillFields, 2> m_geometry;
    std::map<std::pair<int, bool>, std::shared_ptr<const AdaptiveInfillDistanceField>> m_distances;
    std::map<std::tuple<AdaptiveInfillMode, double, double, double, AdaptiveInfillGradient>, std::shared_ptr<const InfillDensityField>> m_fields;
};

} // namespace Slic3r
