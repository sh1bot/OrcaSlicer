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

struct DensityModifierSlice {
    double bottom_z, top_z;
    ExPolygons contours;
};

// Stored in Orca's modifier order, scoped to the unsplit parent print region.
struct DensityModifier {
    int region_id;
    double material_density;
    double blend_width;
    std::vector<DensityModifierSlice> slices;
    std::optional<double> smoothing_sigma, interior_density;
    std::optional<AdaptiveInfillGradient> gradient;
    std::optional<AdaptiveInfillMode> mode;
};

enum class WallDistanceStage { Distance, Encoding, SmoothXY, SmoothZ, MaterialDensity, Modifiers, Relaxation, VoidRefinement, Complete };
// Called synchronously during construction, including within long phases.
// The callback may throw to cancel construction; no partial table is cached.
using WallDistanceProgress = std::function<void(WallDistanceStage, double)>;

// Adaptive material density evaluated from shared object geometry.
class InfillDensityField {
public:
    static constexpr double maximum_smoothing_sigma = 64.;
    explicit InfillDensityField(std::vector<WallSlice> slices, double smoothing_sigma = 0.,
                               const WallDistanceProgress &progress = {}, AdaptiveInfillMode mode = AdaptiveInfillMode::SmoothBlend, double material_density = 0.25,
                               const std::vector<DensityModifier> &modifiers = {}, int region_id = -1, double interior_density = 0.05, AdaptiveInfillGradient gradient = AdaptiveInfillGradient::Linear, AdaptiveInfillFields geometry = {});
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
    bool slice_slabs = false;

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
        const InfillDensitySettings &settings, const WallDistanceProgress &progress = {},
        const std::vector<DensityModifier> &modifiers = {}, int region_id = -1);
    std::shared_ptr<const AdaptiveInfillField> geometry(const std::function<std::vector<WallSlice>()> &load_slices,
        AdaptiveInfillMode mode, bool signed_distances = true,
        const WallDistanceProgress &progress = {}, bool slice_slabs = false);
    const AdaptiveInfillFields &profiles(bool signed_distances, bool slice_slabs = false) const
        { return m_geometry[size_t(signed_distances) + 2 * size_t(slice_slabs)]; }
    void clear();
private:
    std::mutex m_mutex;
    // The object and density evaluators use these same prepared profile arrays.
    std::array<AdaptiveInfillFields, 4> m_geometry;
    std::map<std::tuple<int, bool, bool>, std::shared_ptr<const AdaptiveInfillDistanceField>> m_distances;
    std::map<std::tuple<AdaptiveInfillMode, double, double, int, double, AdaptiveInfillGradient, bool>, std::shared_ptr<const InfillDensityField>> m_fields;
};

} // namespace Slic3r
