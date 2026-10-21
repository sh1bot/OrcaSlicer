#include "InfillDensityField.hpp"
#include "ExPolygon.hpp"
#include "Fill/AdaptiveInfillField.hpp"
#include "PrintConfig.hpp"
#include "Point.hpp"
#include <functional>
#include <memory>
#include <mutex>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <array>
#include <utility>
#include <stdexcept>

namespace Slic3r {
namespace {
void report(const WallDistanceProgress &progress, WallDistanceStage stage, double fraction)
{
    if (progress) progress(stage, fraction);
}
AdaptiveInfillFields density_geometry(InfillDensityCache &cache,
    const std::function<std::vector<WallSlice>()> &load_slices, AdaptiveInfillMode mode,
    const WallDistanceProgress &progress)
{
    AdaptiveInfillFields shapes;
    if (int(mode) < 0 || mode >= AdaptiveInfillMode::Count)
        throw std::invalid_argument("Invalid adaptive infill mode");
    const bool radial = mode == AdaptiveInfillMode::Lobes || mode == AdaptiveInfillMode::NormalX ||
        mode == AdaptiveInfillMode::NormalY || mode == AdaptiveInfillMode::NormalZ;
    const auto selected = radial ? mode : AdaptiveInfillMode::SmoothBlend;
    shapes[size_t(selected)] = cache.geometry(load_slices, selected, false, progress);
    shapes[size_t(mode)] = shapes[size_t(selected)];
    return shapes;
}

} // namespace

struct InfillDensityField::Data {
    AdaptiveInfillFields geometry;
    AdaptiveInfillDistanceField original;
    double reference_density;
    double maximum_ratio;
    AdaptiveInfillMode mode;
    double interior_density;
    AdaptiveInfillGradient gradient;
    explicit Data(std::vector<WallSlice> slices, const WallDistanceProgress &progress,
                   AdaptiveInfillMode mode, double material_density, double interior_density, AdaptiveInfillGradient gradient, AdaptiveInfillFields geometry)
        : geometry(std::move(geometry)), reference_density(material_density), mode(mode),
          interior_density(interior_density), gradient(gradient)
    {
        if (!std::isfinite(material_density) || material_density <= 0. || material_density > 1. ||
            !std::isfinite(interior_density) || interior_density <= 0. || interior_density > 1.) throw std::invalid_argument("Invalid material density settings");
        maximum_ratio = (mode == AdaptiveInfillMode::Disabled ? material_density : std::max(material_density, interior_density)) / material_density;
        const auto selected = mode;
        if (!this->geometry[size_t(selected)]) {
            InfillDensityCache cache;
            this->geometry = density_geometry(cache, [&] { return slices; }, mode, progress);
        }
        original = this->geometry[size_t(selected)]->distances();
        report(progress, WallDistanceStage::Complete, 1.);
    }

    double local_ratio(const Vec3d &point, double raw, AdaptiveInfillMode mode,
        bool graded, double interior_ratio, AdaptiveInfillGradient gradient) const
    {
        if (!graded || raw >= 0.) return 1.;
        const auto &shape = geometry[size_t(mode)];
        if (!shape || shape->empty()) return 1.;
        const auto at_depth = [&](double depth) {
            return adaptive_density_scale(interior_ratio, gradient, depth);
        };
        if (mode == AdaptiveInfillMode::DistanceWarp || mode == AdaptiveInfillMode::SmoothBlend || mode == AdaptiveInfillMode::SteppedShells)
            return at_depth(shape->depth(point));
        AdaptiveInfillField::Radials radials;
        const size_t count = shape->radial(point, radials);
        if (count == 0) return 1.;
        double ratio = 0.;
        for (size_t j = 0; j < count; ++j) ratio += radials[j].weight * at_depth(1. - radials[j].t);
        return ratio;
    }

    double density(const Vec3d &point) const
    {
        if (original.empty()) return reference_density;
        const auto where = original.interpolation(point);
        const double raw = original.sample(where);
        return std::min(1., reference_density * local_ratio(point, raw, mode,
            mode != AdaptiveInfillMode::Disabled && interior_density != reference_density, interior_density / reference_density, gradient));
    }
};

InfillDensityField::InfillDensityField(std::vector<WallSlice> slices, const WallDistanceProgress &progress, AdaptiveInfillMode mode, double material_density,
                                     double interior_density, AdaptiveInfillGradient gradient, AdaptiveInfillFields geometry)
    : m_data(std::make_shared<Data>(std::move(slices), progress, mode, material_density, interior_density, gradient, std::move(geometry))) {}

double InfillDensityField::grid_spacing() const { return m_data->original.grid_spacing(); }
double InfillDensityField::maximum_density_ratio() const { return m_data->maximum_ratio; }

std::vector<InfillDensitySample> InfillDensityField::sample_fields(const std::vector<Vec3d> &points) const
{
    std::vector<InfillDensitySample> result;
    result.reserve(points.size());
    for (const auto &point : points)
        result.push_back({m_data->original.signed_distance(point), m_data->density(point) / m_data->reference_density});
    return result;
}

std::vector<double> InfillDensityField::sample(const std::vector<Vec3d> &points) const
{
    return m_data->original.sample(points);
}

double InfillDensityField::density(const Vec3d &point) const { return m_data->density(point); }

double InfillDensityField::maximum_density() const { return m_data->maximum_ratio * m_data->reference_density; }

std::shared_ptr<const AdaptiveInfillField> InfillDensityCache::geometry(
    const std::function<std::vector<WallSlice>()> &load_slices, AdaptiveInfillMode mode,
    bool signed_distances, const WallDistanceProgress &progress)
{
    const int axis = mode == AdaptiveInfillMode::NormalX ? 0 : mode == AdaptiveInfillMode::NormalY ? 1 : mode == AdaptiveInfillMode::NormalZ ? 2 : -1;
    std::lock_guard<std::mutex> lock(m_mutex);
    auto &field = m_geometry[{mode, signed_distances}];
    if (!field) {
        const auto slices = load_slices();
        BoundingBox bbox;
        for (const auto &slice : slices) bbox.merge(get_extents(*slice.expolygons));
        std::mutex progress_mutex;
        WallDistanceStage stage = WallDistanceStage::Distance;
        double fraction = 0.;
        const auto notify = [&](WallDistanceStage next_stage, double next_fraction) {
            std::lock_guard<std::mutex> guard(progress_mutex);
            stage = next_stage;
            fraction = next_fraction;
            report(progress, stage, fraction);
        };
        const auto cancel = [&] {
            std::lock_guard<std::mutex> guard(progress_mutex);
            report(progress, stage, fraction);
        };
        auto &distances = m_distances[{axis, signed_distances}];
        if (!distances) {
            const auto distance_progress = [&](AdaptiveInfillDistanceField::Stage stage, double fraction) {
                if (stage == AdaptiveInfillDistanceField::Stage::Distance) notify(WallDistanceStage::Distance, fraction);
                else if (stage == AdaptiveInfillDistanceField::Stage::Encoding) notify(WallDistanceStage::Encoding, fraction / 2.);
            };
            distances = std::make_shared<AdaptiveInfillDistanceField>(slices, bbox, cancel, distance_progress, signed_distances, axis);
        }
        notify(WallDistanceStage::Encoding, 0.5);
        field = std::make_shared<AdaptiveInfillField>(slices, bbox, mode, cancel, distances);
    }
    return field;
}

std::shared_ptr<const InfillDensityField> InfillDensityCache::get(
    const std::function<std::vector<WallSlice>()> &load_slices, const InfillDensitySettings &s,
    const WallDistanceProgress &progress)
{
    auto shapes = density_geometry(*this, load_slices, s.mode, progress);
    std::lock_guard<std::mutex> lock(m_mutex);
    auto &field = m_fields[{s.mode, s.material_density, s.interior_density, s.gradient}];
    if (!field)
        field = std::make_shared<InfillDensityField>(std::vector<WallSlice>{}, progress,
            s.mode, s.material_density, s.interior_density, s.gradient, std::move(shapes));
    return field;
}

void InfillDensityCache::clear()
{
    std::lock_guard<std::mutex> lock(m_mutex);
    m_fields.clear();
    m_geometry.clear();
    m_distances.clear();
}

} // namespace Slic3r
