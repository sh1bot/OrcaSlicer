#include "InfillDensityField.hpp"
#include "ExPolygon.hpp"
#include "Fill/AdaptiveInfillField.hpp"
#include <optional>
#include "InfillMaterialDensity.hpp"
#include "opencv2/core/mat.hpp"
#include "opencv2/core/base.hpp"
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
#include <limits>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

namespace Slic3r {
namespace {
constexpr double minimum_density_ratio = 0.01;
constexpr double wall_falloff = 6.;
void report(const WallDistanceProgress &progress, WallDistanceStage stage, double fraction)
{
    if (progress) progress(stage, fraction);
}
// OpenCV filters physical millimetres directly, without quantising distances.
constexpr float maximum_distance = float(AdaptiveInfillDistanceField::far_distance);

cv::Mat extrapolate_border(const cv::Mat &source, int radius_x, int radius_y, float maximum_step)
{
    cv::Mat padded;
    cv::copyMakeBorder(source, padded, radius_y, radius_y, radius_x, radius_x,
                       cv::BORDER_REPLICATE | cv::BORDER_ISOLATED);
    // Linear outward extension preserves planar fields. Bound the slope by
    // the distance field's Lipschitz limit, and saturate at the table cutoff.
    for (int y = radius_y; radius_x > 0 && y < radius_y + source.rows; ++y) {
        float *row = padded.ptr<float>(y);
        const float left = row[radius_x], right = row[radius_x + source.cols - 1];
        const float left_step = std::clamp(left - row[radius_x + 1], 0.f, maximum_step);
        const float right_step = std::clamp(right - row[radius_x + source.cols - 2], 0.f, maximum_step);
        for (int k = 1; k <= radius_x; ++k) {
            row[radius_x - k] = float(std::min(maximum_distance, left + k * left_step));
            row[radius_x + source.cols - 1 + k] = float(std::min(maximum_distance, right + k * right_step));
        }
    }
    const float *top = padded.ptr<float>(radius_y), *next = padded.ptr<float>(radius_y + 1);
    const float *bottom = padded.ptr<float>(radius_y + source.rows - 1);
    const float *previous = padded.ptr<float>(radius_y + source.rows - 2);
    for (int k = 1; k <= radius_y; ++k) {
        float *above = padded.ptr<float>(radius_y - k);
        float *below = padded.ptr<float>(radius_y + source.rows - 1 + k);
        for (int x = 0; x < padded.cols; ++x) {
            above[x] = float(std::min(maximum_distance, top[x] + k * std::clamp(top[x] - next[x], 0.f, maximum_step)));
            below[x] = float(std::min(maximum_distance, bottom[x] + k * std::clamp(bottom[x] - previous[x], 0.f, maximum_step)));
        }
    }
    return padded;
}

void smooth_distances(std::vector<float> &values, const std::array<size_t, 3> &size, double spacing, double sigma,
                      const WallDistanceProgress &progress)
{
    if (sigma == 0.) return;
    report(progress, WallDistanceStage::SmoothXY, 0.);
    const double radius_voxels = std::ceil(3. * sigma / spacing);
    if (radius_voxels > (std::numeric_limits<int>::max() - 1) / 2)
        throw std::invalid_argument("Infill smoothing radius is too large");
    const int radius = int(radius_voxels);
    const cv::Mat kernel = cv::getGaussianKernel(2 * radius + 1, sigma / spacing, CV_32F);
    const float maximum_step = float(std::min(double(maximum_distance), spacing));
    const size_t plane_size = size[0] * size[1];
    // X and Y passes operate on one plane at a time. Only the filter's
    // temporary border is extrapolated; no additional geometry is rasterized.
    for (size_t z = 0; z < size[2]; ++z) {
        cv::Mat plane(int(size[1]), int(size[0]), CV_32FC1, values.data() + z * plane_size);
        cv::Mat padded = extrapolate_border(plane, radius, radius, maximum_step);
        cv::sepFilter2D(padded, padded, CV_32F, kernel, kernel);
        padded(cv::Rect(radius, radius, plane.cols, plane.rows)).copyTo(plane);
        report(progress, WallDistanceStage::SmoothXY, double(z + 1) / size[2]);
    }
    report(progress, WallDistanceStage::SmoothZ, 0.);
    // Treat XY positions as columns for the Z pass. Limit scratch storage by
    // filtering blocks of columns rather than padding an entire second volume.
    const size_t block_width = std::max(size_t(1), size_t(256 * 1024) / (size[2] + 2 * size_t(radius)));
    const cv::Mat identity = cv::Mat::ones(1, 1, CV_32F);
    for (size_t first = 0; first < plane_size; first += block_width) {
        cv::Mat columns(int(size[2]), int(std::min(block_width, plane_size - first)), CV_32FC1,
                        values.data() + first, plane_size * sizeof(float));
        cv::Mat padded = extrapolate_border(columns, 0, radius, maximum_step);
        cv::sepFilter2D(padded, padded, CV_32F, identity, kernel);
        padded(cv::Rect(0, radius, columns.cols, columns.rows)).copyTo(columns);
        report(progress, WallDistanceStage::SmoothZ, double(std::min(first + block_width, plane_size)) / plane_size);
    }
}

std::vector<float> smoothed_distances(const AdaptiveInfillDistanceField &original, double sigma,
    const WallDistanceProgress &progress)
{
    if (sigma == 0. || original.empty()) return {};
    std::vector<float> values(original.node_count());
    for (size_t i = 0; i < values.size(); ++i) {
        if (i % 65536 == 0) report(progress, WallDistanceStage::Encoding, (1. + double(i) / values.size()) / 2.);
        values[i] = float(original.node_value(i));
    }
    smooth_distances(values, original.grid_size(), original.grid_spacing(), sigma, progress);
    return values;
}

AdaptiveInfillFields density_geometry(InfillDensityCache &cache,
    const std::function<std::vector<WallSlice>()> &load_slices, AdaptiveInfillMode mode,
    const std::vector<DensityModifier> &modifiers, int region_id,
    const WallDistanceProgress &progress)
{
    AdaptiveInfillFields shapes;
    shapes[size_t(AdaptiveInfillMode::SmoothBlend)] = cache.geometry(load_slices, AdaptiveInfillMode::SmoothBlend, true, progress);
    const auto prepare = [&](AdaptiveInfillMode mode) {
        if (int(mode) < 0 || mode >= AdaptiveInfillMode::Count)
            throw std::invalid_argument("Invalid adaptive infill mode");
        if (mode == AdaptiveInfillMode::Disabled) return;
        const bool distance = mode == AdaptiveInfillMode::DistanceWarp || mode == AdaptiveInfillMode::SmoothBlend || mode == AdaptiveInfillMode::SteppedShells;
        shapes[size_t(mode)] = distance ? shapes[size_t(AdaptiveInfillMode::SmoothBlend)] : cache.geometry(load_slices, mode, true, progress);
    };
    prepare(mode);
    for (const auto &modifier : modifiers)
        if (modifier.region_id == region_id) prepare(modifier.mode.value_or(mode));
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
    std::vector<float> smoothed;
    std::vector<float> material;
    explicit Data(std::vector<WallSlice> slices, double smoothing_sigma, const WallDistanceProgress &progress,
                   AdaptiveInfillMode mode, double material_density, const std::vector<DensityModifier> &modifiers, int region_id, double interior_density, AdaptiveInfillGradient gradient, AdaptiveInfillFields geometry)
        : geometry(std::move(geometry)), reference_density(material_density), mode(mode),
          interior_density(interior_density), gradient(gradient)
    {
        if (!std::isfinite(material_density) || material_density <= 0. || material_density > 1. ||
            !std::isfinite(interior_density) || interior_density <= 0. || interior_density > 1.) throw std::invalid_argument("Invalid material density settings");
        if (!std::isfinite(smoothing_sigma) || smoothing_sigma < 0. || smoothing_sigma > maximum_smoothing_sigma)
            throw std::invalid_argument("Invalid infill smoothing sigma");
        const bool composited = std::any_of(modifiers.begin(), modifiers.end(), [region_id](const auto &m) { return m.region_id == region_id; });
        const auto density_bound = [&](double density, double sigma, bool graded, double interior) {
            const double boost = InfillDensityField::density_ratio(0., sigma == 0. ? 0. : 2.);
            const double bound = std::min(1., std::max(density, graded ? interior : 0.) + density * (boost - 1.));
            return (composited ? std::max(bound, double(float(bound))) : bound) / reference_density;
        };
        maximum_ratio = std::max(0.01, density_bound(material_density, smoothing_sigma, mode != AdaptiveInfillMode::Disabled, interior_density));
        const auto selected = AdaptiveInfillMode::SmoothBlend;
        if (!this->geometry[size_t(selected)]) {
            InfillDensityCache cache;
            this->geometry = density_geometry(cache, [&] { return slices; }, mode, modifiers, region_id, progress);
        }
        original = this->geometry[size_t(selected)]->distances();
        smoothed = smoothed_distances(original, smoothing_sigma, progress);
        const auto &size = original.grid_size();
        const auto &origin = original.grid_origin();
        const double spacing = original.grid_spacing();
        const size_t count = original.node_count();
        if (count && composited) {
            material.resize(count);
            for (size_t i = 0; i < count; ++i) {
                if (i % 65536 == 0) report(progress, WallDistanceStage::MaterialDensity, double(i) / count);
                const Vec3d point = origin + spacing * Vec3d(double(i % size[0]),
                    double((i / size[0]) % size[1]), double(i / (size[0] * size[1])));
                const double raw = original.node_value(i);
                const double mean = smoothed.empty() ? raw : double(smoothed[i]);
                material[i] = float(std::min(1., material_density * local_ratio(point, raw, mean, mode,
                    mode != AdaptiveInfillMode::Disabled && interior_density != material_density, interior_density / material_density, gradient)));
            }
        }
        size_t modifier_index = 0;
        for (const DensityModifier &modifier : modifiers) {
            if (modifier.region_id != region_id) continue;
            if (!std::isfinite(modifier.material_density) || modifier.material_density < 0. || modifier.material_density > 1. ||
                !std::isfinite(modifier.blend_width) || modifier.blend_width < 0.)
                throw std::invalid_argument("Invalid infill density modifier");
            const double target_sigma = modifier.smoothing_sigma.value_or(smoothing_sigma);
            const double target_interior = modifier.interior_density.value_or(interior_density);
            const auto target_mode = modifier.mode.value_or(mode);
            const bool target_adaptive = target_mode != AdaptiveInfillMode::Disabled && target_interior != modifier.material_density;
            const auto target_gradient = modifier.gradient.value_or(gradient);
            if (!std::isfinite(target_interior) || target_interior <= 0. || target_interior > 1. || !std::isfinite(target_sigma) || target_sigma < 0. || target_sigma > maximum_smoothing_sigma)
                throw std::invalid_argument("Invalid infill modifier density settings");
            maximum_ratio = std::max(maximum_ratio, density_bound(modifier.material_density, target_sigma, target_adaptive, target_interior));
            std::vector<WallSlice> mask_slices;
            for (const auto &slice : modifier.slices)
                mask_slices.push_back({slice.bottom_z, slice.top_z, &slice.contours});
            // Use Orca's sliced modifier outline, including its top and bottom.
            std::mutex mask_progress_mutex;
            AdaptiveInfillDistanceField mask(mask_slices, {}, [&] {
                // Voxel cancellation probes may run concurrently on TBB workers.
                std::lock_guard<std::mutex> lock(mask_progress_mutex);
                report(progress, WallDistanceStage::Modifiers, double(modifier_index) / modifiers.size());
            }, {}, true);
            const bool own_smoothing = target_sigma != smoothing_sigma;
            const auto modifier_values = own_smoothing ? smoothed_distances(original, target_sigma, [&](auto, double) {
                report(progress, WallDistanceStage::Modifiers, double(modifier_index) / modifiers.size());
            }) : std::vector<float>{};
            const auto &target_values = own_smoothing ? modifier_values : smoothed;
            for (size_t first = 0; first < count; first += 4096) {
                report(progress, WallDistanceStage::Modifiers, (modifier_index + double(first) / count) / modifiers.size());
                std::vector<Vec3d> points;
                const size_t last = std::min(count, first + 4096);
                points.reserve(last - first);
                for (size_t i = first; i < last; ++i)
                    points.push_back(origin + spacing * Vec3d(double(i % size[0]), double((i / size[0]) % size[1]), double(i / (size[0] * size[1]))));
                const auto depths = mask.sample(points);
                for (size_t j = 0; j < points.size(); ++j) {
                    const double blend = modifier_blend(-depths[j], modifier.blend_width);
                    if (blend == 0.) continue;
                    const size_t i = first + j;
                    const double raw = original.node_value(i);
                    const double mean = target_values.empty() ? raw : double(target_values[i]);
                    const double density = std::min(1., modifier.material_density == 0. ? 0. : modifier.material_density * local_ratio(points[j], raw, mean, target_mode, target_adaptive, target_interior / modifier.material_density, target_gradient));
                    material[i] = float((1. - blend) * material[i] + blend * density);
                }
            }
            ++modifier_index;
        }
        if (!material.empty()) std::vector<float>().swap(smoothed);
        report(progress, WallDistanceStage::Complete, 1.);
    }

    double local_ratio(const Vec3d &point, double raw, double mean, AdaptiveInfillMode mode,
        bool graded, double interior_ratio, AdaptiveInfillGradient gradient) const
    {
        if (!graded) return InfillDensityField::density_ratio(raw, mean);
        const auto &shape = geometry[size_t(mode)];
        if (!shape || shape->empty()) return InfillDensityField::density_ratio(raw, mean);
        const bool planar = mode == AdaptiveInfillMode::NormalX || mode == AdaptiveInfillMode::NormalY || mode == AdaptiveInfillMode::NormalZ;
        if ((planar ? shape->distances().signed_distance(point) : raw) >= 0.)
            return InfillDensityField::density_ratio(raw, mean);
        const auto at_depth = [&](double depth) {
            if (mean == raw) return adaptive_density_scale(interior_ratio, gradient, depth);
            return InfillDensityField::density_ratio(raw, mean, depth > 0. ? -raw / depth : 0., interior_ratio, gradient);
        };
        if (mode == AdaptiveInfillMode::DistanceWarp || mode == AdaptiveInfillMode::SmoothBlend || mode == AdaptiveInfillMode::SteppedShells)
            return at_depth(shape->depth(point));
        AdaptiveInfillField::Radials radials;
        const size_t count = shape->radial(point, radials);
        if (count == 0) return InfillDensityField::density_ratio(raw, mean);
        double ratio = 0.;
        for (size_t j = 0; j < count; ++j) ratio += radials[j].weight * at_depth(1. - radials[j].t);
        return ratio;
    }

    double density(const Vec3d &point) const
    {
        if (original.empty()) return reference_density;
        const auto where = original.interpolation(point);
        if (!material.empty()) {
            double value = 0.;
            for (size_t j = 0; j < 8; ++j) value += material[where.indices[j]] * where.weights[j];
            return value;
        }
        const double raw = original.sample(where);
        double mean = raw;
        if (!smoothed.empty()) {
            mean = 0.;
            for (size_t j = 0; j < 8; ++j) mean += smoothed[where.indices[j]] * where.weights[j];
            if (where.outside > 0.) mean = std::min(AdaptiveInfillDistanceField::far_distance,
                std::hypot(std::max(0., mean), where.outside));
        }
        return std::min(1., reference_density * local_ratio(point, raw, mean, mode,
            mode != AdaptiveInfillMode::Disabled && interior_density != reference_density, interior_density / reference_density, gradient));
    }
};

InfillDensityField::InfillDensityField(std::vector<WallSlice> slices, double smoothing_sigma, const WallDistanceProgress &progress, AdaptiveInfillMode mode, double material_density,
                                     const std::vector<DensityModifier> &modifiers, int region_id, double interior_density, AdaptiveInfillGradient gradient, AdaptiveInfillFields geometry)
    : m_data(std::make_shared<Data>(std::move(slices), smoothing_sigma, progress, mode, material_density, modifiers, region_id, interior_density, gradient, std::move(geometry))) {}

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

double InfillDensityField::density_ratio(double distance, double sampled_mean, double maximum_depth, double interior_ratio, AdaptiveInfillGradient gradient)
{
    const double correction = std::abs(sampled_mean - distance);
    const double raw_depth = std::max(0., -distance);
    const double depth = std::max(0., raw_depth + (interior_ratio > 1. ? correction : -correction));
    const double t = maximum_depth > 0. ? std::clamp(depth / maximum_depth, 0., 1.) : 0.;
    const double target = adaptive_density_scale(interior_ratio, gradient, t);
    const double boost_depth = maximum_depth == 0. ? correction :
        correction - std::max(0., raw_depth - near_wall_distance);
    const double boost = (1. - minimum_density_ratio) *
        std::expm1(std::clamp(boost_depth, 0., near_wall_distance) / wall_falloff);
    return target + boost;
}

double InfillDensityField::density(const Vec3d &point) const { return m_data->density(point); }

double InfillDensityField::maximum_density() const { return m_data->maximum_ratio * m_data->reference_density; }

std::shared_ptr<const AdaptiveInfillField> InfillDensityCache::geometry(
    const std::function<std::vector<WallSlice>()> &load_slices, AdaptiveInfillMode mode,
    bool signed_distances, const WallDistanceProgress &progress)
{
    const int axis = mode == AdaptiveInfillMode::NormalX ? 0 : mode == AdaptiveInfillMode::NormalY ? 1 : mode == AdaptiveInfillMode::NormalZ ? 2 : -1;
    std::lock_guard<std::mutex> lock(m_mutex);
    auto &field = m_geometry[size_t(signed_distances)][size_t(mode)];
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
    const WallDistanceProgress &progress, const std::vector<DensityModifier> &modifiers, int region_id)
{
    auto shapes = density_geometry(*this, load_slices, s.mode, modifiers, region_id, progress);
    const bool composited = std::any_of(modifiers.begin(), modifiers.end(),
        [region_id](const auto &modifier) { return modifier.region_id == region_id; });
    if (s.sigma == 0. && !composited)
        return std::make_shared<InfillDensityField>(std::vector<WallSlice>{}, 0., progress,
            s.mode, s.material_density, modifiers, region_id, s.interior_density, s.gradient, std::move(shapes));
    std::lock_guard<std::mutex> lock(m_mutex);
    auto &field = m_fields[{s.mode, s.sigma, s.material_density, region_id, s.interior_density, s.gradient}];
    if (!field)
        field = std::make_shared<InfillDensityField>(std::vector<WallSlice>{}, s.sigma, progress,
            s.mode, s.material_density, modifiers, region_id, s.interior_density, s.gradient, std::move(shapes));
    return field;
}

void InfillDensityCache::clear()
{
    std::lock_guard<std::mutex> lock(m_mutex);
    m_fields.clear();
    m_geometry = {};
    m_distances.clear();
}

} // namespace Slic3r
