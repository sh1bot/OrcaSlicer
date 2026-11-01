#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include "../libslic3r.h"
#include "FillBase.hpp"
#include "../BoundingBox.hpp"
#include "../ExPolygon.hpp"
#include "../Point.hpp"
#include "../Polyline.hpp"
#include "../PrintConfig.hpp"
#include <cstdint>

namespace Slic3r {

// Distances used by adaptive infill. Signed sampling is negative inside.
// Slices must be ordered, non-overlapping, and their expolygons valid during construction.
class AdaptiveInfillDistanceField {
public:
    struct Slice { double bottom_z, top_z; const ExPolygons *expolygons; };
    enum class Method { Voxel, SliceSlabs };
    enum class Stage { Distance, Encoding, Complete };
    using Progress = std::function<void(Stage, double)>;
    struct Interpolation {
        std::array<size_t, 8> indices {};
        std::array<double, 8> weights {};
        double outside = 0.;
    };
    static constexpr double far_distance = 1000.;
    AdaptiveInfillDistanceField() = default;
    explicit AdaptiveInfillDistanceField(const std::vector<Slice> &slices, const BoundingBox &bbox = {},
        const std::function<void()> &throw_if_canceled = {}, const Progress &progress = {},
        bool signed_distances = false, int normal_axis = -1, Method method = Method::Voxel);
    Method method() const { return m_method; }
    // Exterior samples are zero when disabled; the voxel backend also skips its exterior transform.
    bool has_signed_distances() const { return m_signed; }
    bool empty() const { return node_count() == 0; }
    double signed_distance(const Vec3d &point) const { return sample(interpolation(point)); }
    double grid_spacing() const { return spacing; }
    const Vec3d &grid_origin() const { return origin; }
    const std::array<size_t, 3> &grid_size() const { return size; }
    size_t node_count() const { return original->distance.empty() ? original->squared.size() : original->distance.size(); }
    // Voxel-transform costs; slab distances are converted to the same grid units.
    float squared_depth(size_t index) const {
        if (!is_inside(index)) return 0.f;
        return original->distance.empty() ? original->squared[index] :
            float(std::pow(decoded[original->distance[index]] / spacing + 0.5, 2));
    }
    bool is_inside(size_t index) const { return original->inside[index / 8] & (uint8_t(1) << (index % 8)); }
    double node_value(size_t index) const {
        if (!m_signed && !is_inside(index)) return 0.;
        const double d = original->distance.empty() ?
            std::min(far_distance, std::max(0., std::sqrt(double(original->squared[index])) - 0.5) * spacing) :
            decoded[original->distance[index]];
        return is_inside(index) ? -d : d;
    }
    // Reuse the same interpolation weights for fields built on this grid.
    Interpolation interpolation(const Vec3d &point) const;
    double sample(const Interpolation &where) const;
    std::vector<double> sample(const std::vector<Vec3d> &points) const;
private:
    void build_slice_slabs(const std::vector<Slice> &slices, const BoundingBox &bbox,
        const std::function<void()> &throw_if_canceled, const Progress &progress);
    static constexpr double distance_scale = 6.;
    Method m_method = Method::Voxel;
    std::array<double, 256> decoded {};
    bool m_signed = false;
    double spacing = 1.;
    Vec3d origin = Vec3d::Zero();
    std::array<size_t, 3> size {{0, 0, 0}};
    struct Field { std::vector<uint8_t> inside, distance; std::vector<float> squared; };
    std::shared_ptr<Field> original = std::make_shared<Field>();
};

// Radial coordinate inside the lobes of the bodies of an object, sampled from its slices: 0 at the center of a
// lobe, 1 at its surface. Lobes are parts of a body separated by a neck, like two spheres united. In the 2D modes,
// every section normal to the axis has its own bodies and lobes, and distances are measured within the section. The
// modes following the distance to the surface use the depth of every point instead; Distance warp also the lobes.
class AdaptiveInfillField
{
public:
    using Slice = AdaptiveInfillDistanceField::Slice;

    struct Radial
    {
        Vec3d  center;
        double t;
        float  weight;
    };

    // Slices sorted by z, in the XY coordinates of the fill and the print Z.
    AdaptiveInfillField(const std::vector<Slice> &slices, const BoundingBox &bbox, AdaptiveInfillMode mode,
                    const std::function<void()> &throw_if_canceled, std::shared_ptr<const AdaptiveInfillDistanceField> distances = {});

    // Lobes blended near the sides between them, from a body or from each of the two sections around a point.
    static constexpr size_t MaxMorph = 4;
    using Radials                    = std::array<Radial, 2 * MaxMorph>;

    // Without a body, as when the object is thinner than the grid cells.
    bool empty() const { return m_bodies.empty(); }

    // Radial coordinates of pt in unscaled coordinates towards the lobe it belongs to, and towards the neighbouring
    // lobes near the sides between them, with weights summing to 1. In the 2D modes, those of the two sections around
    // pt. Returns their count.
    size_t radial(const Vec3d &pt, Radials &out) const;

    // Axis normal to the sections in the 2D modes, -1 in the modes graded in 3D.
    int              axis() const { return m_axis; }
    AdaptiveInfillMode mode() const { return m_mode; }

    // In the modes following the distance to the surface: depth relative to the deepest point of the body, from 0 at
    // the surface to 1.
    double depth(const Vec3d &pt) const;
    const AdaptiveInfillDistanceField &distances() const { return *m_distances; }

private:
    struct Lobe
    {
        Vec3d  center;
        double depth;
        // Distance from the center to the surface on a latitude-longitude grid of directions.
        std::vector<float> reach;
        // Distance warp: mean depth over the ball along each direction, sampled up to the reach.
        std::vector<float> mean_depth;
    };

    struct Body
    {
        size_t first_lobe;
        size_t lobes;
        // Distance from the deepest point to the surface.
        double depth;
    };

    double radial(const Lobe &lobe, const Vec3d &pt) const;
    size_t body_radial(size_t node, const Vec3d &pt, float weight, Radial *out) const;
    // Offset of pt from a center, within the section in the 2D modes.
    Vec3d  offset(const Vec3d &pt, const Vec3d &center) const;
    int    directions() const { return m_axis < 0 ? Polar * Azimuth : Azimuth; }

    // Directions of the reach of a lobe, on a latitude-longitude grid, or a circle in the 2D modes.
    static constexpr int Polar   = 24;
    static constexpr int Azimuth = 48;

    std::shared_ptr<const AdaptiveInfillDistanceField> m_distances;
    AdaptiveInfillMode  m_mode;
    int               m_axis;
    Vec3d             m_origin;
    double            m_cell;
    Vec3i32           m_size;
    // Nearest body of every grid node.
    std::vector<int>  m_body;
    std::vector<Body> m_bodies;
    std::vector<Lobe> m_lobes;
    // In the 2D modes, the nearest section with a body to every section.
    std::vector<int>  m_section;
    // In the modes following the distance to the surface, the depth of every grid node.
    std::vector<float> m_depth;
};

using AdaptiveInfillFieldPtr = std::shared_ptr<const AdaptiveInfillField>;
// A field for every adaptive mode in use, indexed by the mode.
using AdaptiveInfillFields = std::array<AdaptiveInfillFieldPtr, size_t(AdaptiveInfillMode::Count)>;

// Relative density/frequency between the surface (depth 0) and interior (depth 1).
inline double adaptive_density_scale(double ratio, AdaptiveInfillGradient gradient, double depth)
{
    switch (gradient) {
    case AdaptiveInfillGradient::Quadratic: return 1. + (ratio - 1.) * depth * depth;
    case AdaptiveInfillGradient::Exponential: return std::pow(ratio, depth);
    default: return 1. + (ratio - 1.) * depth;
    }
}

struct AdaptiveTpms
{
    // Implicit TPMS equation with a period of 2 PI.
    float (*equation)(float x, float y, float z);
    // Pattern frequencies at the surface and at the center, in radians per mm.
    double               surface_frequency;
    double               interior_frequency;
    AdaptiveInfillGradient gradient;
    // Optional material-density consumer, returning radians per mm in object coordinates.
    std::function<double(const Vec3d &)> frequency_at;
};

// Infill lines in the fill frame, the object frame rotated by -angle; z is the print_z of the layer.
Polylines make_adaptive_tpms(const AdaptiveTpms &tpms, const AdaptiveInfillField &field, BoundingBox bbox,
                             coordf_t z, coordf_t layer_height, coordf_t spacing, float angle);

struct AdaptiveInfillShell
{
    float      density;
    ExPolygons expolygons;
};

// Stepped shells: the parts of an expolygon in the object frame at each density, from the surface inwards; z is the
// middle of the layer.
std::vector<AdaptiveInfillShell> make_adaptive_infill_shells(const AdaptiveInfillField &field, const ExPolygon &expolygon, coordf_t z,
                                        float surface_density, float interior_density, AdaptiveInfillGradient gradient);

// Stepped shells: fills every shell with fill_shell at its density, each shrunk by half a line like a filled region,
// so the lines connected along the boundaries of two shells don't overlap.
void fill_adaptive_infill_shells(const AdaptiveInfillField &field, const ExPolygon &expolygon, coordf_t z, const FillParams &params, coordf_t spacing,
                      const std::function<void(const FillParams &, const ExPolygon &)> &fill_shell);

} // namespace Slic3r
