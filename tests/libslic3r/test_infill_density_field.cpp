#include <cstddef>
#include "libslic3r/Fill/AdaptiveInfillField.hpp"
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"
#include "catch2/generators/catch_generators.hpp"
#include "libslic3r/libslic3r.h"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Polyline.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/BoundingBox.hpp"
#include <algorithm>
#include <limits>
#include <catch2/catch_all.hpp>
#include "libslic3r/InfillDensityField.hpp"
#include "libslic3r/InfillDensityModifier.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include <cmath>
#include <stdexcept>
#include <vector>

using namespace Slic3r;

namespace {
ExPolygon rectangle(double x0, double y0, double x1, double y1)
{
    return ExPolygon(Points {Point::new_scale(x0, y0), Point::new_scale(x1, y0),
                             Point::new_scale(x1, y1), Point::new_scale(x0, y1)});
}
float gyroid(float x, float y, float z)
{
    return std::sin(x) * std::cos(y) + std::sin(y) * std::cos(z) + std::sin(z) * std::cos(x);
}
float diamond(float x, float y, float z)
{
    return std::cos(x) * std::cos(y) * std::cos(z) - std::sin(x) * std::sin(y) * std::sin(z);
}
float fischer_koch(float x, float y, float z)
{
    return std::cos(2 * x) * std::sin(y) * std::cos(z) + std::cos(2 * y) * std::sin(z) * std::cos(x) +
           std::cos(2 * z) * std::sin(x) * std::cos(y);
}
double length(const Polylines &paths)
{
    double total = 0.;
    for (const auto &path : paths) total += unscale<double>(path.length());
    return total;
}
}

TEST_CASE("Shared density scalar queries agree with Voronoi batch queries", "[InfillDensityField]")
{
    const ExPolygons model {rectangle(-20., -20., 20., 20.)};
    const InfillDensityField field({{0., 40., &model}}, 8.);
    const std::vector<Vec3d> points {Vec3d(0., 0., 20.), Vec3d(18., 18., 2.), Vec3d(19., 0., 20.)};
    const auto values = field.sample_fields(points);
    for (size_t i = 0; i < points.size(); ++i) {
        REQUIRE_THAT(field.density(points[i]), Catch::Matchers::WithinAbs(values[i].density * 0.25, 1e-12));
        REQUIRE(field.maximum_density() >= field.density(points[i]));
    }
    REQUIRE(field.density(points[1]) > 0.25);
}

TEST_CASE("Density evaluators share geometry while retaining their own settings", "[InfillDensityField]")
{
    const ExPolygons model {rectangle(-10., -10., 10., 10.)};
    size_t loads = 0;
    const auto slices = [&] { ++loads; return std::vector<WallSlice> {{0., 20., &model}}; };
    InfillDensityCache cache;
    InfillDensitySettings settings;
    const auto first = cache.get(slices, settings);
    REQUIRE_THAT(cache.get(slices, settings)->density(Vec3d(0., 0., 10.)), Catch::Matchers::WithinAbs(first->density(Vec3d(0., 0., 10.)), 1e-12));
    REQUIRE(loads == 1);
    settings.interior_density = 0.1;
    REQUIRE(cache.get(slices, settings) != first);
    REQUIRE(loads == 1);
    REQUIRE(cache.get(slices, settings, {}, {}, 3) != first);
    REQUIRE(loads == 1);
    cache.clear();
    REQUIRE(cache.get(slices, settings) != first);
    REQUIRE(loads == 2);
}

TEST_CASE("Shared material density blends a modifier inside its boundary and preserves empty interiors", "[InfillDensityField]")
{
    const double target = GENERATE(0., 0.5);
    const ExPolygons model {rectangle(-20., -20., 20., 20.)};
    DensityModifier modifier {3, target, 4., {{0., 40., {rectangle(0., -20., 20., 20.)}}}};
    if (target > 0.) modifier.interior_density = target;
    const InfillDensityField field({{0., 40., &model}}, 0., {}, AdaptiveInfillMode::Disabled, 0.1, {modifier}, 3);
    REQUIRE_THAT(field.density(Vec3d(-4., 0., 20.)), Catch::Matchers::WithinRel(0.1, 0.022));
    // Planar voxel boundaries may shift by half a cell; smoothstep has a maximum slope of 1.5.
    const double tolerance = 0.025 * std::max(0.1, target) + 1.5 * std::abs(target - 0.1) * 0.5 * field.grid_spacing() / 4.;
    REQUIRE_THAT(field.density(Vec3d(2., 0., 20.)), Catch::Matchers::WithinAbs((0.1 + target) / 2., tolerance));
    REQUIRE_THAT(field.density(Vec3d(8., 0., 20.)), Catch::Matchers::WithinAbs(target, 0.012));
}

TEST_CASE("Smooth blend accepts constant shared densities without changing lattice phase", "[InfillDensityField][AdaptiveInfillField]")
{
    const auto equation = GENERATE(&gyroid, &diamond, &fischer_koch);
    const double frequency = GENERATE(0.4, 1., 2.5);
    const ExPolygons model {rectangle(-12., -12., 12., 12.)};
    const AdaptiveInfillField field({{0., 24., &model}}, get_extents(model), AdaptiveInfillMode::SmoothBlend, [] {});
    const BoundingBox bounds = get_extents(model);
    const AdaptiveTpms reference {equation, frequency, frequency, AdaptiveInfillGradient::Linear};
    AdaptiveTpms shared {equation, 1., 1., AdaptiveInfillGradient::Linear};
    shared.frequency_at = [frequency](const Vec3d &) { return frequency; };
    const auto expected = make_adaptive_tpms(reference, field, bounds, 12., 0.2, 0.45, 0.3f);
    const auto actual = make_adaptive_tpms(shared, field, bounds, 12., 0.2, 0.45, 0.3f);
    REQUIRE_FALSE(expected.empty());
    REQUIRE(actual.size() == expected.size());
    REQUIRE_THAT(length(actual), Catch::Matchers::WithinAbs(length(expected), 0.02));
}

TEST_CASE("Smooth blend consumes fused density outside the root density range", "[InfillDensityField][AdaptiveInfillField]")
{
    const auto equation = GENERATE(&gyroid, &diamond, &fischer_koch);
    const double target = GENERATE(0., 0.5);
    const ExPolygons model {rectangle(-15., -15., 15., 15.)};
    DensityModifier modifier {3, target, 3., {{0., 30., model}}};
    if (target > 0.) modifier.interior_density = target;
    const InfillDensityField material({{0., 30., &model}}, 0., {}, AdaptiveInfillMode::Disabled, 0.1, {modifier}, 3);
    const AdaptiveInfillField field({{0., 30., &model}}, get_extents(model), AdaptiveInfillMode::SmoothBlend, [] {});
    AdaptiveTpms tpms {equation, 0.4, 0.2, AdaptiveInfillGradient::Linear};
    tpms.frequency_at = [&material](const Vec3d &point) { return 4. * material.density(point); };
    const BoundingBox center(Point::new_scale(-8., -8.), Point::new_scale(8., 8.));
    const auto paths = make_adaptive_tpms(tpms, field, center, 15., 0.2, 0.45, 0.f);
    if (target == 0.) REQUIRE(paths.empty());
    else {
        REQUIRE_FALSE(paths.empty());
        const auto ordinary = make_adaptive_tpms({equation, 0.4, 0.4, AdaptiveInfillGradient::Linear}, field, center, 15., 0.2, 0.45, 0.f);
        REQUIRE(length(paths) > 2. * length(ordinary));
    }
}

TEST_CASE("Density modifier fusion is offered only by compatible parent patterns and modes", "[InfillDensityField]")
{
    const auto pattern = GENERATE(ipVoronoi, ipGyroid, ipTpmsD, ipTpmsFK, ipCubic);
    const auto mode = GENERATE(AdaptiveInfillMode::Disabled, AdaptiveInfillMode::SmoothBlend, AdaptiveInfillMode::Lobes, AdaptiveInfillMode::SteppedShells);
    const double density = GENERATE(0., 15., 100.);
    const bool expected = density > 0. && (pattern == ipVoronoi ||
        (pattern != ipCubic && mode == AdaptiveInfillMode::SmoothBlend && density < 100.));
    REQUIRE(supports_density_field(pattern, mode, density) == expected);
}

TEST_CASE("Smooth blend regional sampling preserves one spatial phase across a modifier", "[InfillDensityField][AdaptiveInfillField]")
{
    constexpr double bead_width = 0.45;
    const auto equation = GENERATE(&gyroid, &diamond, &fischer_koch);
    const ExPolygons model {rectangle(-15., -15., 15., 15.)};
    DensityModifier modifier {3, 0.5, 8., {{0., 30., {rectangle(0., -15., 15., 15.)}}}};
    modifier.interior_density = modifier.material_density;
    const InfillDensityField material({{0., 30., &model}}, 0., {}, AdaptiveInfillMode::Disabled, 0.1, {modifier}, 3);
    const AdaptiveInfillField radial({{0., 30., &model}}, get_extents(model), AdaptiveInfillMode::SmoothBlend, [] {});
    AdaptiveTpms tpms {equation, 0.4, 0.2, AdaptiveInfillGradient::Linear};
    tpms.frequency_at = [&material](const Vec3d &point) { return 4. * material.density(point); };
    const auto whole = make_adaptive_tpms(tpms, radial, get_extents(model), 15., 0.2, bead_width, 0.3f);
    for (bool right : {false, true}) {
        const ExPolygon half = rectangle(right ? 0. : -15., -15., right ? 15. : 0., 15.);
        const ExPolygons window {rectangle(right ? 2. : -12., -12., right ? 12. : -2., 12.)};
        const auto part = make_adaptive_tpms(tpms, radial, get_extents(half), 15., 0.2, bead_width, 0.3f);
        const auto expected = intersection_pl(whole, window);
        const auto actual = intersection_pl(part, window);
        REQUIRE_FALSE(expected.empty());
        REQUIRE_FALSE(actual.empty());
        // Regional marching-square closures may differ; the phase stays within half a bead.
        for (bool reverse : {false, true}) {
            const auto &a = reverse ? expected : actual;
            const auto &b = reverse ? actual : expected;
            for (const auto &line : a)
                for (const Point &point : line.points) {
                    double nearest = std::numeric_limits<double>::infinity();
                    for (const auto &other : b) nearest = std::min(nearest, point.distance_to(point.projection_onto(other)));
                    REQUIRE(unscale<double>(nearest) <= 0.5 * bead_width);
                }
        }
    }
}

TEST_CASE("Default density smoothing is disabled without disabling grading", "[InfillDensityField]")
{
    const ExPolygons model {rectangle(-20., -20., 20., 20.)};
    const InfillDensityField defaults({{0., 40., &model}});
    bool smoothed = false;
    const InfillDensityField explicit_zero({{0., 40., &model}}, 0., [&](auto stage, double) {
        smoothed |= stage == WallDistanceStage::SmoothDistance || stage == WallDistanceStage::SmoothReinforcement;
    });
    REQUIRE_FALSE(smoothed);
    REQUIRE_THAT(InfillDensitySettings{}.sigma, Catch::Matchers::WithinAbs(0., 1e-12));
    for (const Vec3d &p : {Vec3d(0., 0., 20.), Vec3d(19., 19., 1.)})
        REQUIRE_THAT(defaults.density(p), Catch::Matchers::WithinAbs(explicit_zero.density(p), 1e-12));
    REQUIRE(defaults.density(Vec3d(0., 0., 20.)) < defaults.density(Vec3d(19., 19., 1.)));
}

TEST_CASE("Modifier smoothing rejects unbounded kernels", "[InfillDensityField]")
{
    const ExPolygons model {rectangle(-2., -2., 2., 2.)};
    DensityModifier modifier {3, 0.5, 0., {{0., 4., model}}};
    modifier.smoothing_sigma = InfillDensityField::maximum_smoothing_sigma + 1.;
    REQUIRE_THROWS_AS(InfillDensityField({{0., 4., &model}}, 0., {}, AdaptiveInfillMode::Disabled, 0.1, {modifier}, 3), std::invalid_argument);
}

TEST_CASE("Density cache shares raw distances across settings and keeps returned fields alive", "[InfillDensityField]")
{
    const ExPolygons model {rectangle(-4., -4., 4., 4.)};
    InfillDensityCache cache;
    unsigned builds = 0;
    const auto slices = [&] { return std::vector<WallSlice>{{0., 8., &model}}; };
    const auto progress = [&](WallDistanceStage stage, double fraction) {
        if (stage == WallDistanceStage::Distance && fraction == 1.) ++builds;
    };
    InfillDensitySettings settings;
    const auto first = cache.get(slices, settings, progress);
    REQUIRE(builds > 0);
    const unsigned initial = builds;
    settings.material_density = 0.5;
    const auto second = cache.get(slices, settings, progress);
    REQUIRE(second != first);
    REQUIRE(builds == initial);
    const double before = first->density(Vec3d(0., 0., 4.));
    cache.clear();
    REQUIRE_THAT(first->density(Vec3d(0., 0., 4.)), Catch::Matchers::WithinAbs(before, 1e-12));
    cache.get(slices, settings, progress);
    REQUIRE(builds > initial);
}

TEST_CASE("Fused density modifiers inherit or override adaptive density explicitly", "[InfillDensityField][Regression]")
{
    const bool parent_adaptive = GENERATE(false, true);
    const int override_mode = GENERATE(-1, 0, 1);
    const ExPolygons model {rectangle(-20., -20., 20., 20.)};
    DensityModifier modifier {3, 0.3, 0., {{0., 40., model}}};
    modifier.interior_density = 0.05;
    if (override_mode >= 0) modifier.mode = override_mode ? AdaptiveInfillMode::SmoothBlend : AdaptiveInfillMode::Disabled;
    const InfillDensityField field({{0., 40., &model}}, 0., {}, parent_adaptive ? AdaptiveInfillMode::SmoothBlend : AdaptiveInfillMode::Disabled, 0.15, {modifier}, 3);
    const bool graded = modifier.mode.value_or(parent_adaptive ? AdaptiveInfillMode::SmoothBlend : AdaptiveInfillMode::Disabled) != AdaptiveInfillMode::Disabled;
    const double center = field.density(Vec3d(0., 0., 20.));
    const double near_wall = field.density(Vec3d(16., 0., 20.));
    if (graded) {
        REQUIRE(center < near_wall);
        REQUIRE_THAT(center, Catch::Matchers::WithinRel(0.05, 0.04));
    } else {
        REQUIRE_THAT(center, Catch::Matchers::WithinRel(0.3, 0.025));
        REQUIRE_THAT(near_wall, Catch::Matchers::WithinRel(0.3, 0.025));
    }
}

TEST_CASE("Adaptive geometry shares raw distances across modes and density settings", "[InfillDensityField][Regression]")
{
    const ExPolygons model {rectangle(-8., -8., 8., 8.)};
    const auto slices = [&] { return std::vector<WallSlice> {{0., 16., &model}}; };
    InfillDensityCache cache;
    const auto blend = cache.geometry(slices, AdaptiveInfillMode::SmoothBlend);
    const auto warp = cache.geometry(slices, AdaptiveInfillMode::DistanceWarp);
    REQUIRE(blend == cache.geometry(slices, AdaptiveInfillMode::SmoothBlend));
    REQUIRE(&blend->distances() == &warp->distances());
    REQUIRE(blend->distances().has_signed_distances());
    InfillDensitySettings settings;
    settings.mode = AdaptiveInfillMode::SmoothBlend;
    const auto first = cache.get(slices, settings, {}, {}, -1);
    settings.material_density = 0.3;
    REQUIRE(cache.get(slices, settings, {}, {}, -1) != first);
    REQUIRE(blend == cache.geometry(slices, AdaptiveInfillMode::SmoothBlend));
    cache.clear();
    REQUIRE(cache.geometry(slices, AdaptiveInfillMode::SmoothBlend) != blend);
}


TEST_CASE("Voronoi adaptive density follows every shared spatial profile", "[InfillDensityField][Voronoi]")
{
    const auto mode = GENERATE(AdaptiveInfillMode::DistanceWarp, AdaptiveInfillMode::SmoothBlend,
        AdaptiveInfillMode::SteppedShells, AdaptiveInfillMode::Lobes,
        AdaptiveInfillMode::NormalX, AdaptiveInfillMode::NormalY, AdaptiveInfillMode::NormalZ);
    const auto gradient = GENERATE(AdaptiveInfillGradient::Linear, AdaptiveInfillGradient::Quadratic,
        AdaptiveInfillGradient::Exponential);
    const ExPolygons model {rectangle(-16., -10., 16., 10.)};
    const auto slices = [&] { return std::vector<WallSlice>{{0., 24., &model}}; };
    InfillDensityCache cache;
    const InfillDensitySettings settings {mode, 0., 0.3, 0.03, gradient};
    const auto material = cache.get(slices, settings);
    const auto shape = cache.geometry(slices, mode);
    const auto raw = cache.geometry(slices, AdaptiveInfillMode::SmoothBlend);
    for (const Vec3d &point : {Vec3d(0., 0., 12.), Vec3d(8., 3., 7.), Vec3d(-5., -6., 18.)}) {
        // Compare at a density-grid node, before interpolation can mix nonlinear values.
        const auto &grid = raw->distances();
        const Vec3d node = grid.grid_origin() + grid.grid_spacing() *
            ((point - grid.grid_origin()) / grid.grid_spacing()).array().round().matrix();
        double scale = 0.;
        if (mode == AdaptiveInfillMode::DistanceWarp || mode == AdaptiveInfillMode::SmoothBlend || mode == AdaptiveInfillMode::SteppedShells)
            scale = adaptive_density_scale(0.1, gradient, shape->depth(node));
        else {
            AdaptiveInfillField::Radials radials;
            const auto count = shape->radial(node, radials);
            REQUIRE(count > 0);
            for (size_t i = 0; i < count; ++i)
                scale += radials[i].weight * adaptive_density_scale(0.1, gradient, 1. - radials[i].t);
        }
        REQUIRE_THAT(material->density(node), Catch::Matchers::WithinRel(0.3 * scale, 0.025));
    }
    REQUIRE(cache.get(slices, settings)->density(Vec3d(0., 0., 12.)) == material->density(Vec3d(0., 0., 12.)));
    auto other = settings;
    other.mode = AdaptiveInfillMode::Disabled;
    REQUIRE(cache.get(slices, other) != material);
    REQUIRE_THAT(cache.get(slices, other)->density(Vec3d(0., 0., 12.)), Catch::Matchers::WithinRel(0.3, 0.025));
}

TEST_CASE("Voronoi distance application modes produce the same continuous density", "[InfillDensityField][Voronoi]")
{
    const ExPolygons model {rectangle(-16., -10., 16., 10.), rectangle(24., -6., 36., 6.)};
    const auto slices = [&] { return std::vector<WallSlice>{{0., 24., &model}}; };
    InfillDensityCache cache;
    InfillDensitySettings settings {AdaptiveInfillMode::SmoothBlend, 0., 0.3, 0.03, AdaptiveInfillGradient::Quadratic};
    const auto reference = cache.get(slices, settings);
    for (const auto mode : {AdaptiveInfillMode::DistanceWarp, AdaptiveInfillMode::SteppedShells}) {
        settings.mode = mode;
        const auto field = cache.get(slices, settings);
        for (const Vec3d &point : {Vec3d(0., 0., 12.), Vec3d(8., 3., 7.), Vec3d(30., 0., 12.)})
            REQUIRE_THAT(field->density(point), Catch::Matchers::WithinAbs(reference->density(point), 1e-12));
    }
    REQUIRE_THAT(reference->density(Vec3d(0., 0., 12.)), Catch::Matchers::WithinRel(0.03, 0.025));
    REQUIRE_THAT(reference->density(Vec3d(30., 0., 12.)), Catch::Matchers::WithinRel(0.03, 0.025));
}

TEST_CASE("Fused modifiers retain their selected adaptive spatial mode", "[InfillDensityField][Voronoi]")
{
    const auto mode = GENERATE(AdaptiveInfillMode::Lobes, AdaptiveInfillMode::NormalX,
        AdaptiveInfillMode::NormalY, AdaptiveInfillMode::NormalZ);
    const bool inherit = GENERATE(false, true);
    const ExPolygons model {rectangle(-16., -10., 16., 10.)};
    const auto slices = [&] { return std::vector<WallSlice>{{0., 24., &model}}; };
    DensityModifier modifier {3, 0.3, 0., {{0., 24., model}}};
    modifier.interior_density = 0.03;
    if (!inherit) modifier.mode = mode;
    InfillDensityCache parent_cache, reference_cache;
    InfillDensitySettings settings {inherit ? mode : AdaptiveInfillMode::Disabled, 0., 0.15, 0.05, AdaptiveInfillGradient::Linear};
    const auto field = parent_cache.get(slices, settings, {}, {modifier}, 3);
    settings.mode = mode;
    settings.material_density = 0.3;
    settings.interior_density = 0.03;
    const auto reference = reference_cache.get(slices, settings);
    for (const Vec3d &point : {Vec3d(0., 0., 12.), Vec3d(8., 3., 7.), Vec3d(-5., -6., 18.)})
        REQUIRE_THAT(field->density(point), Catch::Matchers::WithinRel(reference->density(point), 2. * std::numeric_limits<float>::epsilon()));
}

TEST_CASE("Normal density profiles stay constant along their selected axis", "[InfillDensityField][Voronoi]")
{
    const int axis = GENERATE(0, 1, 2);
    const auto mode = axis == 0 ? AdaptiveInfillMode::NormalX : axis == 1 ? AdaptiveInfillMode::NormalY : AdaptiveInfillMode::NormalZ;
    const ExPolygons model {rectangle(-16., -10., 16., 10.)};
    InfillDensityField field({{0., 24., &model}}, 0., {}, mode, 0.3, {}, -1, 0.03);
    const Vec3d center(0., 0., 12.), extent(32., 20., 24.);
    Vec3d along = center, across = center;
    along[axis] += 0.25 * extent[axis];
    across[(axis + 1) % 3] += 0.3 * extent[(axis + 1) % 3];
    REQUIRE_THAT(field.density(along), Catch::Matchers::WithinAbs(field.density(center), 0.001));
    REQUIRE(field.density(across) > field.density(center) + 0.03);
}

TEST_CASE("Infill density caches keep slice-slab and voxel backends separate", "[InfillDensityField]")
{
    const ExPolygons model {rectangle(-10., -10., 10., 10.)};
    const auto slices = [&] { return std::vector<WallSlice>{{0., 20., &model}}; };
    InfillDensityCache cache;
    const auto voxel = cache.geometry(slices, AdaptiveInfillMode::SmoothBlend, true);
    const auto slab = cache.geometry(slices, AdaptiveInfillMode::SmoothBlend, true, {}, true);
    REQUIRE(voxel->distances().method() == AdaptiveInfillDistanceField::Method::Voxel);
    REQUIRE(slab->distances().method() == AdaptiveInfillDistanceField::Method::SliceSlabs);
    REQUIRE(cache.geometry(slices, AdaptiveInfillMode::SmoothBlend, true) == voxel);
    REQUIRE(cache.geometry(slices, AdaptiveInfillMode::SmoothBlend, true, {}, true) == slab);
    InfillDensitySettings settings;
    settings.mode = AdaptiveInfillMode::SmoothBlend;
    const auto first = cache.get(slices, settings);
    settings.slice_slabs = true;
    const auto second = cache.get(slices, settings);
    REQUIRE(first != second);
    REQUIRE(cache.get(slices, settings)->density(Vec3d(0., 0., 4.)) == second->density(Vec3d(0., 0., 4.)));
    const auto planar = cache.geometry(slices, AdaptiveInfillMode::NormalZ, true, {}, true);
    REQUIRE(planar->distances().method() == AdaptiveInfillDistanceField::Method::Voxel);
}

TEST_CASE("Adaptive density evaluates the upstream gradient at the query point", "[InfillDensityField]")
{
    const ExPolygons model {rectangle(-10., -10., 10., 10.)};
    const auto slices = [&] { return std::vector<WallSlice>{{0., 20., &model}}; };
    InfillDensityCache cache;
    InfillDensitySettings settings;
    settings.mode = AdaptiveInfillMode::SmoothBlend;
    settings.material_density = 0.3;
    settings.interior_density = 0.06;
    settings.gradient = AdaptiveInfillGradient::Exponential;
    const auto field = cache.get(slices, settings);
    const auto shape = cache.geometry(slices, settings.mode);
    const Vec3d point(1.234, 3.456, 8.765);
    const double expected = settings.material_density * adaptive_density_scale(
        settings.interior_density / settings.material_density, settings.gradient, shape->depth(point));
    REQUIRE_THAT(field->density(point), Catch::Matchers::WithinRel(expected, 1e-12));
}

TEST_CASE("Uniform material density retains its floating point value", "[InfillDensityField]")
{
    const ExPolygons model {rectangle(-10., -10., 10., 10.)};
    const auto slices = [&] { return std::vector<WallSlice>{{0., 20., &model}}; };
    InfillDensityCache cache;
    InfillDensitySettings settings;
    settings.material_density = 0.123456789;
    const auto field = cache.get(slices, settings);
    REQUIRE_THAT(field->density(Vec3d(1.234, 3.456, 8.765)),
        Catch::Matchers::WithinRel(settings.material_density, 1e-12));
}

TEST_CASE("Zero corner sigma preserves the upstream adaptive profile near floors", "[InfillDensityField][Regression]")
{
    const auto mode = GENERATE(AdaptiveInfillMode::SmoothBlend, AdaptiveInfillMode::Lobes,
        AdaptiveInfillMode::NormalX, AdaptiveInfillMode::NormalY, AdaptiveInfillMode::NormalZ);
    const ExPolygons model {rectangle(-10., -10., 10., 10.)};
    const auto slices = [&] { return std::vector<WallSlice>{{0., 20., &model}}; };
    InfillDensityCache cache;
    InfillDensitySettings settings;
    settings.mode = mode;
    settings.material_density = 0.3;
    settings.interior_density = 0.06;
    settings.gradient = AdaptiveInfillGradient::Exponential;
    const auto field = cache.get(slices, settings);
    const auto shape = cache.geometry(slices, mode);
    const Vec3d point(1.234, 3.456, 0.25);
    const auto at_depth = [&](double depth) {
        return adaptive_density_scale(settings.interior_density / settings.material_density, settings.gradient, depth);
    };
    double ratio = 1.;
    if (shape->distances().signed_distance(point) < 0.) {
        if (mode == AdaptiveInfillMode::SmoothBlend) ratio = at_depth(shape->depth(point));
        else {
            AdaptiveInfillField::Radials radials;
            const size_t count = shape->radial(point, radials);
            REQUIRE(count > 0);
            ratio = 0.;
            for (size_t i = 0; i < count; ++i) ratio += radials[i].weight * at_depth(1. - radials[i].t);
        }
    }
    REQUIRE_THAT(field->density(point), Catch::Matchers::WithinRel(settings.material_density * ratio, 1e-12));
}

TEST_CASE("Object profile arrays expose the same geometry used by density evaluation", "[InfillDensityField][Regression]")
{
    const ExPolygons model {rectangle(-8., -8., 8., 8.)};
    size_t loads = 0;
    const auto slices = [&] { ++loads; return std::vector<WallSlice>{{0., 16., &model}}; };
    InfillDensityCache cache;
    const auto &profiles = cache.profiles(true);
    const auto shape = cache.geometry(slices, AdaptiveInfillMode::SmoothBlend);
    REQUIRE(profiles[size_t(AdaptiveInfillMode::SmoothBlend)] == shape);
    const size_t initial_loads = loads;
    InfillDensitySettings settings;
    settings.mode = AdaptiveInfillMode::SmoothBlend;
    cache.get(slices, settings);
    REQUIRE(loads == initial_loads);
    REQUIRE(profiles[size_t(settings.mode)] == shape);
    cache.clear();
    REQUIRE_FALSE(profiles[size_t(settings.mode)]);
    const auto rebuilt = cache.geometry(slices, settings.mode);
    REQUIRE(rebuilt != shape);
    REQUIRE(profiles[size_t(settings.mode)] == rebuilt);
}
