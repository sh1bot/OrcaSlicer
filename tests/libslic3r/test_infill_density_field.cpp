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
#include "libslic3r/ClipperUtils.hpp"
#include <cmath>
#include <vector>

using namespace Slic3r;

namespace {
ExPolygon rectangle(double x0, double y0, double x1, double y1)
{
    return ExPolygon(Points {Point::new_scale(x0, y0), Point::new_scale(x1, y0),
                             Point::new_scale(x1, y1), Point::new_scale(x0, y1)});
}
}


TEST_CASE("Shared density cache reuses only fields with matching settings", "[InfillDensityField]")
{
    const ExPolygons model {rectangle(-10., -10., 10., 10.)};
    size_t loads = 0;
    const auto slices = [&] { ++loads; return std::vector<WallSlice> {{0., 20., &model}}; };
    InfillDensityCache cache;
    InfillDensitySettings settings;
    const auto first = cache.get(slices, settings);
    REQUIRE(cache.get(slices, settings) == first);
    REQUIRE(loads == 1);
    settings.interior_density = 0.1;
    REQUIRE(cache.get(slices, settings) != first);
    REQUIRE(loads == 1);
    cache.clear();
    REQUIRE(cache.get(slices, settings) != first);
    REQUIRE(loads == 2);
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


TEST_CASE("Adaptive geometry shares raw distances across modes and density settings", "[InfillDensityField][Regression]")
{
    const ExPolygons model {rectangle(-8., -8., 8., 8.)};
    const auto slices = [&] { return std::vector<WallSlice> {{0., 16., &model}}; };
    InfillDensityCache cache;
    const auto blend = cache.geometry(slices, AdaptiveInfillMode::SmoothBlend);
    const auto warp = cache.geometry(slices, AdaptiveInfillMode::DistanceWarp);
    REQUIRE(blend == cache.geometry(slices, AdaptiveInfillMode::SmoothBlend));
    REQUIRE(&blend->distances() == &warp->distances());
    REQUIRE_FALSE(blend->distances().has_signed_distances());
    InfillDensitySettings settings;
    settings.mode = AdaptiveInfillMode::SmoothBlend;
    const auto first = cache.get(slices, settings);
    settings.material_density = 0.3;
    REQUIRE(cache.get(slices, settings) != first);
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
    const InfillDensitySettings settings {mode, 0.3, 0.03, gradient};
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
    REQUIRE(cache.get(slices, settings) == material);
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
    InfillDensitySettings settings {AdaptiveInfillMode::SmoothBlend, 0.3, 0.03, AdaptiveInfillGradient::Quadratic};
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


TEST_CASE("Normal density profiles stay constant along their selected axis", "[InfillDensityField][Voronoi]")
{
    const int axis = GENERATE(0, 1, 2);
    const auto mode = axis == 0 ? AdaptiveInfillMode::NormalX : axis == 1 ? AdaptiveInfillMode::NormalY : AdaptiveInfillMode::NormalZ;
    const ExPolygons model {rectangle(-16., -10., 16., 10.)};
    InfillDensityField field({{0., 24., &model}}, {}, mode, 0.3, 0.03);
    const Vec3d center(0., 0., 12.), extent(32., 20., 24.);
    Vec3d along = center, across = center;
    along[axis] += 0.25 * extent[axis];
    across[(axis + 1) % 3] += 0.3 * extent[(axis + 1) % 3];
    REQUIRE_THAT(field.density(along), Catch::Matchers::WithinAbs(field.density(center), 0.001));
    REQUIRE(field.density(across) > field.density(center) + 0.03);
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
