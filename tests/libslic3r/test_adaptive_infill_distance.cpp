#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "catch2/matchers/catch_matchers.hpp"
#include "libslic3r/Fill/AdaptiveInfillField.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/PrintConfig.hpp"
#include <cmath>
#include <functional>
#include <memory>
#include <cstddef>
#include <stdexcept>
#include <vector>

using namespace Slic3r;

namespace {
ExPolygon rectangle(double x0, double y0, double x1, double y1)
{
    return ExPolygon {Point::new_scale(x0, y0), Point::new_scale(x1, y0),
        Point::new_scale(x1, y1), Point::new_scale(x0, y1)};
}
}

TEST_CASE("Signed voxel distances preserve the original interior transform", "[AdaptiveInfillDistance]")
{
    const int axis = GENERATE(-1, 0, 1, 2);
    const ExPolygons model {rectangle(-8., -6., 8., 6.)};
    const std::vector<AdaptiveInfillField::Slice> slices {{0., 16., &model}};
    const BoundingBox bounds = model.front().contour.bounding_box();
    const AdaptiveInfillDistanceField interior(slices, bounds, {}, {}, false, axis);
    const AdaptiveInfillDistanceField signed_field(slices, bounds, {}, {}, true, axis);
    REQUIRE_FALSE(interior.has_signed_distances());
    REQUIRE(signed_field.has_signed_distances());
    REQUIRE(interior.node_count() == signed_field.node_count());
    for (size_t i = 0; i < interior.node_count(); ++i) {
        REQUIRE(interior.is_inside(i) == signed_field.is_inside(i));
        REQUIRE(interior.squared_depth(i) == signed_field.squared_depth(i));
    }
}

TEST_CASE("Signed voxel distances include walls roofs floors and holes", "[AdaptiveInfillDistance]")
{
    ExPolygons model {rectangle(-10., -10., 10., 10.)};
    model.front().holes.push_back(rectangle(-2., -2., 2., 2.).contour);
    const AdaptiveInfillDistanceField field({{0., 20., &model}}, {}, {}, {}, true);
    REQUIRE(field.signed_distance(Vec3d(6., 0., 10.)) < 0.);
    REQUIRE(field.signed_distance(Vec3d(0., 0., 10.)) > 0.);
    REQUIRE(field.signed_distance(Vec3d(11., 0., 10.)) > 0.);
    REQUIRE(field.signed_distance(Vec3d(6., 0., -1.)) > 0.);
    REQUIRE(field.signed_distance(Vec3d(6., 0., 21.)) > 0.);
    REQUIRE_THAT(field.signed_distance(Vec3d(6., 0., 1.)), Catch::Matchers::WithinAbs(-1., field.grid_spacing()));
}

TEST_CASE("Signed voxel distances keep sloping and changing slice profiles", "[AdaptiveInfillDistance]")
{
    const ExPolygons lower {rectangle(-8., -8., 8., 8.)}, upper {rectangle(-4., -4., 4., 4.)};
    const AdaptiveInfillDistanceField field({{0., 8., &lower}, {8., 16., &upper}}, {}, {}, {}, true);
    REQUIRE(field.signed_distance(Vec3d(6., 0., 4.)) < 0.);
    REQUIRE(field.signed_distance(Vec3d(6., 0., 12.)) > 0.);
    REQUIRE(field.signed_distance(Vec3d(0., 0., 12.)) < 0.);
}

TEST_CASE("Voxel grids coarsen large objects and retain finite exterior queries", "[AdaptiveInfillDistance]")
{
    const ExPolygons model {rectangle(-1000., -1000., 1000., 1000.)};
    const AdaptiveInfillDistanceField field({{0., 2000., &model}}, {}, {}, {}, true);
    REQUIRE(field.node_count() <= size_t(1 << 20));
    REQUIRE(field.grid_spacing() > 0.5);
    REQUIRE(std::isfinite(field.signed_distance(Vec3d(3000., 0., 1000.))));
    REQUIRE(field.signed_distance(Vec3d(3000., 0., 1000.)) > 0.);
}

TEST_CASE("Voxel distance construction reports progress and supports cancellation", "[AdaptiveInfillDistance]")
{
    const ExPolygons model {rectangle(-8., -8., 8., 8.)};
    double previous = 0.;
    bool complete = false;
    const AdaptiveInfillDistanceField field({{0., 16., &model}}, {}, {}, [&](auto stage, double fraction) {
        REQUIRE(fraction >= previous);
        previous = fraction;
        complete = stage == AdaptiveInfillDistanceField::Stage::Complete;
    }, true);
    REQUIRE(complete);
    REQUIRE_THROWS_AS(AdaptiveInfillDistanceField({{0., 16., &model}}, {}, [] { throw std::runtime_error("cancel"); }, {}, true), std::runtime_error);
}

TEST_CASE("Signed distances do not change adaptive depth or radial profiles", "[AdaptiveInfillDistance][AdaptiveInfillField]")
{
    const auto mode = GENERATE(AdaptiveInfillMode::DistanceWarp, AdaptiveInfillMode::SmoothBlend,
        AdaptiveInfillMode::SteppedShells, AdaptiveInfillMode::Lobes,
        AdaptiveInfillMode::NormalX, AdaptiveInfillMode::NormalY, AdaptiveInfillMode::NormalZ);
    const int axis = mode == AdaptiveInfillMode::NormalX ? 0 : mode == AdaptiveInfillMode::NormalY ? 1 : mode == AdaptiveInfillMode::NormalZ ? 2 : -1;
    const ExPolygons model {rectangle(-8., -6., 8., 6.)};
    const BoundingBox bounds = model.front().contour.bounding_box();
    const std::vector<AdaptiveInfillField::Slice> slices {{0., 16., &model}};
    const AdaptiveInfillField original(slices, bounds, mode, [] {});
    const auto signed_grid = std::make_shared<AdaptiveInfillDistanceField>(slices, bounds,
        std::function<void()>{}, AdaptiveInfillDistanceField::Progress{}, true, axis);
    const AdaptiveInfillField extended(slices, bounds, mode, [] {}, signed_grid);
    REQUIRE_FALSE(original.distances().has_signed_distances());
    for (const Vec3d &point : {Vec3d(0., 0., 8.), Vec3d(4., 2., 3.), Vec3d(-5., -3., 12.)}) {
        if (mode == AdaptiveInfillMode::DistanceWarp || mode == AdaptiveInfillMode::SmoothBlend || mode == AdaptiveInfillMode::SteppedShells) {
            REQUIRE_THAT(original.depth(point), Catch::Matchers::WithinAbs(extended.depth(point), 1e-12));
            continue;
        }
        AdaptiveInfillField::Radials a, b;
        const auto count = original.radial(point, a);
        REQUIRE(count == extended.radial(point, b));
        for (size_t i = 0; i < count; ++i) {
            REQUIRE_THAT(a[i].t, Catch::Matchers::WithinAbs(b[i].t, 1e-12));
            REQUIRE_THAT(a[i].weight, Catch::Matchers::WithinAbs(b[i].weight, 1e-12));
            REQUIRE_THAT((a[i].center - b[i].center).norm(), Catch::Matchers::WithinAbs(0., 1e-12));
        }
    }
}

TEST_CASE("Empty adaptive geometry exposes an empty distance grid", "[AdaptiveInfillDistance]")
{
    const AdaptiveInfillField field({}, {}, AdaptiveInfillMode::SmoothBlend, [] {});
    REQUIRE(field.empty());
    REQUIRE(field.distances().empty());
    REQUIRE(field.distances().sample(std::vector<Vec3d>{}).empty());
}
