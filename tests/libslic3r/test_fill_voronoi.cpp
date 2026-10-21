#include <algorithm>
#include <catch2/catch_all.hpp>
#include <math.h>
#include <tbb/parallel_for.h>

#include <Eigen/Core>
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/InfillDensityField.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polyline.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "catch2/catch_test_macros.hpp"
#include "libslic3r/Fill/VoronoiPointCloud.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "catch2/generators/catch_generators.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"
#include "libslic3r/libslic3r.h"
#include "libslic3r/Fill/VoronoiInfillCgal.hpp"
#include "libslic3r/VoronoiCloudMethod.hpp"
#include "libslic3r/Fill/VoronoiInfill.hpp"
#include "libslic3r/Fill/VoronoiRouting.hpp"
#include "libslic3r/Fill/VoronoiWallDistance.hpp"
#include "libslic3r/Fill/AdaptiveInfillField.hpp"
#include "libslic3r/InfillMaterialDensity.hpp"
#include "libslic3r/Fill/VoronoiLattice.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/EdgeGrid.hpp"
#include <atomic>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <memory>
#include <set>
#include <utility>
#include <vector>
#include <stdexcept>
#include <cstdint>
#include <future>

using namespace Slic3r;

namespace {

ExPolygon rectangle(double x0, double y0, double x1, double y1)
{
    return ExPolygon{Point::new_scale(x0, y0), Point::new_scale(x1, y0),
                     Point::new_scale(x1, y1), Point::new_scale(x0, y1)};
}

double path_length(const Polylines &paths)
{
    double length = 0.;
    for (const Polyline &path : paths)
        for (size_t i = 1; i < path.points.size(); ++i)
            length += (path.points[i].cast<double>() - path.points[i - 1].cast<double>()).norm();
    return length;
}

using Segments = std::set<std::pair<Point, Point>>;

Segments segments(const Polylines &paths)
{
    Segments result;
    for (const Polyline &path : paths) {
        REQUIRE(path.points.size() >= 2);
        for (size_t i = 1; i < path.points.size(); ++i) {
            Point a = path.points[i - 1], b = path.points[i];
            REQUIRE(a != b);
            if (b < a)
                std::swap(a, b);
            REQUIRE(result.emplace(a, b).second);
        }
    }
    return result;
}

class FixedProvider : public Voronoi::PointCloudProvider
{
public:
    double spacing(double s, double h, double d) const override
        { return Voronoi::PoissonPointCloud::nominal_spacing(s, h, d); }
    std::vector<Vec3d> sites;
    std::vector<Vec3d> points_in(const BoundingBoxf3 &bounds, double) const override
    {
        std::vector<Vec3d> result;
        for (const Vec3d &site : sites)
            if ((site.array() >= bounds.min.array()).all() && (site.array() < bounds.max.array()).all())
                result.push_back(site);
        return result;
    }
    std::optional<BoundingBoxf3> extent() const override
    {
        BoundingBoxf3 box;
        for (const Vec3d &site : sites) box.merge(site);
        return box;
    }
};

class FixedSites : public VoronoiInfill
{
public:
    FixedSites() : FixedSites(std::make_shared<FixedProvider>()) {}
    std::vector<Vec3d> &sites;
private:
    explicit FixedSites(std::shared_ptr<FixedProvider> provider) : VoronoiInfill(provider), sites(provider->sites) {}
};
using UniformSites = Voronoi::PoissonPointCloud;

} // namespace

TEST_CASE("Voronoi infill uses the requested material density", "[FillVoronoi]")
{
    const double density = GENERATE(0.05, 0.10, 0.20, 0.40);
    const double z = GENERATE(0., 13.);
    const double spacing = 0.45;
    const ExPolygon region = rectangle(0., 0., 200., 200.);
    VoronoiInfill infill;
    const Polylines paths = infill.fill_geometry(region.contour.bounding_box(), region, {spacing, density, z, spacing, 0.}).walls;
    REQUIRE_FALSE(paths.empty());
    const double actual = path_length(paths) * scale_(spacing) / region.area();
    REQUIRE_THAT(actual, Catch::Matchers::WithinRel(density, 0.10));
    segments(paths); // Every shared wall is printed once.
}

TEST_CASE("Voronoi faces move with height according to three dimensional distance", "[FillVoronoi]")
{
    const double z = GENERATE(0., 10., 20.);
    const ExPolygon region = rectangle(-20., -20., 20., 20.);
    FixedSites infill;
    infill.sites = {Vec3d(-10., 0., 0.), Vec3d(10., 0., 20.)};
    const Polylines paths = infill.fill_geometry(region.contour.bounding_box(), region, {0.45, 0.2, z, 0.45, 0.}).walls;
    REQUIRE(paths.size() == 1);
    // The 3D bisector is x+z=10. Its section moves left as the layer rises.
    REQUIRE_THAT(unscale_(paths.front().first_point().x()), Catch::Matchers::WithinAbs(10. - z, EPSILON));
    REQUIRE_THAT(unscale_(paths.front().last_point().x()), Catch::Matchers::WithinAbs(10. - z, EPSILON));
}

TEST_CASE("Voronoi cells appear and disappear as the slice passes their height", "[FillVoronoi]")
{
    const ExPolygon region = rectangle(-20., -20., 20., 20.);
    FixedSites infill;
    infill.sites = {Vec3d(-10., 0., 0.), Vec3d(10., 0., 0.), Vec3d(0., 0., 20.)};
    const BoundingBox bounds = region.contour.bounding_box();
    REQUIRE(infill.fill_geometry(bounds, region, {0.45, 0.2, 0., 0.45, 0.}).walls.size() == 1);
    const Polylines middle = infill.fill_geometry(bounds, region, {0.45, 0.2, 10., 0.45, 0.}).walls;
    REQUIRE(middle.size() == 2);
    for (const Polyline &path : middle)
        REQUIRE_THAT(std::abs(unscale_(path.first_point().x())), Catch::Matchers::WithinAbs(5., EPSILON));
    REQUIRE(infill.fill_geometry(bounds, region, {0.45, 0.2, 20., 0.45, 0.}).walls.empty());
}

TEST_CASE("Voronoi slice walls separate the nearest three dimensional sites", "[FillVoronoi]")
{
    const double z = GENERATE(-3., 4., 12., 20.);
    const ExPolygon region = rectangle(-20., -20., 20., 20.);
    FixedSites infill;
    infill.sites = {Vec3d(-12., -8., -6.), Vec3d(12., -8., 2.), Vec3d(0., 14., 10.), Vec3d(-2., -1., 25.)};
    const Polylines paths = infill.fill_geometry(region.contour.bounding_box(), region, {0.45, 0.2, z, 0.45, 0.}).walls;
    REQUIRE_FALSE(paths.empty());
    segments(paths);
    for (const Polyline &path : paths) {
        for (size_t i = 1; i < path.points.size(); ++i) {
            const Vec2d xy = (path.points[i - 1].cast<double>() + path.points[i].cast<double>()) * (0.5 * SCALING_FACTOR);
            const Vec3d midpoint(xy.x(), xy.y(), z);
            std::vector<double> distances;
            for (const Vec3d &site : infill.sites)
                distances.push_back((site - midpoint).squaredNorm());
            std::sort(distances.begin(), distances.end());
            REQUIRE_THAT(distances[1], Catch::Matchers::WithinAbs(distances[0], EPSILON));
        }
    }
}

TEST_CASE("Voronoi layers change while a moving point band remains reproducible", "[FillVoronoi]")
{
    class CountingProvider : public Voronoi::PoissonPointCloud {
    public:
        mutable std::vector<BoundingBoxf3> queries;
        std::vector<Vec3d> points_in_cells(const Vec3i64 &first, const Vec3i64 &last, double spacing) const override
        {
            const BoundingBoxf3 bounds(first.cast<double>() * cell_size(spacing), last.cast<double>() * cell_size(spacing));
            queries.push_back(bounds);
            return PoissonPointCloud::points_in_cells(first, last, spacing);
        }
    };
    const ExPolygon region = rectangle(-50., -50., 50., 50.);
    const BoundingBox bounds = region.contour.bounding_box();
    auto provider = std::make_shared<CountingProvider>();
    VoronoiInfill infill(provider);
    const Segments first = segments(infill.fill_geometry(bounds, region, {0.45, 0.15, 0.1, 0.45, 0.}).walls);
    REQUIRE_FALSE(provider->queries.empty());
    REQUIRE(first != segments(infill.fill_geometry(bounds, region, {0.45, 0.15, 0.2, 0.45, 0.}).walls));
    infill.fill_geometry(bounds, region, {0.45, 0.15, 1000., 0.45, 0.});
    // Jumping past the active range starts locally, without ingesting the gap.
    REQUIRE(provider->queries.back().min.z() > 900.);
    REQUIRE(first == segments(infill.fill_geometry(bounds, region, {0.45, 0.15, 0.1, 0.45, 0.}).walls));
}

TEST_CASE("Voronoi infill repeats at the same height across independently filled regions", "[FillVoronoi]")
{
    const ExPolygon whole = rectangle(-50., -50., 50., 50.);
    const ExPolygon part = rectangle(-20., -10., 30., 40.);
    const BoundingBox bounds = whole.contour.bounding_box();
    VoronoiInfill first, second;
    const Polylines full = first.fill_geometry(bounds, whole, {0.45, 0.15, 0., 0.45, 0.}).walls;
    const Polylines partial = second.fill_geometry(bounds, part, {0.45, 0.15, 0., 0.45, 0.}).walls;
    REQUIRE_FALSE(full.empty());
    REQUIRE(segments(full) == segments(second.fill_geometry(bounds, whole, {0.45, 0.15, 0., 0.45, 0.}).walls));
    REQUIRE(segments(intersection_pl(full, part)) == segments(partial));
    // Changing a cached density and changing it back must not change the pattern.
    first.fill_geometry(bounds, whole, {0.45, 0.30, 0., 0.45, 0.});
    REQUIRE(segments(full) == segments(first.fill_geometry(bounds, whole, {0.45, 0.15, 0., 0.45, 0.}).walls));
}

TEST_CASE("Voronoi infill stays inside concave surfaces and outside holes", "[FillVoronoi]")
{
    ExPolygon region{Point::new_scale(0., 0.), Point::new_scale(80., 0.),
                     Point::new_scale(80., 30.), Point::new_scale(30., 30.),
                     Point::new_scale(30., 80.), Point::new_scale(0., 80.)};
    Polygon hole = rectangle(5., 5., 20., 20.).contour;
    hole.reverse();
    region.holes.push_back(std::move(hole));
    VoronoiInfill infill;
    const Polylines paths = infill.fill_geometry(region.contour.bounding_box(), region, {0.45, 0.20, 0., 0.45, 0.}).walls;
    REQUIRE_FALSE(paths.empty());
    REQUIRE(path_length(diff_pl(paths, region)) < SCALED_EPSILON);
    segments(paths);
}

TEST_CASE("Voronoi infill reports invalid parameters and accepts empty requests", "[FillVoronoi]")
{
    const ExPolygon region = rectangle(0., 0., 10., 10.);
    VoronoiInfill infill;
    const double density = GENERATE(-0.1, std::numeric_limits<double>::denorm_min(),
        std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN());
    REQUIRE_THROWS_AS(infill.fill_geometry(region.contour.bounding_box(), region, {0.45, density, 0., 0.45, 0.}), std::invalid_argument);
    REQUIRE_THROWS_AS(infill.fill_geometry(region.contour.bounding_box(), region, {0., 0.2, 0., 0., 0.}), std::invalid_argument);
    REQUIRE(infill.fill_geometry(region.contour.bounding_box(), region, {0.45, 0., 0., 0.45, 0.}).walls.empty());
    REQUIRE(infill.fill_geometry(region.contour.bounding_box(), ExPolygon{}, {0.45, 0.2, 0., 0.45, 0.}).walls.empty());
}

TEST_CASE("Voronoi bisectors reach the boundary for collinear sites", "[FillVoronoi]")
{
    const size_t count = GENERATE(size_t(2), size_t(3));
    const ExPolygon region = rectangle(-20., -20., 20., 20.);
    FixedSites infill;
    infill.sites = {Vec3d(-10., 0., 0.), Vec3d(10., 0., 0.)};
    if (count == 3)
        infill.sites.push_back(Vec3d(0., 0., 0.));
    const Polylines paths = infill.fill_geometry(region.contour.bounding_box(), region, {0.45, 0.2, 0., 0.45, 0.}).walls;
    REQUIRE(paths.size() == count - 1);
    REQUIRE_THAT(path_length(paths), Catch::Matchers::WithinAbs(scale_(40. * (count - 1)), SCALED_EPSILON));
    for (const Polyline &path : paths) {
        REQUIRE(path.first_point().x() == path.last_point().x());
        REQUIRE(std::abs(path.first_point().y()) == scale_(20.));
        REQUIRE(std::abs(path.last_point().y()) == scale_(20.));
    }
}

TEST_CASE("Voronoi rays and cocircular sites produce unique connected walls", "[FillVoronoi]")
{
    const size_t count = GENERATE(size_t(3), size_t(4));
    const ExPolygon region = rectangle(-20., -20., 20., 20.);
    FixedSites infill;
    infill.sites = {Vec3d(-10., -10., 0.), Vec3d(10., -10., 0.), Vec3d(10., 10., 0.)};
    if (count == 4)
        infill.sites.push_back(Vec3d(-10., 10., 0.));
    const Polylines paths = infill.fill_geometry(region.contour.bounding_box(), region, {0.45, 0.2, 0., 0.45, 0.}).walls;
    REQUIRE(segments(paths).size() == count);
    REQUIRE(path_length(diff_pl(paths, region)) < SCALED_EPSILON);
    // Each wall separates its two nearest sites, including the infinite rays.
    for (const Polyline &path : paths) {
        for (size_t i = 1; i < path.points.size(); ++i) {
            const Vec2d xy = (path.points[i - 1].cast<double>() + path.points[i].cast<double>()) * (0.5 * SCALING_FACTOR);
            const Vec3d midpoint(xy.x(), xy.y(), 0.);
            std::vector<double> distances;
            for (const Vec3d &site : infill.sites)
                distances.push_back((site - midpoint).squaredNorm());
            std::sort(distances.begin(), distances.end());
            REQUIRE_THAT(distances[1], Catch::Matchers::WithinAbs(distances[0], EPSILON));
        }
    }
}

TEST_CASE("Voronoi infill accepts duplicate sites and large translated coordinates", "[FillVoronoi]")
{
    const double translation = GENERATE(-5000., 0., 5000.);
    const double separation = GENERATE(10., 10000.);
    const ExPolygon region = rectangle(translation - 20., -20., translation + 20., 20.);
    FixedSites infill;
    infill.sites = {Vec3d(translation - separation, 0., 0.), Vec3d(translation + separation, 0., 0.),
                   Vec3d(translation - separation, 0., 0.)};
    const Polylines paths = infill.fill_geometry(region.contour.bounding_box(), region, {0.45, 0.2, 0., 0.45, 0.}).walls;
    REQUIRE(paths.size() == 1);
    REQUIRE(paths.front().first_point().x() == scale_(translation));
    REQUIRE(paths.front().last_point().x() == scale_(translation));
    REQUIRE_THAT(path_length(paths), Catch::Matchers::WithinAbs(scale_(40.), SCALED_EPSILON));
}

TEST_CASE("Voronoi infill produces no walls with fewer than two distinct sites", "[FillVoronoi]")
{
    const size_t count = GENERATE(size_t(0), size_t(1), size_t(2));
    const ExPolygon region = rectangle(-20., -20., 20., 20.);
    FixedSites infill;
    infill.sites.assign(count, Vec3d(0., 0., 0.));
    REQUIRE(infill.fill_geometry(region.contour.bounding_box(), region, {0.45, 0.2, 0., 0.45, 0.}).walls.empty());
}

TEST_CASE("Voronoi sites are uniformly distributed in three dimensions", "[FillVoronoi]")
{
    const BoundingBoxf3 bounds(Vec3d(-100., -100., -100.), Vec3d(100., 100., 100.));
    UniformSites sampler;
    const std::vector<Vec3d> sites = sampler.points_in(bounds, 10.);
    REQUIRE_THAT(double(sites.size()), Catch::Matchers::WithinRel(8000., 0.05));
    std::array<size_t, 64> bins{};
    for (const Vec3d &p : sites) {
        REQUIRE(bounds.contains(p));
        const size_t x = size_t((p.x() + 100.) / 50.);
        const size_t y = size_t((p.y() + 100.) / 50.);
        const size_t z = size_t((p.z() + 100.) / 50.);
        ++bins[z * 16 + y * 4 + x];
    }
    double chi_squared = 0.;
    const double expected = double(sites.size()) / bins.size();
    for (size_t bin : bins)
        chi_squared += std::pow(double(bin) - expected, 2.) / expected;
    REQUIRE(chi_squared < 110.); // 64 uniform bins, 63 degrees of freedom.
}

TEST_CASE("Voronoi shallow face bands overlap adjacent layers by half", "[FillVoronoi]")
{
    const ExPolygon region = rectangle(-20., -20., 20., 20.);
    const BoundingBox bounds = region.contour.bounding_box();
    FixedSites infill;
    // Bisector x+5*z=25 moves one millimeter per 0.2 mm layer.
    infill.sites = {Vec3d(-1., 0., 0.), Vec3d(1., 0., 10.)};
    const auto below = infill.fill_geometry(bounds, region, {0.4, 0.2, 4.8, 0.45, 0.2});
    const auto current = infill.fill_geometry(bounds, region, {0.4, 0.2, 5., 0.45, 0.2});
    const auto above = infill.fill_geometry(bounds, region, {0.4, 0.2, 5.2, 0.45, 0.2});
    REQUIRE(path_length(current.walls) < SCALED_EPSILON);
    REQUIRE(current.bands.size() == 1);
    REQUIRE(below.bands.size() == 1);
    REQUIRE(above.bands.size() == 1);
    const ExPolygons polygon {current.bands.front().polygon};
    const double area = current.bands.front().polygon.area();
    for (const auto *neighbor : {&below, &above}) {
        const ExPolygons overlap = intersection_ex(polygon, ExPolygons{neighbor->bands.front().polygon});
        REQUIRE(overlap.size() == 1);
        REQUIRE_THAT(overlap.front().area() / area, Catch::Matchers::WithinAbs(0.5, 1e-6));
    }
    REQUIRE(current.bands.front().supported_side.x() > 0.);
}

TEST_CASE("Voronoi rounds shallow floors outward to whole rows with at least half layer overlap", "[FillVoronoi]")
{
    const double shift = GENERATE(0.23, 0.3, 0.5, 1.07);
    const double offset = GENERATE(0., 37., -63.);
    const double height = 0.2, spacing = 0.4070796327;
    const ExPolygon region = rectangle(offset - 20., -20., offset + 20., 20.);
    FixedSites infill;
    infill.sites = {Vec3d(offset - 1., 0., 5. - shift / height), Vec3d(offset + 1., 0., 5. + shift / height)};
    const auto bounds = region.contour.bounding_box();
    const auto below = infill.fill_geometry(bounds, region, {spacing, 0.2, 5. - height, 0.45, height});
    const auto current = infill.fill_geometry(bounds, region, {spacing, 0.2, 5., 0.45, height});
    const auto above = infill.fill_geometry(bounds, region, {spacing, 0.2, 5. + height, 0.45, height});
    REQUIRE(current.bands.size() == 1);
    const auto &band = current.bands.front();
    const double rounded_width = std::ceil(2. * shift / spacing) * spacing;
    REQUIRE_THAT(unscale_(band.polygon.contour.bounding_box().max.x()), Catch::Matchers::WithinAbs(offset + shift, 1e-5));
    REQUIRE_THAT(unscale_(band.polygon.contour.bounding_box().min.x()),
                 Catch::Matchers::WithinAbs(offset + shift - rounded_width, 1e-5));
    REQUIRE_THAT(band.row_origin, Catch::Matchers::WithinAbs(offset + shift, 1e-9));
    for (const auto *neighbor : {&below, &above}) {
        REQUIRE(neighbor->bands.size() == 1);
        const ExPolygons overlap = intersection_ex(ExPolygons{band.polygon}, ExPolygons{neighbor->bands.front().polygon});
        REQUIRE(overlap.size() == 1);
        REQUIRE(overlap.front().area() / band.polygon.area() >= 0.5 - 1e-6);
    }
}

TEST_CASE("Voronoi faces retain a single line when bead overlap is sufficient", "[FillVoronoi]")
{
    const double width = GENERATE(0.45, 0.8);
    const ExPolygon region = rectangle(-20., -20., 20., 20.);
    FixedSites infill;
    infill.sites = {Vec3d(-10., 0., 0.), Vec3d(10., 0., 20.)};
    const auto layer = infill.fill_geometry(region.contour.bounding_box(), region, {0.4, 0.2, 10., width, 0.2});
    REQUIRE(layer.bands.empty());
    REQUIRE(layer.walls.size() == 1);
    // At the same height a narrower bead needs a band instead; cached geometry
    // must include extrusion dimensions, not just the layer Z coordinate.
    const auto narrow = infill.fill_geometry(region.contour.bounding_box(), region, {0.2, 0.2, 10., 0.3, 0.2});
    REQUIRE(path_length(narrow.walls) < SCALED_EPSILON);
    REQUIRE(narrow.bands.size() == 1);
}

TEST_CASE("Voronoi horizontal faces fill once and stop at neighboring cells", "[FillVoronoi]")
{
    const ExPolygon region = rectangle(-20., -20., 20., 20.);
    const BoundingBox bounds = region.contour.bounding_box();
    FixedSites infill;
    infill.sites = {Vec3d(0., 0., 0.), Vec3d(0., 0., 1.), Vec3d(4., 0., 0.5),
                    Vec3d(-4., 0., 0.5), Vec3d(0., 4., 0.5), Vec3d(0., -4., 0.5)};
    REQUIRE(infill.fill_geometry(bounds, region, {0.4, 0.2, 0.25, 0.45, 0.25}).bands.empty());
    REQUIRE(infill.fill_geometry(bounds, region, {0.4, 0.2, 0.75, 0.45, 0.25}).bands.empty());
    const auto floor = infill.fill_geometry(bounds, region, {0.4, 0.2, 0.5, 0.45, 0.25});
    REQUIRE(floor.bands.size() == 1);
    REQUIRE_THAT(floor.bands.front().polygon.area() / (scale_(1.) * scale_(1.)),
                 Catch::Matchers::WithinAbs(3.9375 * 3.9375, 1e-5));
}

TEST_CASE("Voronoi bands respect holes and do not duplicate walls or polygon areas", "[FillVoronoi]")
{
    ExPolygon region = rectangle(-20., -20., 20., 20.);
    Polygon hole = rectangle(-2., -2., 2., 2.).contour;
    hole.reverse();
    region.holes.push_back(hole);
    FixedSites infill;
    infill.sites = {Vec3d(-1., 0., -10.), Vec3d(1., 0., 10.), Vec3d(0., -1., -10.), Vec3d(0., 1., 10.)};
    const auto layer = infill.fill_geometry(region.contour.bounding_box(), region, {0.4, 0.2, 0., 0.45, 0.2});
    REQUIRE_FALSE(layer.bands.empty());
    ExPolygons covered;
    for (const auto &band : layer.bands) {
        REQUIRE(diff_ex(ExPolygons{band.polygon}, ExPolygons{region}).empty());
        REQUIRE(intersection_ex(ExPolygons{band.polygon}, covered).empty());
        covered.push_back(band.polygon);
    }
    // Supporting walls on floor boundaries remain printable. Only band
    // interiors beyond the fixed-point cleanup tolerance exclude walls.
    ExPolygons interiors;
    for (const auto &band : layer.bands)
        append(interiors, offset_ex(band.polygon, -float(scale_(0.00001))));
    REQUIRE(intersection_pl(layer.walls, interiors).empty());
}

TEST_CASE("Voronoi shallow face compensation preserves statistical material density", "[FillVoronoi]")
{
    const double height = GENERATE(0.12, 0.2, 0.3);
    const double density = GENERATE(0.1, 0.2);
    const ExPolygon region = rectangle(0., 0., 100., 100.);
    VoronoiInfill infill;
    const auto layer = infill.fill_geometry(region.contour.bounding_box(), region, {0.4, density, 5., 0.45, height});
    double material = path_length(layer.walls) * scale_(0.4);
    for (const auto &band : layer.bands)
        material += band.polygon.area();
    REQUIRE_FALSE(layer.bands.empty());
    REQUIRE_THAT(material / region.area(), Catch::Matchers::WithinRel(density, 0.15));
}

TEST_CASE("Voronoi site placement stays fixed when extrusion dimensions change", "[FillVoronoi]")
{
    class RecordingProvider : public Voronoi::PoissonPointCloud {
    public:
        mutable std::vector<double> spacings;
        std::vector<Vec3d> points_in_cells(const Vec3i64 &first, const Vec3i64 &last, double spacing) const override
        {
            spacings.push_back(spacing);
            return PoissonPointCloud::points_in_cells(first, last, spacing);
        }
    };
    const ExPolygon region = rectangle(-20., -20., 20., 20.);
    auto provider = std::make_shared<RecordingProvider>();
    VoronoiInfill infill(provider);
    const BoundingBox bounds = region.contour.bounding_box();
    const double nominal = Voronoi::PoissonPointCloud::nominal_spacing(0.4, 0.2, 0.2);
    infill.fill_geometry(bounds, region, {0.4, 0.2, 0.3, 0.6, 0.3, 0.2, -1., 1, nominal});
    REQUIRE_FALSE(provider->spacings.empty());
    const double initial_spacing = provider->spacings.front();
    infill.fill_geometry(bounds, region, {0.47, 0.2, 20., 0.45, 0.12, 0.2, -1., 1, nominal});
    for (double spacing : provider->spacings)
        REQUIRE_THAT(spacing, Catch::Matchers::WithinAbs(initial_spacing, EPSILON));
}

TEST_CASE("Voronoi overlapping bands have reproducible ownership across independent fills", "[FillVoronoi]")
{
    const ExPolygon region = rectangle(-20., -20., 20., 20.);
    const BoundingBox bounds = region.contour.bounding_box();
    VoronoiInfill first, second, unrelated;
    const auto expected = first.fill_geometry(bounds, region, {0.4, 0.2, 5., 0.45, 0.2});
    unrelated.fill_geometry(bounds, region, {0.4, 0.3, 13., 0.45, 0.3});
    const auto actual = second.fill_geometry(bounds, region, {0.4, 0.2, 5., 0.45, 0.2});
    REQUIRE_FALSE(expected.bands.empty());
    REQUIRE(actual.bands.size() == expected.bands.size());
    for (size_t i = 0; i < expected.bands.size(); ++i) {
        REQUIRE_THAT((expected.bands[i].supported_side - actual.bands[i].supported_side).norm(),
                     Catch::Matchers::WithinAbs(0., EPSILON));
        ExPolygons difference = diff_ex(ExPolygons{expected.bands[i].polygon}, ExPolygons{actual.bands[i].polygon});
        append(difference, diff_ex(ExPolygons{actual.bands[i].polygon}, ExPolygons{expected.bands[i].polygon}));
        double area = 0.;
        for (const auto &polygon : difference)
            area += polygon.area();
        REQUIRE(area < scale_(EPSILON) * scale_(1.));
    }
}

TEST_CASE("Voronoi snaps short edges without transitive vertex drift", "[FillVoronoi]")
{
    const Polylines walls {Polyline(Points{Point(0, 0), Point(6, 0), Point(12, 0), Point(0, 100)})};
    const auto graph = Voronoi::build_wall_graph(walls, 8);
    REQUIRE(graph.indices.at(Point(6, 0)) == graph.indices.at(Point(0, 0)));
    REQUIRE(graph.indices.at(Point(12, 0)) != graph.indices.at(Point(0, 0)));
    REQUIRE(graph.edges.size() == 2);
    for (const auto &entry : graph.indices)
        REQUIRE((entry.first.cast<double>() - graph.vertices[entry.second].point.cast<double>()).norm() <= 8.);
}

TEST_CASE("Voronoi prints straight through rounded four way and six way crossings", "[FillVoronoi]")
{
    const size_t degree = GENERATE(size_t(4), size_t(6));
    const bool rounded = GENERATE(false, true);
    Polylines walls;
    for (size_t i = 0; i < degree; ++i) {
        const double angle = 2. * M_PI * double(i) / double(degree);
        const Point center = rounded && i % 2 ? Point(3, 0) : Point(0, 0);
        walls.emplace_back(Points{center, Point::new_scale(20. * std::cos(angle), 20. * std::sin(angle))});
    }
    if (rounded) walls.emplace_back(Points{Point(0, 0), Point(3, 0)});
    const auto graph = Voronoi::build_wall_graph(walls);
    const auto connections = Voronoi::pair_junctions(graph, {}, 0);
    const auto trails = Voronoi::trace_trails(graph, connections, 0);
    REQUIRE(trails.size() == degree / 2);
    REQUIRE(segments(trails).size() == degree);
    for (const Polyline &trail : trails) {
        REQUIRE(trail.points.size() == 3);
        const Vec2d a = (trail.points[0] - trail.points[1]).cast<double>();
        const Vec2d b = (trail.points[2] - trail.points[1]).cast<double>();
        REQUIRE(a.dot(b) / (a.norm() * b.norm()) < -0.99999);
    }
}

TEST_CASE("Voronoi cycles actual three way connections using persistent generating sites", "[FillVoronoi]")
{
    const BoundingBox bounds(Point::new_scale(-20., -20.), Point::new_scale(20., 20.));
    std::vector<Vec3d> sites {Vec3d(0., 2., 0.), Vec3d(-std::sqrt(3.), -1., 0.), Vec3d(std::sqrt(3.), -1., 0.)};
    const bool rounded_metadata = GENERATE(false, true);
    const bool duplicate_site = GENERATE(false, true);
    if (duplicate_site) sites.push_back(sites.front());
    std::set<Point> unpaired;
    Segments expected;
    for (size_t layer_id = 0; layer_id < 3; ++layer_id) {
        // Reordered point insertion and a changing Z plane must not change
        // which physical branch occupies each position in the three-layer cycle.
        std::rotate(sites.begin(), sites.begin() + 1, sites.end());
        auto section = Voronoi::section(sites, bounds, .2 * layer_id);
        const Point junction_point = section.junctions.front().point;
        if (rounded_metadata) section.junctions.front().point += Point(1, 0);
        const auto graph = Voronoi::build_wall_graph(section.walls);
        const auto connections = Voronoi::pair_junctions(graph, section.junctions, layer_id);
        auto trails = Voronoi::trace_trails(graph, connections, layer_id);
        REQUIRE(trails.size() == 2);
        if (layer_id == 0) expected = segments(trails);
        REQUIRE(segments(trails) == expected);
        REQUIRE(section.junctions.size() == 1);
        const size_t center = graph.indices.at(junction_point);
        size_t open = 0;
        for (size_t e : graph.vertices[center].edges) {
            const size_t side = graph.edges[e][0] == center ? 0 : 1;
            if (connections[e][side] == Voronoi::no_connection) {
                ++open;
                unpaired.insert(graph.vertices[graph.edges[e][1 - side]].point);
            }
        }
        REQUIRE(open == 1);
        // Test actual traced paths, not only the intended pairing table.
        size_t endpoints_at_center = 0;
        for (const auto &trail : trails) {
            endpoints_at_center += trail.first_point() == graph.vertices[center].point;
            endpoints_at_center += trail.last_point() == graph.vertices[center].point;
        }
        REQUIRE(endpoints_at_center == 1);
    }
    REQUIRE(unpaired.size() == 3);
}

TEST_CASE("Voronoi tracing keeps parallel edges and changes closed trail starts", "[FillVoronoi]")
{
    const Polylines walls {Polyline(Points{Point(0, 0), Point(100, 0)}),
                          Polyline(Points{Point(0, 0), Point(100, 0)})};
    const auto graph = Voronoi::build_wall_graph(walls);
    const auto connections = Voronoi::pair_junctions(graph, {}, 0);
    const auto first = Voronoi::trace_trails(graph, connections, 0);
    const auto second = Voronoi::trace_trails(graph, connections, 1);
    REQUIRE(first.size() == 1);
    REQUIRE(second.size() == 1);
    REQUIRE(first.front().points.size() == 3);
    REQUIRE(first.front().first_point() != second.front().first_point());
    REQUIRE(first.front().first_point() == first.front().last_point());
    REQUIRE_THAT(path_length(first), Catch::Matchers::WithinAbs(200., 1e-9));
    REQUIRE_THAT(path_length(second), Catch::Matchers::WithinAbs(200., 1e-9));
}

TEST_CASE("Voronoi orders whole floor groups without reversing their supported sweep", "[FillVoronoi]")
{
    const std::vector<Polylines> groups {
        {Polyline(Point(100, 0), Point(1, 0))},
        {Polyline(Point(10, 0), Point(20, 0)), Polyline(Point(21, 0), Point(90, 0))},
        {Polyline(Point(200, 0), Point(210, 0))}
    };
    std::vector<Voronoi::PathGroup> input {{{}, false}};
    for (const auto &group : groups) input.push_back({group, false});
    const auto ordered = Voronoi::order_groups(std::move(input), Point(0, 0));
    REQUIRE(ordered.size() == 3);
    std::set<size_t> seen;
    for (const auto &group : ordered) {
        const auto found = std::find_if(groups.begin(), groups.end(), [&](const auto &expected) {
            return expected.front().first_point() == group.paths.front().first_point();
        });
        REQUIRE(found != groups.end());
        REQUIRE(seen.insert(size_t(found - groups.begin())).second);
        REQUIRE(group.paths.size() == found->size());
        for (size_t j = 0; j < group.paths.size(); ++j)
            REQUIRE(group.paths[j].points == (*found)[j].points);
    }
    REQUIRE(Voronoi::order_groups({}, Point(0, 0)).empty());
}

TEST_CASE("Voronoi mixes floor groups and reversible wall trails using their actual exits", "[FillVoronoi]")
{
    const Point floor_exit = GENERATE(Point(90, 10), Point(10, 10));
    const Polylines floor {Polyline(Point(10, 1), Point(80, 0)), Polyline(Point(80, 1), floor_exit)};
    const Polyline wall(Point(0, 0), Point(10, 0));
    const Polyline next_wall(floor_exit + Point(20, 0), floor_exit);
    std::vector<Voronoi::PathGroup> input {{{next_wall}, true}, {floor, false}, {{wall}, true}};
    const auto ordered = Voronoi::order_groups(input, Point(0, 0));
    REQUIRE(ordered.size() == 3);
    REQUIRE(ordered[0].reversible);
    REQUIRE_FALSE(ordered[1].reversible);
    REQUIRE(ordered[2].reversible);
    REQUIRE(ordered[0].paths.front().points == wall.points);
    REQUIRE(ordered[1].paths.size() == floor.size());
    for (size_t i = 0; i < floor.size(); ++i)
        REQUIRE(ordered[1].paths[i].points == floor[i].points);
    REQUIRE(ordered[2].paths.front().first_point() == floor_exit);
    REQUIRE(ordered[2].paths.front().last_point() == next_wall.first_point());
    // A floor can also be the first group; walls have no blanket priority.
    const Polyline distant_wall(Point(1000, 0), Point(1100, 0));
    const auto floor_first = Voronoi::order_groups({{floor, false}, {{distant_wall}, true}}, floor.front().first_point());
    REQUIRE_FALSE(floor_first.front().reversible);
    REQUIRE(floor_first.back().paths.front().first_point() == distant_wall.first_point());
    REQUIRE_THAT((floor_first.back().paths.front().first_point() - floor_first.front().paths.back().last_point()).cast<double>().norm(),
        Catch::Matchers::WithinAbs((distant_wall.first_point() - floor_exit).cast<double>().norm(), 1e-9));
    // A reversed proxy must still expand from the supported side.
    const auto only_floor = Voronoi::order_groups({{floor, false}}, floor_exit);
    REQUIRE(only_floor.size() == 1);
    REQUIRE(only_floor.front().paths.size() == floor.size());
    for (size_t i = 0; i < floor.size(); ++i)
        REQUIRE(only_floor.front().paths[i].points == floor[i].points);
}

TEST_CASE("Voronoi keeps supporting walls between horizontal floor faces", "[FillVoronoi]")
{
    FixedSites infill;
    for (int x = -1; x <= 4; ++x)
        for (int y = -1; y <= 4; ++y)
            for (int z = -1; z <= 2; ++z)
                infill.sites.emplace_back(20. * x, 20. * y, 20. * z);
    const ExPolygon region = rectangle(0., 0., 50., 50.);
    const auto layer = infill.fill_geometry(region.contour.bounding_box(), region, {.4, .05, 10.1, .45, .2});
    REQUIRE(layer.bands.size() > 1);
    for (int axis = 0; axis < 2; ++axis)
        for (double position : {10., 30.})
            for (double along : {1., 20., 40., 49.}) {
                bool covered = false;
                for (const Polyline &wall : layer.walls)
                    for (size_t i = 1; i < wall.points.size(); ++i) {
                        const Point &a = wall.points[i - 1], &b = wall.points[i];
                        if (std::abs(a[axis] - scale_(position)) <= 1 && std::abs(b[axis] - scale_(position)) <= 1 &&
                            scale_(along) >= std::min(a[1 - axis], b[1 - axis]) &&
                            scale_(along) <= std::max(a[1 - axis], b[1 - axis]))
                            covered = true;
                    }
                REQUIRE(covered);
            }
}

TEST_CASE("Voronoi point queries preserve the field when partitioned", "[FillVoronoi]")
{
    Voronoi::PoissonPointCloud provider;
    const BoundingBoxf3 whole(Vec3d(-13.2, -8.1, -9.4), Vec3d(14.7, 12.6, 11.3));
    auto keys = [](const std::vector<Vec3d> &points) {
        std::set<std::array<double, 3>> result;
        for (const Vec3d &p : points)
            REQUIRE(result.insert({p.x(), p.y(), p.z()}).second);
        return result;
    };
    const auto expected = keys(provider.points_in(whole, 3.));
    REQUIRE_FALSE(expected.empty());
    const int axis = GENERATE(0, 1, 2);
    BoundingBoxf3 first = whole, second = whole;
    first.max[axis] = second.min[axis] = 0.37;
    auto points = provider.points_in(first, 3.);
    auto rest = provider.points_in(second, 3.);
    points.insert(points.end(), rest.begin(), rest.end());
    REQUIRE(keys(points) == expected);
}

TEST_CASE("Voronoi coverage includes boundary corners without internal nodes", "[FillVoronoi]")
{
    const auto bounds = rectangle(-2., -3., 2., 3.).contour.bounding_box();
    const auto section = Voronoi::section({Vec3d(0., 0., 4.)}, bounds, 0.);
    REQUIRE(section.walls.empty());
    REQUIRE_THAT(section.covering_radius, Catch::Matchers::WithinAbs(std::sqrt(29.), EPSILON));
    REQUIRE_FALSE(std::isfinite(Voronoi::section({}, bounds, 0.).covering_radius));
}

TEST_CASE("Voronoi coverage bounds nearest distance throughout clipped cells", "[FillVoronoi]")
{
    const auto bounds = rectangle(-7., -5., 9., 8.).contour.bounding_box();
    const std::vector<Vec3d> sites {Vec3d(-5., -3., -2.), Vec3d(6., -2., 3.),
                                  Vec3d(-3., 6., 1.), Vec3d(5., 7., -4.)};
    const double z = GENERATE(-3., 0., 5.);
    const auto section = Voronoi::section(sites, bounds, z);
    REQUIRE_FALSE(section.walls.empty());
    for (double x = -7.; x <= 9.; x += 0.25)
        for (double y = -5.; y <= 8.; y += 0.25) {
            double nearest = std::numeric_limits<double>::infinity();
            for (const auto &site : sites)
                nearest = std::min(nearest, (Vec3d(x, y, z) - site).norm());
            REQUIRE(nearest <= section.covering_radius + EPSILON);
        }
}

TEST_CASE("Voronoi adaptive coverage crosses gaps and preserves an exterior field", "[FillVoronoi]")
{
    auto provider = std::make_shared<FixedProvider>();
    provider->sites = {Vec3d(-100., 0., -20.), Vec3d(100., 0., 20.),
                       Vec3d(0., 100., 0.), Vec3d(0., -100., 0.)};
    Voronoi::ActivePointCloud cloud(provider);
    const auto bounds = rectangle(-5., -5., 5., 5.).contour.bounding_box();
    for (double z : {0., 1., 20., 0.}) {
        Voronoi::Section actual;
        REQUIRE(cloud.update(bounds, z, 0.5, 1., [&](const auto &sites) {
            actual = Voronoi::section(sites, bounds, z);
            return actual.covering_radius;
        }));
        const auto expected = Voronoi::section(provider->sites, bounds, z);
        REQUIRE(segments(actual.walls) == segments(expected.walls));
        REQUIRE_THAT(actual.covering_radius, Catch::Matchers::WithinAbs(expected.covering_radius, EPSILON));
    }
}

TEST_CASE("Voronoi retirement preserves future slices and rounded floors", "[FillVoronoi]")
{
    auto provider = std::make_shared<FixedProvider>();
    // Structured sites include near-horizontal faces and points well below the slice.
    for (double z = -20.; z <= 40.; z += 4.)
        for (double x : {-6., 0., 6.})
            for (double y : {-6., 0., 6.})
                provider->sites.emplace_back(x + 0.03 * z, y, z);
    Voronoi::ActivePointCloud cloud(provider);
    const auto bounds = rectangle(-4., -4., 4., 4.).contour.bounding_box();
    const double width = 0.45, height = 0.2, row_spacing = 0.4;
    const double margin = height * (1. + 2. * row_spacing / width);
    for (double z : {0., 2., 4., 10., 20., 0.}) {
        REQUIRE(cloud.update(bounds, z, margin, 3., [&](const auto &sites) {
            return Voronoi::section(sites, bounds, z).covering_radius;
        }));
        REQUIRE(cloud.sites().size() < provider->sites.size());
        REQUIRE(segments(Voronoi::section(cloud.sites(), bounds, z).walls) ==
                segments(Voronoi::section(provider->sites, bounds, z).walls));
        Voronoi::CellGeometry actual_mesh;
        actual_mesh.update(cloud.sites());
        const auto actual = actual_mesh.shallow_bands(bounds, z, width, height, row_spacing);
        Voronoi::CellGeometry expected_mesh;
        expected_mesh.update(provider->sites);
        const auto expected = expected_mesh.shallow_bands(bounds, z, width, height, row_spacing);
        ExPolygons actual_polygons, expected_polygons;
        for (const auto &band : actual) actual_polygons.push_back(band.polygon);
        for (const auto &band : expected) expected_polygons.push_back(band.polygon);
        REQUIRE(diff_ex(actual_polygons, expected_polygons).empty());
        REQUIRE(diff_ex(expected_polygons, actual_polygons).empty());
        for (const auto &site : cloud.sites())
            REQUIRE(site.z() >= cloud.retired_below());
    }
}

TEST_CASE("Voronoi finite empty clouds terminate without geometry", "[FillVoronoi]")
{
    auto provider = std::make_shared<FixedProvider>();
    VoronoiInfill infill(provider);
    const auto region = rectangle(-5., -5., 5., 5.);
    REQUIRE(infill.fill_geometry(region.contour.bounding_box(), region, {0.45, 0.2, 0., 0.45, 0.}).walls.empty());
}

TEST_CASE("Voronoi wall distance favors walls and both corner signs", "[FillVoronoi]")
{
    using Cloud = Voronoi::WallDistancePointCloud;
    REQUIRE(InfillDensityField::density_ratio(-1., -1., 32.) > InfillDensityField::density_ratio(-15., -15., 32.));
    REQUIRE(InfillDensityField::density_ratio(-4., -3., 32.) > InfillDensityField::density_ratio(-4., -4., 32.));
    REQUIRE(InfillDensityField::density_ratio(-4., -5., 32.) > InfillDensityField::density_ratio(-4., -4., 32.));
    REQUIRE_THAT(InfillDensityField::density_ratio(-1000., -1000., 32.), Catch::Matchers::WithinAbs(0.01, EPSILON));
    REQUIRE_THAT(InfillDensityField::density_ratio(0., 0., 32.), Catch::Matchers::WithinAbs(1., EPSILON));
    REQUIRE_THAT(InfillDensityField::density_ratio(-6., -6., 32.), Catch::Matchers::WithinAbs(std::pow(0.01, 6. / 32.), EPSILON));
    REQUIRE(InfillDensityField::density_ratio(-50., -50., 32.) < 0.02);
    REQUIRE(InfillDensityField::density_ratio(-1000., -1000., 32.) > 0.);
    REQUIRE_THAT(InfillDensityField::density_ratio(5., 5., 32.), Catch::Matchers::WithinAbs(1., EPSILON));
}

TEST_CASE("Voronoi wall thinning preserves reproducibility and site intensity", "[FillVoronoi]")
{
    auto constant = [](const std::vector<Vec3d> &points) { return std::vector<InfillDensitySample>(points.size(), {-6., InfillDensityField::density_ratio(-6., -6., 32.)}); };
    Voronoi::PoissonPointCloud uniform;
    Voronoi::WallDistancePointCloud thinned(constant);
    const BoundingBoxf3 box(Vec3d(-20., -20., -20.), Vec3d(20., 20., 20.));
    auto keys = [](const std::vector<Vec3d> &points) {
        std::set<std::array<double, 3>> result;
        for (const Vec3d &p : points) REQUIRE(result.insert({p.x(), p.y(), p.z()}).second);
        return result;
    };
    const double maximum_ratio = InfillDensityField::density_ratio(0., 2., 32.);
    const auto source = keys(uniform.points_in(box, 1. / maximum_ratio));
    const auto expected = keys(thinned.points_in(box, 1.));
    REQUIRE_FALSE(expected.empty());
    REQUIRE_THAT(double(expected.size()) / source.size(), Catch::Matchers::WithinRel(std::pow(InfillDensityField::density_ratio(-6., -6., 32.) / maximum_ratio, 3), 0.15));
    for (const auto &point : expected) REQUIRE(source.count(point) == 1);
    auto first = box, second = box;
    first.max.z() = second.min.z() = 0.37;
    auto points = thinned.points_in(second, 1.);
    auto rest = thinned.points_in(first, 1.);
    points.insert(points.end(), rest.begin(), rest.end());
    REQUIRE(keys(points) == expected);
}

TEST_CASE("Voronoi wall cloud makes smaller cells near model walls", "[FillVoronoi]")
{
    const ExPolygons outlines {rectangle(-40., -40., 40., 40.)};
    auto field = std::make_shared<InfillDensityField>(std::vector<WallSlice>{{-100., 100., &outlines}});
    auto provider = std::make_shared<Voronoi::WallDistancePointCloud>([field](const auto &points) { return field->sample_fields(points); });
    VoronoiInfill infill(provider);
    const auto bounds = outlines.front().contour.bounding_box();
    const auto walls = infill.fill_geometry(bounds, outlines.front(), {0.45, 0.25, 0., 0.45, 0.}).walls;
    REQUIRE_FALSE(walls.empty());
    const auto near_wall = rectangle(30., -30., 38., 30.);
    const auto center = rectangle(-4., -30., 4., 30.);
    REQUIRE(path_length(intersection_pl(walls, near_wall)) > path_length(intersection_pl(walls, center)));
    REQUIRE(segments(infill.fill_geometry(bounds, outlines.front(), {0.45, 0.25, 5., 0.45, 0.}).walls) != segments(walls));
    REQUIRE(segments(infill.fill_geometry(bounds, outlines.front(), {0.45, 0.25, 0., 0.45, 0.}).walls) == segments(walls));
}

TEST_CASE("Voronoi wall clouds become denser near both floors and roofs", "[FillVoronoi]")
{
    const ExPolygons outlines {rectangle(-50., -50., 50., 50.)};
    auto field = std::make_shared<InfillDensityField>(std::vector<WallSlice>{{0., 50., &outlines}});
    Voronoi::WallDistancePointCloud cloud([field](const auto &points) { return field->sample_fields(points); });
    const BoundingBoxf3 bottom(Vec3d(-10., -10., 0.), Vec3d(10., 10., 5.));
    const BoundingBoxf3 center(Vec3d(-10., -10., 22.5), Vec3d(10., 10., 27.5));
    const BoundingBoxf3 top(Vec3d(-10., -10., 45.), Vec3d(10., 10., 50.));
    const auto floor_points = cloud.points_in(bottom, 1.);
    const auto center_points = cloud.points_in(center, 1.);
    const auto roof_points = cloud.points_in(top, 1.);
    REQUIRE(floor_points.size() > 3 * center_points.size());
    REQUIRE(roof_points.size() > 3 * center_points.size());
    REQUIRE(cloud.points_in(bottom, 1.) == floor_points);
}

TEST_CASE("Voronoi cloud setup is reused across parallel layers and resets on reslicing", "[FillVoronoi]")
{
    const ExPolygons outlines {rectangle(-20., -20., 20., 20.)};
    std::atomic<size_t> loads {0};
    const auto slices = [&] {
        ++loads;
        return std::vector<WallSlice>{{0., 20., &outlines}};
    };
    Voronoi::PointCloudSettings settings;
    settings.density.mode = AdaptiveInfillMode::SmoothBlend;
    Voronoi::WallDistanceCloudCache cache;
    std::vector<std::future<std::shared_ptr<const Voronoi::PointCloudProvider>>> workers;
    for (int i = 0; i < 8; ++i)
        workers.push_back(std::async(std::launch::async, [&] { return cache.get(slices, settings); }));
    const auto cloud = workers.front().get();
    for (size_t i = 1; i < workers.size(); ++i) REQUIRE(workers[i].get() == cloud);
    REQUIRE(cache.get(slices, settings) == cloud);
    REQUIRE(loads == 1);
    const BoundingBoxf3 region(Vec3d(-10., -10., 0.), Vec3d(10., 10., 3.));
    const auto points = cloud->points_in(region, 4.);
    REQUIRE_FALSE(points.empty());
    cache.clear();
    const auto regenerated = cache.get(slices, settings);
    REQUIRE(regenerated != cloud);
    REQUIRE(loads == 2);
    REQUIRE(regenerated->points_in(region, 4.) == points);
}

TEST_CASE("Voronoi Gaussian smoothing preserves flat-wall material density", "[FillVoronoi]")
{
    const ExPolygons outlines {rectangle(-30., -30., 30., 30.)};
    InfillDensityField field({{0., 60., &outlines}}, 2.);
    for (const auto &sample : field.sample_fields({Vec3d(0., 0., -2.), Vec3d(0., 0., 0.),
                Vec3d(0., 0., 2.), Vec3d(0., 0., 6.), Vec3d(28., 0., 30.)}))
        REQUIRE_THAT(sample.density, Catch::Matchers::WithinRel(
            InfillDensityField::density_ratio(sample.distance, sample.distance, 30., 0.2, AdaptiveInfillGradient::Linear), 0.04));
}

TEST_CASE("Voronoi Gaussian smoothing increases density at convex and concave corners", "[FillVoronoi]")
{
    const ExPolygons outlines {rectangle(-20., -20., 20., 20.)};
    InfillDensityField convex({{0., 40., &outlines}}, 2.);
    const auto outer = convex.sample_fields({Vec3d(14., 14., 20.), Vec3d(14., 0., 20.),
                                             Vec3d(14., 0., 34.), Vec3d(0., 14., 34.)});
    using Cloud = Voronoi::WallDistancePointCloud;
    REQUIRE(outer[0].density > outer[1].density);
    REQUIRE_THAT(outer[1].density, Catch::Matchers::WithinRel(InfillDensityField::density_ratio(outer[1].distance, outer[1].distance, 20., 0.2, AdaptiveInfillGradient::Linear), 0.04));
    REQUIRE_THAT(outer[2].density, Catch::Matchers::WithinRel(outer[0].density, 0.04));
    REQUIRE_THAT(outer[3].density, Catch::Matchers::WithinRel(outer[0].density, 0.04));
    ExPolygon hollow = rectangle(-30., -30., 30., 30.);
    hollow.holes.push_back(rectangle(-10., -10., 10., 10.).contour);
    hollow.holes.back().reverse();
    const ExPolygons hole {hollow};
    InfillDensityField concave({{0., 40., &hole}}, 2.);
    const auto inner = concave.sample_fields({Vec3d(14., 14., 20.)}).front();
    InfillDensityField unadjusted({{0., 40., &hole}}, 0.);
    REQUIRE(inner.density > unadjusted.sample_fields({Vec3d(14., 14., 20.)}).front().density);
}

TEST_CASE("Voronoi separable smoothing matches a three dimensional Gaussian", "[FillVoronoi]")
{
    const ExPolygons outlines {rectangle(-20., -20., 20., 20.)};
    InfillDensityField field({{0., 40., &outlines}}, 2.);
    const Vec3d center(14., 14., 34.);
    std::vector<Vec3d> points;
    std::vector<double> weights;
    for (int x = -6; x <= 6; ++x)
        for (int y = -6; y <= 6; ++y)
            for (int z = -6; z <= 6; ++z) {
                points.push_back(center + Vec3d(x, y, z));
                weights.push_back(std::exp(-(x * x + y * y + z * z) / 8.));
            }
    const auto distances = field.sample(points);
    double expected = 0., weight_sum = 0.;
    for (size_t i = 0; i < distances.size(); ++i) {
        expected += weights[i] * distances[i];
        weight_sum += weights[i];
    }
    REQUIRE_THAT(field.sample_fields({center}).front().density,
                 Catch::Matchers::WithinRel(InfillDensityField::density_ratio(field.sample({center}).front(), expected / weight_sum, 20., 0.2, AdaptiveInfillGradient::Linear), 0.04));
}

TEST_CASE("Voronoi Gaussian smoothing remains planar when its radius exceeds the halo", "[FillVoronoi]")
{
    const ExPolygons outlines {rectangle(-50., -50., 50., 50.)};
    InfillDensityField field({{0., 100., &outlines}}, 8.);
    const auto samples = field.sample_fields({Vec3d(0., 0., -2.), Vec3d(0., 0., 0.),
                                               Vec3d(0., 0., 4.), Vec3d(46., 0., 50.)});
    for (const auto &sample : samples)
        REQUIRE_THAT(sample.density, Catch::Matchers::WithinRel(InfillDensityField::density_ratio(sample.distance, sample.distance, 50., 0.2, AdaptiveInfillGradient::Linear), 0.06));
}

TEST_CASE("Voronoi Gaussian smoothing accepts zero width and rejects invalid widths", "[FillVoronoi]")
{
    const ExPolygons outlines {rectangle(-10., -10., 10., 10.)};
    InfillDensityField field({{0., 20., &outlines}}, 0.);
    const AdaptiveInfillField reference({{0., 20., &outlines}}, get_extents(outlines), AdaptiveInfillMode::SmoothBlend, [] {});
    for (const Vec3d &point : {Vec3d(0., 0., 10.), Vec3d(0., 0., -1.), Vec3d(9.5, 0., 10.)})
        REQUIRE_THAT(field.sample_fields({point}).front().density,
            Catch::Matchers::WithinRel(adaptive_density_scale(0.2, AdaptiveInfillGradient::Linear, reference.depth(point)), 0.025));
    REQUIRE_THROWS_AS(InfillDensityField({{0., 20., &outlines}}, -1.), std::invalid_argument);
    REQUIRE_THROWS_AS(InfillDensityField({{0., 20., &outlines}}, InfillDensityField::maximum_smoothing_sigma + 1.), std::invalid_argument);
    REQUIRE_THROWS_AS(InfillDensityField({{0., 20., &outlines}}, std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
}

TEST_CASE("Voronoi thinning samples prepared material density once per candidate", "[FillVoronoi]")
{
    size_t queried = 0;
    Voronoi::WallDistancePointCloud cloud([&](const auto &points) {
        queried += points.size();
        return std::vector<InfillDensitySample>(points.size(), {-6., InfillDensityField::density_ratio(-6., -5., 32.)});
    });
    const BoundingBoxf3 box(Vec3d(-10., -10., -10.), Vec3d(10., 10., 10.));
    const auto candidates = Voronoi::PoissonPointCloud().points_in(box, 3. / InfillDensityField::density_ratio(0., 2., 32.));
    cloud.points_in(box, 3.);
    REQUIRE(queried == candidates.size());
}

TEST_CASE("Voronoi gradients retain surface and interior targets", "[FillVoronoi]")
{
    using Cloud = Voronoi::WallDistancePointCloud;
    REQUIRE_THAT(InfillDensityField::density_ratio(2., 2., 32.), Catch::Matchers::WithinAbs(1., EPSILON));
    REQUIRE_THAT(InfillDensityField::density_ratio(-2., -2., 32.), Catch::Matchers::WithinAbs(std::pow(0.01, 2. / 32.), EPSILON));
    REQUIRE_THAT(InfillDensityField::density_ratio(-32., -32., 32.), Catch::Matchers::WithinAbs(0.01, EPSILON));
    REQUIRE(InfillDensityField::density_ratio(-8., -8., 32.) > InfillDensityField::density_ratio(-16., -16., 32.));
}

TEST_CASE("Voronoi corners boost density above the setting without boosting flat walls", "[FillVoronoi]")
{
    using Cloud = Voronoi::WallDistancePointCloud;
    for (double distance : {0., 1., 2., -1., -2.})
        REQUIRE_THAT(InfillDensityField::density_ratio(distance, distance, 32.), Catch::Matchers::WithinAbs(std::pow(0.01, std::max(0., -distance) / 32.), EPSILON));
    const double maximum = 0.01 + 0.99 * std::exp(2. / 6.);
    REQUIRE_THAT(InfillDensityField::density_ratio(-2., -1., 32.), Catch::Matchers::WithinAbs(std::pow(0.01, 1. / 32.) + 0.99 * std::expm1(1. / 6.), EPSILON));
    REQUIRE_THAT(InfillDensityField::density_ratio(-2., -3., 32.), Catch::Matchers::WithinAbs(InfillDensityField::density_ratio(-2., -1., 32.), EPSILON));
    REQUIRE_THAT(InfillDensityField::density_ratio(-1., -20., 32.), Catch::Matchers::WithinAbs(maximum, EPSILON));

    const BoundingBoxf3 box(Vec3d(-20., -20., -20.), Vec3d(20., 20., 20.));
    const auto uniform = Voronoi::PoissonPointCloud().points_in(box, 2.);
    auto count = [&](InfillDensitySample sample) {
        Cloud cloud([sample](const auto &points) { return std::vector<InfillDensitySample>(points.size(), sample); });
        return double(cloud.points_in(box, 2.).size()) / uniform.size();
    };
    REQUIRE_THAT(count({-1., 1.}), Catch::Matchers::WithinRel(1., 0.1));
    REQUIRE_THAT(count({-2., InfillDensityField::density_ratio(-2., -1., 32.)}), Catch::Matchers::WithinRel(std::pow(InfillDensityField::density_ratio(-2., -1., 32.), 3), 0.1));
    REQUIRE_THAT(count({-1., maximum}), Catch::Matchers::WithinRel(std::pow(maximum, 3), 0.1));
}

TEST_CASE("Uniform Voronoi density retains corner boosts", "[FillVoronoi]")
{
    using Cloud = Voronoi::WallDistancePointCloud;
    for (double depth : {0., 2., 12., 1000.}) {
        REQUIRE_THAT(InfillDensityField::density_ratio(-depth, -depth, 0.), Catch::Matchers::WithinAbs(1., EPSILON));
        REQUIRE(InfillDensityField::density_ratio(-depth, -depth + 1., 0.) > 1.);
    }
    REQUIRE(InfillDensityField::density_ratio(-12., -12., 64.) > InfillDensityField::density_ratio(-12., -12., 32.));
    const ExPolygons outlines {rectangle(-10., -10., 10., 10.)};
    Voronoi::WallDistanceCloudCache cache;
    size_t loads = 0;
    const auto slices = [&] { ++loads; return std::vector<WallSlice>{{0., 20., &outlines}}; };
    const auto plain = cache.get(slices, {{AdaptiveInfillMode::Disabled, 0., 0.25, 0.05, AdaptiveInfillGradient::Linear}, 0., VoronoiCloudMethod::Random, -1}, {}, {});
    const BoundingBoxf3 box(Vec3d(-5., -5., 0.), Vec3d(5., 5., 10.));
    REQUIRE(plain->points_in(box, 2.) == Voronoi::PoissonPointCloud().points_in(box, 2.));
    REQUIRE(cache.get(slices, {{AdaptiveInfillMode::Disabled, 0., 0.25, 0.05, AdaptiveInfillGradient::Linear}, 0., VoronoiCloudMethod::Random, -1}, {}, {}) == plain);
    REQUIRE(loads == 0);
    const auto graded = cache.get(slices, {{AdaptiveInfillMode::SmoothBlend, 0., 0.25, 0.05, AdaptiveInfillGradient::Linear}, 0., VoronoiCloudMethod::Random, -1}, {}, {});
    REQUIRE(graded != plain);
    REQUIRE(cache.get(slices, {{AdaptiveInfillMode::SmoothBlend, 0., 0.25, 0.05, AdaptiveInfillGradient::Linear}, 0., VoronoiCloudMethod::Random, -1}, {}, {}) == graded);
    REQUIRE(loads == 1);
    REQUIRE_THROWS_AS(Cloud({}, -1.), std::invalid_argument);
}

TEST_CASE("Voronoi distance construction reports progress through completion", "[FillVoronoi]")
{
    const ExPolygons outlines{rectangle(0., 0., 20., 20.)};
    const double sigma = GENERATE(0., 2.);
    auto previous_stage = WallDistanceStage::Distance;
    double previous_fraction = 0.;
    bool saw_smoothing = false;
    InfillDensityField field({{0., 20., &outlines}}, sigma,
        [&](WallDistanceStage stage, double fraction) {
            REQUIRE(int(stage) >= int(previous_stage));
            REQUIRE(fraction >= 0.);
            REQUIRE(fraction <= 1.);
            if (stage == previous_stage) REQUIRE(fraction >= previous_fraction);
            previous_stage = stage;
            previous_fraction = fraction;
            if (stage == WallDistanceStage::SmoothXY || stage == WallDistanceStage::SmoothZ)
                saw_smoothing = true;
        });
    REQUIRE(previous_stage == WallDistanceStage::Complete);
    REQUIRE_THAT(previous_fraction, Catch::Matchers::WithinAbs(1., 1e-12));
    REQUIRE(saw_smoothing == (sigma > 0.));
    REQUIRE(field.sample({Vec3d(10., 10., 10.)}).front() < 0.);
}

TEST_CASE("Canceled Voronoi smoothing leaves the cache usable", "[FillVoronoi]")
{
    const ExPolygons outlines{rectangle(0., 0., 20., 20.)};
    Voronoi::WallDistanceCloudCache cache;
    const auto slices = [&] { return std::vector<WallSlice>{{0., 20., &outlines}}; };
    REQUIRE_THROWS_AS(cache.get(slices, {{AdaptiveInfillMode::SmoothBlend, 2., 0.25, 0.05, AdaptiveInfillGradient::Linear}, 0., VoronoiCloudMethod::Random, -1}, [](WallDistanceStage stage, double) {
            if (stage == WallDistanceStage::SmoothXY) throw std::runtime_error("cancel");
        }, {}), std::runtime_error);
    bool completed = false;
    const auto provider = cache.get(slices, {{AdaptiveInfillMode::SmoothBlend, 2., 0.25, 0.05, AdaptiveInfillGradient::Linear}, 0., VoronoiCloudMethod::Random, -1}, [&](WallDistanceStage stage, double fraction) {
            if (stage == WallDistanceStage::Complete && fraction == 1.) completed = true;
        }, {});
    REQUIRE(completed);
    REQUIRE(provider != nullptr);
    REQUIRE(cache.get(slices, {{AdaptiveInfillMode::SmoothBlend, 2., 0.25, 0.05, AdaptiveInfillGradient::Linear}, 0., VoronoiCloudMethod::Random, -1}, [](auto, double) { throw std::runtime_error("cached fields must not rebuild"); }, {}) == provider);
}






TEST_CASE("FCC has twelve nearest neighbours and a vertical body diagonal", "[FillVoronoi][Fcc]")
{
    const Voronoi::LatticePointCloud cloud(VoronoiCloudMethod::Fcc, {}, 1.);
    const auto sites = cloud.points_in(BoundingBoxf3(Vec3d::Constant(-8.), Vec3d::Constant(8.)), 3.);
    std::vector<Vec3d> neighbours;
    for (const auto &p : sites)
        if (p.norm() > 1e-9) neighbours.push_back(p);
    std::sort(neighbours.begin(), neighbours.end(), [](const auto &a, const auto &b) { return a.norm() < b.norm(); });
    REQUIRE(neighbours.size() > 12);
    const double edge = std::pow(2., 1. / 6.) * 3.;
    int vertical_faces = 0;
    for (size_t i = 0; i < 12; ++i) {
        REQUIRE_THAT(neighbours[i].norm(), Catch::Matchers::WithinAbs(edge, 1e-9));
        const double cosine = std::abs(neighbours[i].z()) / edge;
        if (cosine < 1e-9) ++vertical_faces;
        else REQUIRE_THAT(cosine, Catch::Matchers::WithinAbs(std::sqrt(2. / 3.), 1e-9));
    }
    REQUIRE(vertical_faces == 6);
    REQUIRE(neighbours[12].norm() > edge * 1.1);
}

TEST_CASE("Cubic sites have six orthogonal neighbours and a vertical body diagonal", "[FillVoronoi][Lattice]")
{
    const Voronoi::LatticePointCloud cloud(VoronoiCloudMethod::Cubic, {}, 1.);
    const auto sites = cloud.points_in(BoundingBoxf3(Vec3d::Constant(-8.), Vec3d::Constant(8.)), 3.);
    std::vector<Vec3d> neighbours;
    for (const auto &p : sites)
        if (p.norm() > 1e-9) neighbours.push_back(p);
    std::sort(neighbours.begin(), neighbours.end(), [](const auto &a, const auto &b) { return a.norm() < b.norm(); });
    REQUIRE(neighbours.size() > 6);
    for (size_t i = 0; i < 6; ++i) {
        REQUIRE_THAT(neighbours[i].norm(), Catch::Matchers::WithinAbs(3., 1e-9));
        REQUIRE_THAT(std::abs(neighbours[i].z()) / 3., Catch::Matchers::WithinAbs(1. / std::sqrt(3.), 1e-9));
        for (size_t j = i + 1; j < 6; ++j) {
            const double cosine = neighbours[i].dot(neighbours[j]) / 9.;
            REQUIRE_THAT(cosine, Catch::Matchers::WithinAbs(cosine < -0.5 ? -1. : 0., 1e-9));
        }
    }
    REQUIRE_THAT(neighbours[6].norm(), Catch::Matchers::WithinAbs(3. * std::sqrt(2.), 1e-9));
}

TEST_CASE("FCC thinning follows cubic tetrahedra and coarser checkerboards", "[FillVoronoi][Fcc]")
{
    const double retention = GENERATE(0.75, 0.5, 0.25, 0.125, 0.015625);
    const double maximum_depth = 6., spacing = 3.;
    const double near = InfillDensityField::near_wall_distance;
    const double maximum = InfillDensityField::density_ratio(0., near, 0.);
    const double target = maximum * std::cbrt(retention);
    // Find a constant distance sample giving this retention through the public
    // density function, independently of its particular maximum_depth curve.
    double low = 0., high = 64.;
    for (int step = 0; step < 48; ++step) {
        const double depth = 0.5 * (low + high);
        if (InfillDensityField::density_ratio(-std::max(near, depth), -depth, maximum_depth) > target)
            low = depth;
        else high = depth;
    }
    const double depth = 0.5 * (low + high);
    const auto samples = [=](const auto &points) {
        return std::vector<InfillDensitySample>(points.size(), {-std::max(near, depth), InfillDensityField::density_ratio(-std::max(near, depth), -depth, maximum_depth)});
    };
    const auto saturated = [=](const auto &points) {
        return std::vector<InfillDensitySample>(points.size(), {0., InfillDensityField::density_ratio(0., near, maximum_depth)});
    };
    const BoundingBoxf3 box(Vec3d::Constant(-12.), Vec3d::Constant(12.));
    const auto selected = Voronoi::LatticePointCloud(VoronoiCloudMethod::Fcc, samples, maximum).points_in(box, spacing);
    const auto candidates = Voronoi::LatticePointCloud(VoronoiCloudMethod::Fcc, saturated, maximum).points_in(box, spacing);
    Eigen::Matrix3d orientation;
    orientation.row(0) = Vec3d(1., -1., 0.).normalized();
    orientation.row(1) = Vec3d(1., 1., -2.).normalized();
    orientation.row(2) = Vec3d(1., 1., 1.).normalized();
    const double half_cube = std::pow(2., 1. / 6.) * spacing / maximum / std::sqrt(2.);
    const auto coordinates = [&](const Vec3d &p) {
        const Vec3d untilted = orientation.transpose() * p;
        return std::array<int64_t, 3>{std::llround(untilted.x() / half_cube), std::llround(untilted.y() / half_cube),
                                      std::llround(untilted.z() / half_cube)};
    };
    std::set<std::array<int64_t, 3>> actual, expected;
    for (const auto &p : selected) REQUIRE(actual.insert(coordinates(p)).second);
    for (const auto &p : candidates) {
        const auto q = coordinates(p);
        const bool u = q[0] % 2 != 0, v = q[1] % 2 != 0, w = q[2] % 2 != 0;
        const int64_t scale = retention == 0.015625 ? 4 : 2;
        const bool keep = retention == 0.75 ? !(u && v) : retention == 0.5 ? !u :
                          retention == 0.25 ? !u && !v && !w :
                          q[0] % scale == 0 && q[1] % scale == 0 && q[2] % scale == 0 &&
                          (q[0] + q[1] + q[2]) / scale % 2 == 0;
        if (keep) expected.insert(q);
    }
    REQUIRE_FALSE(expected.empty());
    REQUIRE(actual == expected);
}

TEST_CASE("Lattice thinning follows density and preserves partitioned queries", "[FillVoronoi][Lattice]")
{
    const auto method = GENERATE(VoronoiCloudMethod::Fcc, VoronoiCloudMethod::Cubic);
    const bool adaptive = GENERATE(false, true);
    const double maximum = InfillDensityField::density_ratio(0., 2., 0.);
    const auto sample = [adaptive](const std::vector<Vec3d> &points) {
        return std::vector<InfillDensitySample>(points.size(), {-6., InfillDensityField::density_ratio(-6., -6., adaptive ? 32. : 0.)});
    };
    const Voronoi::LatticePointCloud cloud(method, sample, maximum);
    const BoundingBoxf3 box(Vec3d(-50.3, -49.7, -48.8), Vec3d(51.2, 50.4, 49.9));
    const auto whole = cloud.points_in(box, 3.);
    REQUIRE_FALSE(whole.empty());
    auto keys = [](const auto &points) {
        std::set<std::array<double, 3>> result;
        for (const auto &p : points) REQUIRE(result.insert({p.x(), p.y(), p.z()}).second);
        return result;
    };
    const int axis = GENERATE(0, 1, 2);
    BoundingBoxf3 first = box, second = box;
    first.max[axis] = second.min[axis] = 0.37;
    auto parts = cloud.points_in(first, 3.);
    append(parts, cloud.points_in(second, 3.));
    REQUIRE(keys(parts) == keys(whole));
    const auto candidates = Voronoi::LatticePointCloud(method, {}, 1.).points_in(box, 3. / maximum);
    const double probability = std::pow(InfillDensityField::density_ratio(-6., -6., adaptive ? 32. : 0.) / maximum, 3.);
    REQUIRE_THAT(double(whole.size()) / candidates.size(), Catch::Matchers::WithinRel(probability, 0.08));
    REQUIRE(cloud.points_in(box, 3.) == whole);
}

TEST_CASE("The candidate density bound covers both sparse interiors and dense walls", "[FillVoronoi][MaterialDensity]")
{
    const ExPolygons model {rectangle(-40., -40., 40., 40.)};
    const bool adaptive = true;
    const double sigma = GENERATE(0., 8.);
    const double density = 0.3;
    InfillDensityField field({{-40., 40., &model}}, sigma, {}, adaptive ? AdaptiveInfillMode::SmoothBlend : AdaptiveInfillMode::Disabled, density);
    std::vector<Vec3d> points;
    for (int x = 0; x <= 40; ++x) points.emplace_back(x, 0., 0.);
    const auto samples = field.sample_fields(points);
    const double boost = InfillDensityField::density_ratio(0., sigma == 0. ? 0. : 2., 32.);
    const double expected = boost;
    REQUIRE_THAT(field.maximum_density_ratio(), Catch::Matchers::WithinRel(expected, 1e-12));
    REQUIRE(samples.front().density < samples.back().density);
    for (const auto &sample : samples)
        REQUIRE(field.maximum_density_ratio() >= sample.density - 1e-12);
}

TEST_CASE("Fused material density uses an inside-only smoothstep including modifier ceilings", "[FillVoronoi][MaterialDensity]")
{
    const ExPolygons model {rectangle(-20., -20., 20., 20.)};
    DensityModifier modifier {3, 0.5, 4., {{8., 32., {rectangle(0., -20., 20., 20.)}}}};
    modifier.interior_density = modifier.material_density;
    InfillDensityField field({{0., 40., &model}}, 0., {}, AdaptiveInfillMode::Disabled, 0.1, {modifier}, 3);
    for (double depth : {-2., 0., 1., 2., 3., 4., 8.}) {
        const double t = std::clamp(depth / 4., 0., 1.);
        const double expected = 0.1 + 0.4 * t * t * (3. - 2. * t);
        const auto samples = field.sample_fields({Vec3d(depth, 0., 20.), Vec3d(10., 0., 8. + depth)});
        // Allow for the planar boundary and interpolation error in the blend mask.
        const double tolerance = 1.5 * 0.4 * 0.5 * field.grid_spacing() / 4. + 0.025 * 0.5;
        for (const auto &sample : samples)
            REQUIRE_THAT(0.1 * sample.density, Catch::Matchers::WithinAbs(expected, tolerance));
    }
    REQUIRE(field.maximum_density_ratio() >= 5.);
    InfillDensityField other_region({{0., 40., &model}}, 0., {}, AdaptiveInfillMode::Disabled, 0.1, {modifier}, 7);
    REQUIRE_THAT(other_region.sample_fields({Vec3d(10., 0., 20.)}).front().density,
                 Catch::Matchers::WithinRel(1., 0.022));
}

TEST_CASE("Fused modifiers compose material density in order", "[FillVoronoi][MaterialDensity]")
{
    const ExPolygons model {rectangle(-20., -20., 20., 20.)};
    DensityModifier first {3, 0.5, 0., {{0., 40., {rectangle(0., -20., 20., 20.)}}}};
    first.interior_density = first.material_density;
    auto second = first;
    second.material_density = 0.2;
    second.interior_density = second.material_density;
    InfillDensityField ordered({{0., 40., &model}}, 0., {}, AdaptiveInfillMode::Disabled, 0.1, {first, second}, 3);
    const auto sample = ordered.sample_fields({Vec3d(10., 0., 20.)}).front();
    REQUIRE_THAT(sample.density * 0.1, Catch::Matchers::WithinRel(0.2, 0.022));
    first.blend_width = second.blend_width = 4.;
    InfillDensityField blended({{0., 40., &model}}, 0., {}, AdaptiveInfillMode::Disabled, 0.1, {first, second}, 3);
    InfillDensityField first_only({{0., 40., &model}}, 0., {}, AdaptiveInfillMode::Disabled, 0.1, {first}, 3);
    const Vec3d point(2., 0., 20.);
    const double first_value = first_only.density(point);
    const double blend = (first_value - 0.1) / (first.material_density - 0.1);
    const double expected = first_value * (1. - blend) + second.material_density * blend;
    REQUIRE_THAT(blended.density(point), Catch::Matchers::WithinRel(expected, 0.05));
}

TEST_CASE("Fused point clouds accommodate higher modifier density and empty interiors", "[FillVoronoi][MaterialDensity]")
{
    const auto method = GENERATE(VoronoiCloudMethod::Random, VoronoiCloudMethod::Fcc, VoronoiCloudMethod::Cubic);
    const double target = GENERATE(0., 0.5);
    const ExPolygons model {rectangle(-20., -20., 20., 20.)};
    DensityModifier modifier {3, target, 0., {{0., 40., {rectangle(0., -20., 20., 20.)}}}};
    if (target > 0.) modifier.interior_density = target;
    Voronoi::WallDistanceCloudCache cache;
    const auto slices = [&] { return std::vector<WallSlice>{{0., 40., &model}}; };
    const auto field = cache.get(slices, {{AdaptiveInfillMode::Disabled, 0., 0.1, 0.05, AdaptiveInfillGradient::Linear}, 0., method, 3}, {}, {modifier});
    const BoundingBoxf3 box(Vec3d(4., -5., 15.), Vec3d(14., 5., 25.));
    const auto points = field->points_in(box, 3.);
    if (target == 0.) REQUIRE(points.empty());
    else {
        const auto uniform = cache.get(slices, {{AdaptiveInfillMode::Disabled, 0., 0.25, 0.05, AdaptiveInfillGradient::Linear}, 0., method, -1}, {}, {})->points_in(box, 3.);
        REQUIRE(points.size() > 50 * uniform.size());
    }
    const auto left = field->points_in(BoundingBoxf3(box.min, Vec3d(9., 5., 25.)), 3.);
    const auto right = field->points_in(BoundingBoxf3(Vec3d(9., -5., 15.), box.max), 3.);
    REQUIRE(left.size() + right.size() == points.size());
}

TEST_CASE("Modifier blend weights use smoothstep entirely inside the boundary", "[FillVoronoi][MaterialDensity]")
{
    REQUIRE_THAT(Slic3r::modifier_blend(-1., 4.), Catch::Matchers::WithinAbs(0., EPSILON));
    REQUIRE_THAT(Slic3r::modifier_blend(0., 4.), Catch::Matchers::WithinAbs(0., EPSILON));
    REQUIRE_THAT(Slic3r::modifier_blend(1., 4.), Catch::Matchers::WithinAbs(0.15625, EPSILON));
    REQUIRE_THAT(Slic3r::modifier_blend(2., 4.), Catch::Matchers::WithinAbs(0.5, EPSILON));
    REQUIRE_THAT(Slic3r::modifier_blend(3., 4.), Catch::Matchers::WithinAbs(0.84375, EPSILON));
    REQUIRE_THAT(Slic3r::modifier_blend(4., 4.), Catch::Matchers::WithinAbs(1., EPSILON));
    REQUIRE_THAT(Slic3r::modifier_blend(8., 4.), Catch::Matchers::WithinAbs(1., EPSILON));
    REQUIRE_THAT(Slic3r::modifier_blend(0., 0.), Catch::Matchers::WithinAbs(1., EPSILON));
}

TEST_CASE("Fused modifiers control interior density and corner smoothing", "[FillVoronoi][MaterialDensity]")
{
    const double parent_sigma = GENERATE(0., 8.);
    const ExPolygons model {rectangle(-20., -20., 20., 20.)};
    DensityModifier modifier {3, 0.1, 0., {{0., 40., model}}};
    modifier.interior_density = modifier.material_density;
    modifier.smoothing_sigma = 0.;
    InfillDensityField uniform({{0., 40., &model}}, parent_sigma, {}, AdaptiveInfillMode::SmoothBlend, 0.1, {modifier}, 3);
    const auto center = uniform.sample_fields({Vec3d(0., 0., 20.)}).front();
    REQUIRE_THAT(center.density, Catch::Matchers::WithinRel(1., 0.022));
    const Vec3d corner(14., 14., 20.);
    const auto unadjusted = uniform.sample_fields({corner}).front();
    modifier.smoothing_sigma = 8.;
    InfillDensityField smoothed({{0., 40., &model}}, parent_sigma, {}, AdaptiveInfillMode::SmoothBlend, 0.1, {modifier}, 3);
    const auto adjusted = smoothed.sample_fields({corner}).front();
    REQUIRE(adjusted.density > unadjusted.density + 0.1);
    REQUIRE_THAT(adjusted.distance, Catch::Matchers::WithinAbs(unadjusted.distance, EPSILON));
}

TEST_CASE("Incremental Voronoi geometry agrees with fresh builds while the slab moves", "[FillVoronoi][Incremental]")
{
    const auto provider = std::make_shared<Voronoi::PoissonPointCloud>();
    const auto bounds = rectangle(-6., -6., 6., 6.).contour.bounding_box();
    Voronoi::CellGeometry mesh;
    size_t removed = 0;
    for (double z : {0.1, 1.1, 2.1, 3.1, 5.1, 8.1, 12.1, 0.1}) {
        const auto points = provider->points_in(BoundingBoxf3(Vec3d(-14., -14., z - 7.), Vec3d(14., 14., z + 7.)), 4.);
        const auto changed = mesh.update(points);
        removed += changed.second;
        const auto unchanged = mesh.update(points);
        REQUIRE(unchanged.first == 0);
        REQUIRE(unchanged.second == 0);
        const auto bands = mesh.shallow_bands(bounds, z, 0.45, 0.2, 0.4);
        Voronoi::CellGeometry reference_mesh;
        reference_mesh.update(points);
        const auto reference = reference_mesh.shallow_bands(bounds, z, 0.45, 0.2, 0.4);
        ExPolygons got, want;
        for (const auto &band : bands) got.push_back(band.polygon);
        for (const auto &band : reference) want.push_back(band.polygon);
        ExPolygons difference = diff_ex(got, want);
        append(difference, diff_ex(want, got));
        double area = 0.;
        for (const auto &polygon : difference) area += polygon.area() * SCALING_FACTOR * SCALING_FACTOR;
        REQUIRE(area < 0.0001);
    }
    REQUIRE(removed > 0);
}

TEST_CASE("Copied Voronoi networks own their storage through dimension changes", "[FillVoronoi][Incremental]")
{
    const auto bounds = rectangle(-10., -10., 10., 10.).contour.bounding_box();
    Voronoi::CellGeometry mesh;
    std::vector<Vec3d> points {Vec3d(-5., 0., 0.), Vec3d(5., 0., 0.), Vec3d(0., 5., 0.), Vec3d(0., 0., 10.)};
    const auto contours = [&](const Voronoi::CellGeometry &geometry) {
        std::vector<Points> result;
        for (const auto &band : geometry.shallow_bands(bounds, 4.9, 0.45, 0.2, 0.4))
            result.push_back(band.polygon.contour.points);
        return result;
    };
    mesh.update(points);
    Voronoi::CellGeometry copy = mesh, assigned;
    assigned = mesh;
    const auto expected = contours(copy);
    REQUIRE_FALSE(expected.empty());
    while (!points.empty()) {
        points.pop_back();
        mesh.update(points);
        Voronoi::CellGeometry fresh;
        fresh.update(points);
        REQUIRE(contours(mesh) == contours(fresh));
        REQUIRE(contours(copy) == expected);
        REQUIRE(contours(assigned) == expected);
    }
    REQUIRE(contours(mesh).empty());
}

TEST_CASE("Candidate density bounds include modifier settings and ignore other regions", "[FillVoronoi][MaterialDensity]")
{
    const ExPolygons model {rectangle(-10., -10., 10., 10.)};
    const double density = GENERATE(0., 0.7, 0.9);
    const double sigma = GENERATE(0., 8.);
    DensityModifier modifier {3, density, 0., {{0., 20., model}}};
    modifier.smoothing_sigma = sigma;
    const DensityModifier unrelated {7, 1., 0., {{0., 20., model}}};
    InfillDensityField field({{0., 20., &model}}, 0., {}, AdaptiveInfillMode::SmoothBlend, 0.1, {modifier, unrelated}, 3);
    const double boost = InfillDensityField::density_ratio(0., sigma == 0. ? 0. : 2., 32.);
    const double maximum = std::max(0.1, std::min(1., density * boost));
    const double expected = std::max(maximum, double(float(maximum))) / 0.1;
    REQUIRE_THAT(field.maximum_density_ratio(), Catch::Matchers::WithinRel(expected, 1e-12));
}

TEST_CASE("Voronoi uses the shared adaptive density gradients over body-relative depth", "[FillVoronoi]")
{
    using Cloud = Voronoi::WallDistancePointCloud;
    const auto gradient = GENERATE(AdaptiveInfillGradient::Linear, AdaptiveInfillGradient::Quadratic, AdaptiveInfillGradient::Exponential);
    const double middle = gradient == AdaptiveInfillGradient::Linear ? 0.6 :
        gradient == AdaptiveInfillGradient::Quadratic ? 0.8 : std::sqrt(0.2);
    REQUIRE_THAT(InfillDensityField::density_ratio(0., 0., 40., 0.2, gradient), Catch::Matchers::WithinAbs(1., EPSILON));
    REQUIRE_THAT(InfillDensityField::density_ratio(-20., -20., 40., 0.2, gradient), Catch::Matchers::WithinAbs(middle, EPSILON));
    REQUIRE_THAT(InfillDensityField::density_ratio(-40., -40., 40., 0.2, gradient), Catch::Matchers::WithinAbs(0.2, EPSILON));
    REQUIRE_THAT(InfillDensityField::density_ratio(-1000., -1000., 40., 0.2, gradient), Catch::Matchers::WithinAbs(0.2, EPSILON));
    REQUIRE_THAT(InfillDensityField::density_ratio(-20., -20., 0., 0.2, gradient), Catch::Matchers::WithinAbs(1., EPSILON));
    REQUIRE(InfillDensityField::density_ratio(-20., -18., 40., 2., gradient) >= InfillDensityField::density_ratio(-20., -20., 40., 2., gradient));
}

TEST_CASE("Voronoi prepared interior density is absolute and each connected body reaches it", "[FillVoronoi]")
{
    const ExPolygons model {rectangle(-50., -10., -30., 10.), rectangle(-10., -30., 50., 30.)};
    const std::vector<WallSlice> slices {{0., 80., &model}};
    const AdaptiveInfillField reference(slices, get_extents(model), AdaptiveInfillMode::SmoothBlend, [] {});
    const std::vector<Vec3d> points {Vec3d(-40., 0., 40.), Vec3d(20., 0., 40.)};
    for (double surface : {0.1, 0.25}) {
        InfillDensityField field({{0., 80., &model}}, 0., {}, AdaptiveInfillMode::SmoothBlend, surface);
        const auto values = field.sample_fields(points);
        for (size_t i = 0; i < points.size(); ++i) {
            const double depth = reference.depth(points[i]);
            REQUIRE(depth > 0.95);
            const double expected = surface + (0.05 - surface) * depth;
            REQUIRE_THAT(values[i].density * surface, Catch::Matchers::WithinAbs(expected, 0.002));
        }
    }
}

TEST_CASE("Voronoi unsmoothed material density follows TPMS depth when the model is scaled", "[FillVoronoi]")
{
    for (double half : {10., 30.}) {
        const ExPolygons model {rectangle(-half, -half, half, half)};
        InfillDensityField field({{0., 2. * half, &model}}, 0.);
        const auto values = field.sample_fields({Vec3d(0., 0., 0.), Vec3d(0., 0., half / 2.), Vec3d(0., 0., half)});
        REQUIRE_THAT(values[0].density, Catch::Matchers::WithinRel(1., 0.04));
        REQUIRE_THAT(values[1].density, Catch::Matchers::WithinRel(0.6, 0.04));
        REQUIRE_THAT(values[2].density, Catch::Matchers::WithinRel(0.2, 0.04));
    }
}

TEST_CASE("Fused Voronoi interior targets and gradients override the parent without exceeding the candidate bound", "[FillVoronoi][Modifiers]")
{
    const ExPolygons model {rectangle(-50., -50., 50., 50.)};
    DensityModifier modifier {3, 0.1, 0., {{0., 100., {rectangle(-5., -5., 5., 5.)}}}};
    modifier.interior_density = 0.9;
    modifier.gradient = AdaptiveInfillGradient::Quadratic;
    InfillDensityField field({{0., 100., &model}}, 0., {}, AdaptiveInfillMode::SmoothBlend, 0.1, {modifier}, 3);
    const auto values = field.sample_fields({Vec3d(0., 0., 25.), Vec3d(0., 0., 50.)});
    REQUIRE_THAT(values[0].density * 0.1, Catch::Matchers::WithinRel(0.3, 0.04));
    REQUIRE_THAT(values[1].density * 0.1, Catch::Matchers::WithinRel(0.9, 0.04));
    REQUIRE(values[1].density <= field.maximum_density_ratio());
}

TEST_CASE("Unbounded empty Voronoi providers stop expanding and report worker errors", "[FillVoronoi]")
{
    class EmptyProvider : public Voronoi::PointCloudProvider {
    public:
    double spacing(double s, double h, double d) const override
        { return Voronoi::PoissonPointCloud::nominal_spacing(s, h, d); }
        mutable unsigned queries = 0;
        std::vector<Vec3d> points_in(const BoundingBoxf3 &, double) const override
            { ++queries; return {}; }
    };
    auto provider = std::make_shared<EmptyProvider>();
    Voronoi::ActivePointCloud cloud(provider);
    const BoundingBox bounds = rectangle(-1., -1., 1., 1.).contour.bounding_box();
    // Catch outside TBB: slicing must receive the original, readable exception.
    REQUIRE_THROWS_WITH(tbb::parallel_for(size_t(0), size_t(1), [&](size_t) {
        cloud.update(bounds, 0., 0., 1., [](const auto &) {
            return std::numeric_limits<double>::infinity();
        });
    }), "Voronoi point sampling found no control points; check infill density and modifiers");
    REQUIRE(provider->queries <= 19);

    Voronoi::ActivePointCloud finite(std::make_shared<FixedProvider>());
    REQUIRE_FALSE(finite.update(bounds, 0., 0., 1., [](const auto &) {
        return std::numeric_limits<double>::infinity();
    }));
}

TEST_CASE("An all-zero fused modifier cannot expand the active cloud indefinitely", "[FillVoronoi]")
{
    const ExPolygons model {rectangle(-1., -1., 1., 1.)};
    DensityModifier modifier {3, 0., 0., {{-20., 20., {rectangle(-20., -20., 20., 20.)}}}};
    auto field = std::make_shared<InfillDensityField>(std::vector<WallSlice>{{-1., 1., &model}},
        0., WallDistanceProgress{}, AdaptiveInfillMode::Disabled, 0.15, std::vector<DensityModifier>{modifier}, 3);
    auto provider = std::make_shared<Voronoi::WallDistancePointCloud>([field](const auto &points) {
        return field->sample_fields(points);
    }, field->maximum_density() / 0.15);
    Voronoi::ActivePointCloud cloud(provider);
    REQUIRE_THROWS_WITH(cloud.update(get_extents(model), 0., 0., 0.5, [](const auto &) {
        return std::numeric_limits<double>::infinity();
    }), "Voronoi point sampling found no control points; check infill density and modifiers");
    REQUIRE(cloud.sites().empty());
}

TEST_CASE("Active cloud retains sites ingested above its current slab", "[FillVoronoi]")
{
    class RecordingCloud : public Voronoi::PoissonPointCloud {
    public:
        mutable std::vector<BoundingBoxf3> queries;
        std::vector<Vec3d> points_in_cells(const Vec3i64 &first, const Vec3i64 &last, double spacing) const override {
            queries.emplace_back(first.cast<double>() * cell_size(spacing), last.cast<double>() * cell_size(spacing));
            return PoissonPointCloud::points_in_cells(first, last, spacing);
        }
    };
    auto provider = std::make_shared<RecordingCloud>();
    Voronoi::ActivePointCloud cloud(provider);
    const BoundingBox bounds(Point::new_scale(-2., -2.), Point::new_scale(2., 2.));
    REQUIRE(cloud.update(bounds, 0., 0., 2., [](const auto &) { return 0.1; }));
    const double top = cloud.known_bounds().max.z();
    provider->queries.clear();
    REQUIRE(cloud.update(bounds, 0.05, 0., 2., [](const auto &) { return 0.1; }));
    REQUIRE(cloud.known_bounds().max.z() >= top);
    for (const auto &query : provider->queries)
        REQUIRE(query.min.z() >= top);
}

TEST_CASE("Closed Voronoi trails scatter seams rather than advancing around the loop", "[FillVoronoi]")
{
    Points ring;
    for (int i = 0; i < 64; ++i) {
        const double angle = 2. * M_PI * i / 64.;
        ring.emplace_back(coord_t(10000. * std::cos(angle)), coord_t(10000. * std::sin(angle)));
    }
    Polylines walls;
    for (size_t i = 0; i < ring.size(); ++i) walls.emplace_back(ring[i], ring[(i + 1) % ring.size()]);
    const auto graph = Voronoi::build_wall_graph(walls);
    const auto connections = Voronoi::pair_junctions(graph, {}, 0);
    std::set<Point> starts;
    double largest_jump = 0.;
    Point previous;
    for (size_t layer = 0; layer < 32; ++layer) {
        const auto trails = Voronoi::trace_trails(graph, connections, layer);
        REQUIRE(trails.size() == 1);
        const Point start = trails.front().first_point();
        starts.insert(start);
        if (layer) largest_jump = std::max(largest_jump, (start - previous).cast<double>().norm());
        previous = start;
    }
    REQUIRE(starts.size() > 8);
    REQUIRE(largest_jump > 10000.);
}

TEST_CASE("Default Voronoi density keeps a uniform cloud without loading geometry", "[FillVoronoi][Regression]")
{
    Voronoi::WallDistanceCloudCache cache;
    size_t loads = 0;
    const auto provider = cache.get([&] { ++loads; return std::vector<WallSlice>{}; });
    REQUIRE(loads == 0);
    const BoundingBoxf3 box(Vec3d(-5., -5., 0.), Vec3d(5., 5., 10.));
    REQUIRE(provider->points_in(box, 2.) == Voronoi::PoissonPointCloud().points_in(box, 2.));
}

TEST_CASE("Voronoi point providers cache each adaptive spatial mode independently", "[FillVoronoi][Regression]")
{
    const ExPolygons model {rectangle(-12., -8., 12., 8.)};
    const auto slices = [&] { return std::vector<WallSlice>{{0., 20., &model}}; };
    Voronoi::WallDistanceCloudCache cache;
    Voronoi::PointCloudSettings settings;
    settings.density.mode = AdaptiveInfillMode::NormalX;
    const auto first = cache.get(slices, settings);
    settings.density.mode = AdaptiveInfillMode::NormalY;
    const auto second = cache.get(slices, settings);
    REQUIRE(second != first);
    settings.density.mode = AdaptiveInfillMode::NormalX;
    REQUIRE(cache.get(slices, settings) == first);
}

TEST_CASE("Whole generation cubes are visited once and retired before the next batch", "[FillVoronoi][Incremental]")
{
    class RecordedCloud : public Voronoi::PoissonPointCloud {
    public:
        mutable std::set<std::array<int64_t, 3>> visited;
        mutable size_t calls = 0;
        mutable bool retired = false;
        mutable bool saw_retirement = false;
        std::vector<Vec3d> points_in_cells(const Vec3i64 &first, const Vec3i64 &last, double spacing) const override
        {
            ++calls;
            if (retired) { saw_retirement = true; retired = false; }
            for (int64_t x = first.x(); x < last.x(); ++x)
                for (int64_t y = first.y(); y < last.y(); ++y)
                    for (int64_t z = first.z(); z < last.z(); ++z)
                        REQUIRE(visited.insert({x, y, z}).second);
            return PoissonPointCloud::points_in_cells(first, last, spacing);
        }
    };
    auto provider = std::make_shared<RecordedCloud>();
    Voronoi::ActivePointCloud cloud(provider);
    const auto bounds = rectangle(-2., -2., 2., 2.).contour.bounding_box();
    Voronoi::CellGeometry mesh;
    const auto radius = [&](const auto &sites) {
        mesh.update(sites);
        return 1.;
    };
    const auto retire = [&](double cutoff) {
        mesh.remove_below(cutoff);
        if (std::isfinite(cutoff)) provider->retired = true;
    };
    REQUIRE(cloud.update(bounds, 0.1, 0., 2., radius, retire));
    const size_t initial_calls = provider->calls, revision = cloud.revision();
    REQUIRE(cloud.update(bounds, 0.2, 0., 2., radius, retire));
    REQUIRE(provider->calls == initial_calls);
    REQUIRE(cloud.revision() == revision);
    for (double z = 0.4; z < 30.; z += 0.2)
        REQUIRE(cloud.update(bounds, z, 0., 2., radius, retire));
    REQUIRE(provider->saw_retirement);
    REQUIRE(cloud.retired_below() > 0.);
}

TEST_CASE("Batch Delaunay retirement agrees with a fresh retained cloud", "[FillVoronoi][Incremental]")
{
    Voronoi::PoissonPointCloud provider;
    auto sites = provider.points_in(BoundingBoxf3(Vec3d(-8., -8., -8.), Vec3d(8., 8., 12.)), 3.);
    Voronoi::CellGeometry mesh;
    mesh.update(sites);
    Voronoi::CellGeometry copy = mesh;
    const size_t before = sites.size();
    sites.erase(std::remove_if(sites.begin(), sites.end(), [](const auto &p) { return p.z() < 0.; }), sites.end());
    REQUIRE(mesh.remove_below(0.) == before - sites.size());
    REQUIRE(mesh.update(sites) == std::make_pair(size_t(0), size_t(0)));
    Voronoi::CellGeometry fresh;
    fresh.update(sites);
    const auto bounds = rectangle(-4., -4., 4., 4.).contour.bounding_box();
    ExPolygons got, expected;
    for (const auto &band : mesh.shallow_bands(bounds, 5., 0.45, 0.2, 0.4)) got.push_back(band.polygon);
    for (const auto &band : fresh.shallow_bands(bounds, 5., 0.45, 0.2, 0.4)) expected.push_back(band.polygon);
    REQUIRE(diff_ex(got, expected).empty());
    REQUIRE(diff_ex(expected, got).empty());
    REQUIRE(copy.remove_below(std::numeric_limits<double>::infinity()) == before);
    REQUIRE(copy.remove_below(0.) == 0);
}

TEST_CASE("A wider slice slab does not reopen safely retired generation cells", "[FillVoronoi][Incremental]")
{
    Voronoi::ActivePointCloud cloud;
    const auto bounds = rectangle(-2., -2., 2., 2.).contour.bounding_box();
    const auto radius = [](const auto &) { return 1.; };
    for (double z = 0.; z < 8.4; z += 0.2)
        REQUIRE(cloud.update(bounds, z, 0., 2., radius));
    const double cutoff = cloud.retired_below();
    REQUIRE(cutoff > 0.);
    size_t probes = 0;
    REQUIRE(cloud.update(bounds, 14., 5., 2., [&](const auto &) {
        if (++probes > 10) throw std::runtime_error("Coverage did not converge");
        return 1.;
    }));
    REQUIRE(cloud.retired_below() >= cutoff);
}

TEST_CASE("Adaptive streaming never resamples rejected candidates", "[FillVoronoi][Incremental]")
{
    std::set<std::array<double, 3>> visited;
    auto provider = std::make_shared<Voronoi::WallDistancePointCloud>([&](const auto &points) {
        std::vector<InfillDensitySample> samples;
        for (const auto &p : points) {
            REQUIRE(visited.insert({p.x(), p.y(), p.z()}).second);
            samples.push_back({-1., 0.5});
        }
        return samples;
    }, 1.);
    Voronoi::ActivePointCloud cloud(provider);
    const auto bounds = rectangle(-2., -2., 2., 2.).contour.bounding_box();
    for (double z = 0.; z < 20.; z += 0.2)
        REQUIRE(cloud.update(bounds, z, 0., 2., [](const auto &) { return 1.; }));
    REQUIRE_FALSE(visited.empty());
}

TEST_CASE("Active cloud uses the provider's generation cell size", "[FillVoronoi][Incremental]")
{
    class WideCells : public Voronoi::PointCloudProvider {
    public:
        double spacing(double, double, double) const override { return 1.; }
        double cell_size(double spacing) const override { return 7. * spacing; }
        std::vector<Vec3d> points_in(const BoundingBoxf3 &, double) const override { return {}; }
        std::vector<Vec3d> points_in_cells(const Vec3i64 &first, const Vec3i64 &last, double spacing) const override {
            std::vector<Vec3d> points;
            for (int64_t x = first.x(); x < last.x(); ++x)
                for (int64_t y = first.y(); y < last.y(); ++y)
                    for (int64_t z = first.z(); z < last.z(); ++z)
                        points.emplace_back((Vec3d(double(x), double(y), double(z)).array() + 0.5).matrix() * cell_size(spacing));
            return points;
        }
    };
    Voronoi::ActivePointCloud cloud(std::make_shared<WideCells>());
    const auto bounds = rectangle(-1., -1., 1., 1.).contour.bounding_box();
    std::vector<double> cutoffs;
    const auto retire = [&](double cutoff) { if (std::isfinite(cutoff)) cutoffs.push_back(cutoff); };
    const auto radius = [](const auto &) { return 1.; };
    REQUIRE(cloud.update(bounds, 0.1, 0., 1., radius, retire));
    REQUIRE_THAT(cloud.known_bounds().max.z(), Catch::Matchers::WithinAbs(7., EPSILON));
    const size_t revision = cloud.revision();
    REQUIRE(cloud.update(bounds, 5., 0., 1., radius, retire));
    REQUIRE(cloud.revision() == revision);
    REQUIRE(cutoffs.empty());
    REQUIRE(cloud.update(bounds, 6.5, 0., 1., radius, retire));
    REQUIRE(cutoffs.size() == 1);
    REQUIRE_THAT(cutoffs.front(), Catch::Matchers::WithinAbs(0., EPSILON));
    REQUIRE_THAT(cloud.known_bounds().max.z(), Catch::Matchers::WithinAbs(14., EPSILON));
    REQUIRE(cloud.update(bounds, 9., 0., 1., radius, retire));
    REQUIRE(cutoffs.size() == 1);
    REQUIRE(cloud.update(bounds, 13.5, 0., 1., radius, retire));
    REQUIRE(cutoffs.size() == 2);
    REQUIRE_THAT(cutoffs.back(), Catch::Matchers::WithinAbs(7., EPSILON));
    REQUIRE_THAT(cloud.known_bounds().max.z(), Catch::Matchers::WithinAbs(21., EPSILON));
    for (const auto &point : cloud.sites()) REQUIRE(point.z() >= cutoffs.back());
}
