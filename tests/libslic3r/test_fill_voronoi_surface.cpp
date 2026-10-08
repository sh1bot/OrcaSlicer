#include <catch2/catch_all.hpp>
#include <array>
#include <algorithm>
#include <cmath>
#include <vector>
#include <limits>
#include <memory>
#include <utility>

#include "libslic3r/Fill/FillVoronoi.hpp"
#include "libslic3r/Fill/VoronoiWallDistance.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/Surface.hpp"
#include "libslic3r/ExtrusionEntityCollection.hpp"

using namespace Slic3r;

#include "libslic3r/BridgeDetector.hpp"

TEST_CASE("Voronoi routes whole native floor zigzags from the previous path endpoint", "[FillVoronoi]")
{
    FillParams params;
    params.flow = Flow(0.45f, 0.2f, 0.4f);
    FillVoronoi filler;
    std::vector<Voronoi::Band> bands;
    for (double y : {100., 0., 50.})
        bands.push_back({ExPolygon{Point::new_scale(-2., y), Point::new_scale(0., y),
                                  Point::new_scale(0., y + 10.), Point::new_scale(-2., y + 10.)},
                         Vec2d(1., 0.), 0.});
    std::vector<Polylines> original;
    for (const auto &band : bands) original.push_back(filler.fill_band(band, params));
    const Point start = original[1].front().first_point();
    Voronoi::Layer layer;
    layer.bands = bands;
    const auto ordered = Voronoi::order_groups(filler.fill_paths(std::move(layer), params), start);
    REQUIRE(ordered.size() == original.size());
    for (size_t i = 0; i < ordered.size(); ++i) {
        const auto &expected = original[std::array<size_t, 3>{1, 2, 0}[i]];
        REQUIRE(ordered[i].paths.size() == expected.size());
        for (size_t j = 0; j < expected.size(); ++j)
            REQUIRE(ordered[i].paths[j].points == expected[j].points);
    }
    const auto travel = [&start](const std::vector<Polylines> &groups) {
        double distance = 0.;
        Point current = start;
        for (const auto &group : groups) {
            distance += (group.front().first_point() - current).cast<double>().norm();
            current = group.back().last_point();
        }
        return distance;
    };
    std::vector<Polylines> ordered_paths;
    for (const auto &group : ordered) ordered_paths.push_back(group.paths);
    REQUIRE(travel(ordered_paths) < travel(original));
}

TEST_CASE("Voronoi mixes native odd and even floor zigzags between wall trails", "[FillVoronoi]")
{
    const size_t rows = GENERATE(3, 4);
    FillParams params;
    params.flow = Flow(0.45f, 0.2f, 0.4f);
    const double width = rows * params.flow.spacing();
    Voronoi::Band band{ExPolygon{Point::new_scale(-width, -10.), Point::new_scale(0., -10.),
                                Point::new_scale(0., 10.), Point::new_scale(-width, 10.)}, Vec2d(1., 0.), 0.};
    FillVoronoi filler;
    const Polylines floor = filler.fill_band(band, params);
    REQUIRE(floor.size() == 1);
    const Point entry = floor.front().first_point(), exit = floor.back().last_point();
    REQUIRE((entry.y() == exit.y()) == (rows % 2 == 0));
    const Point start = entry + Point::new_scale(0., -50.);
    Voronoi::Layer layer;
    layer.bands.push_back(band);
    layer.walls = {Polyline(exit + Point::new_scale(0., 50.), exit), Polyline(start, entry)};
    const auto ordered = Voronoi::order_groups(filler.fill_paths(std::move(layer), params), start);
    REQUIRE(ordered.size() == 3);
    REQUIRE(ordered[0].reversible);
    REQUIRE_FALSE(ordered[1].reversible);
    REQUIRE(ordered[2].reversible);
    REQUIRE(ordered[0].paths.front().last_point() == entry);
    REQUIRE(ordered[1].paths.front().points == floor.front().points);
    REQUIRE(ordered[2].paths.front().first_point() == exit);
}

TEST_CASE("Voronoi fills the full width of narrow floor bands", "[FillVoronoi]")
{
    const double width = GENERATE(0.46, 0.6, 1., 2.2);
    const double height = GENERATE(0.12, 0.2, 0.3);
    FillParams params;
    params.flow = Flow(0.45f, float(height), 0.4f);
    const double rounded_width = std::ceil(width / params.flow.spacing()) * params.flow.spacing();
    const double supported_edge = width - 0.4;
    Voronoi::Band band{ExPolygon{Point::new_scale(supported_edge - rounded_width, -10.), Point::new_scale(supported_edge, -10.),
                                Point::new_scale(supported_edge, 10.), Point::new_scale(supported_edge - rounded_width, 10.)},
                       Vec2d(1., 0.), supported_edge};
    FillVoronoi filler;
    filler.set_bounding_box(BoundingBox(Point::new_scale(-20., -20.), Point::new_scale(20., 20.)));
    auto filled = filler.fill_band(band, params);
    REQUIRE_FALSE(filled.empty());
    REQUIRE(filled.size() == 1);
    REQUIRE(diff_pl(filled, band.polygon).empty());
    const Polyline &zigzag = filled.front();
    double previous_x = std::numeric_limits<double>::max();
    double previous_direction = 0.;
    size_t row_count = 0;
    Polyline first_row;
    for (size_t i = 1; i < zigzag.points.size(); ++i) {
        const Point &a = zigzag.points[i - 1], &b = zigzag.points[i];
        const double direction = unscale_(b.y() - a.y());
        if (std::abs(direction) < 19.)
            continue;
        REQUIRE(a.x() <= previous_x);
        if (row_count > 0) {
            REQUIRE_THAT(unscale_(previous_x - a.x()), Catch::Matchers::WithinAbs(params.flow.spacing(), 1e-5));
            REQUIRE(direction * previous_direction < 0.);
        } else {
            first_row = Polyline(a, b);
        }
        previous_x = a.x();
        previous_direction = direction;
        ++row_count;
    }
    REQUIRE(row_count == size_t(std::ceil(width / params.flow.spacing())));
    Voronoi::Band previous = band;
    previous.polygon.translate(Point::new_scale(width / 2., 0.));
    previous.row_origin += width / 2.;
    auto previous_fill = filler.fill_band(previous, params);
    ExtrusionEntityCollection prior_extrusion;
    extrusion_entities_append_paths(prior_extrusion.entities, std::move(previous_fill), erInternalInfill,
                                    params.flow.mm3_per_mm(), params.flow.width(), params.flow.height());
    ExtrusionPath first(erInternalInfill, params.flow.mm3_per_mm(), params.flow.width(), params.flow.height());
    first.polyline = Polyline3(first_row);
    const ExPolygons unsupported = diff_ex(first.polygons_covered_by_width(), prior_extrusion.polygons_covered_by_width());
    double unsupported_area = 0.;
    for (const auto &polygon : unsupported)
        unsupported_area += polygon.area();
    REQUIRE(unsupported_area / band.polygon.area() < 0.001);
    ExtrusionEntityCollection extrusion;
    extrusion_entities_append_paths(extrusion.entities, std::move(filled), erInternalInfill,
                                    params.flow.mm3_per_mm(), params.flow.width(), params.flow.height());
    const ExPolygons uncovered = diff_ex(ExPolygons{band.polygon}, extrusion.polygons_covered_by_width());
    double area = 0.;
    for (const auto &polygon : uncovered)
        area += polygon.area();
    REQUIRE(area / band.polygon.area() < 0.001);
}

TEST_CASE("Voronoi floor zigzags remain within bands containing holes", "[FillVoronoi]")
{
    FillParams params;
    params.flow = Flow(0.45f, 0.2f, 0.4f);
    const double pitch = params.flow.spacing();
    Voronoi::Band band{ExPolygon{Point::new_scale(-8. * pitch, -10.), Point::new_scale(0., -10.),
                                Point::new_scale(0., 10.), Point::new_scale(-8. * pitch, 10.)}, Vec2d(1., 0.), 0.};
    band.polygon.holes.emplace_back(Points{Point::new_scale(-6. * pitch, -2.), Point::new_scale(-6. * pitch, 2.),
                                           Point::new_scale(-2. * pitch, 2.), Point::new_scale(-2. * pitch, -2.)});
    FillVoronoi filler;
    const Polylines paths = filler.fill_band(band, params);
    REQUIRE_FALSE(paths.empty());
    REQUIRE(diff_pl(paths, band.polygon).empty());
    bool has_join = false;
    for (const Polyline &path : paths) {
        REQUIRE(path.first_point().x() >= path.last_point().x());
        if (path.points.size() > 3)
            has_join = true;
    }
    REQUIRE(has_join);
}

TEST_CASE("Voronoi fills shallow bands with continuous coverage and configured flow", "[FillVoronoi]")
{
    const double height = GENERATE(0.12, 0.2, 0.3);
    ExPolygon region{Point::new_scale(0., 0.), Point::new_scale(80., 0.),
                     Point::new_scale(80., 80.), Point::new_scale(0., 80.)};
    region.holes.emplace_back(Points{Point::new_scale(30., 30.), Point::new_scale(30., 50.),
                                    Point::new_scale(50., 50.), Point::new_scale(50., 30.)});
    FillParams params;
    params.density = 0.2f;
    params.pattern = ipVoronoi;
    params.flow = Flow(0.45f, float(height), 0.4f);
    params.extrusion_role = erInternalInfill;
    params.using_internal_flow = true;
    FillVoronoi filler;
    filler.set_bounding_box(region.contour.bounding_box());
    filler.layer_id = 30;
    filler.z = 6.;
    filler.angle = 0.f;
    filler.spacing = params.flow.spacing();
    Surface surface(stInternal, region);
    const Polylines paths = filler.fill_surface(&surface, params);
    REQUIRE_FALSE(paths.empty());
    REQUIRE(diff_pl(paths, region).empty());

    // Compare the deposited paths with the geometric bands. The real adapter
    // must fill their interiors, rather than merely tracing polygon edges.
    VoronoiInfill geometry;
    const ExPolygons inset = offset_ex(region, -0.5f * float(scale_(filler.spacing)));
    REQUIRE(inset.size() == 1);
    const auto layer = geometry.fill_layer(filler.bounding_box, inset.front(), filler.spacing,
                                          params.density, filler.z, params.flow.width(), params.flow.height(), height, params.flow.spacing());
    REQUIRE_FALSE(layer.bands.empty());
    ExtrusionEntityCollection output;
    filler.fill_surface_extrusion(&surface, params, output.entities);
    REQUIRE(output.entities.size() == 1);
    const auto *ordered = dynamic_cast<const ExtrusionEntityCollection *>(output.entities.front());
    REQUIRE(ordered != nullptr);
    REQUIRE(ordered->no_sort);
    REQUIRE_FALSE(ordered->can_reverse());
    const Polygons deposited = output.polygons_covered_by_width();
    double band_area = 0., uncovered_area = 0.;
    for (const auto &band : layer.bands) {
        band_area += band.polygon.area();
        for (const ExPolygon &uncovered : diff_ex(ExPolygons{band.polygon}, deposited))
            uncovered_area += uncovered.area();
        // The long rows and their short contour joins both belong to this
        // floor; joins may run in a different direction from the rows.
        // Allow fixed-point rounding of vertices on diagonal clip boundaries.
        REQUIRE(diff_pl(filler.fill_band(band, params), offset_ex(band.polygon, float(SCALED_EPSILON))).empty());
    }
    // Clipped tapered tips can be narrower than a row; the rectangular
    // regression above requires complete coverage across the floor width.
    REQUIRE(uncovered_area / band_area < 0.02);
    const double deposited_volume = output.total_volume();
    const double region_volume = region.area() * SCALING_FACTOR * SCALING_FACTOR * height;
    REQUIRE_THAT(deposited_volume / region_volume, Catch::Matchers::WithinRel(double(params.density), 0.2));

    // Rows keep their order but may be printed from either end.
    const auto flattened = output.flatten();
    REQUIRE_FALSE(flattened.entities.empty());
    for (const ExtrusionEntity *entity : flattened.entities) {
        REQUIRE(entity->can_reverse());
        const auto *path = dynamic_cast<const ExtrusionPath *>(entity);
        REQUIRE(path != nullptr);
        REQUIRE_THAT(double(path->width), Catch::Matchers::WithinAbs(double(params.flow.width()), 1e-6));
        REQUIRE_THAT(path->mm3_per_mm, Catch::Matchers::WithinAbs(params.flow.mm3_per_mm(), 1e-9));
    }
}

TEST_CASE("Worker-local Voronoi caches survive temporary fillers and isolate clones", "[FillVoronoi][Incremental]")
{
    class CountingCloud : public Voronoi::PoissonPointCloud {
    public:
        mutable size_t queries = 0;
        std::vector<Vec3d> points_in(const BoundingBoxf3 &box, double spacing) const override
            { ++queries; return PoissonPointCloud::points_in(box, spacing); }
    };
    auto provider = std::make_shared<CountingCloud>();
    VoronoiInfillCache cache;
    FillParams params;
    params.density = 0.25f;
    params.flow = Flow(0.45f, 0.2f, 0.4f);
    params.using_internal_flow = true;
    const ExPolygon outline {Point::new_scale(-10., -10.), Point::new_scale(10., -10.),
                             Point::new_scale(10., 10.), Point::new_scale(-10., 10.)};
    const Surface surface(stInternal, outline);
    const auto configure = [&](FillVoronoi &fill) {
        fill.set_point_cloud(provider, &cache);
        fill.set_bounding_box(outline.contour.bounding_box());
        fill.spacing = params.flow.spacing();
        fill.angle = 0.f;
        fill.z = 1.1;
        fill.layer_id = 5;
    };
    Polylines expected;
    {
        FillVoronoi temporary;
        configure(temporary);
        expected = temporary.fill_surface(&surface, params);
    }
    REQUIRE(provider->queries > 0);
    REQUIRE(cache.clouds.size() == 1);
    const size_t initial_queries = provider->queries;
    FillVoronoi next;
    configure(next);
    REQUIRE(next.fill_surface(&surface, params) == expected);
    REQUIRE(provider->queries == initial_queries);
    std::unique_ptr<Fill> clone(next.clone());
    clone->z = 1000.;
    REQUIRE_FALSE(clone->fill_surface(&surface, params).empty());
    const size_t after_clone = provider->queries;
    REQUIRE(next.fill_surface(&surface, params) == expected);
    REQUIRE(provider->queries == after_clone);
}

TEST_CASE("Voronoi wall distribution reduces material through native extrusion", "[FillVoronoi]")
{
    const ExPolygons outlines {ExPolygon{Point::new_scale(-40., -40.), Point::new_scale(40., -40.),
                                        Point::new_scale(40., 40.), Point::new_scale(-40., 40.)}};
    auto field = std::make_shared<Voronoi::CubicWallDistance>(std::vector<Voronoi::WallSlice>{{-100., 100., &outlines}});
    FillVoronoi uniform, biased;
    const Voronoi::WallDistanceSamples distance = [field](const auto &points) { return field->sample_fields(points); };
    std::shared_ptr<const Voronoi::PointCloudProvider> cloud = std::make_shared<Voronoi::WallDistancePointCloud>(distance);
    biased.set_point_cloud(std::move(cloud));
    FillParams params;
    params.density = 0.2f;
    params.flow = Flow(0.45f, 0.2f, 0.4f);
    params.using_internal_flow = true;
    for (auto *filler : {&uniform, &biased}) {
        filler->set_bounding_box(outlines.front().contour.bounding_box());
        filler->spacing = params.flow.spacing();
        filler->angle = 0.f;
        filler->z = 4.;
        filler->layer_id = 20;
    }
    Surface surface(stInternal, outlines.front());
    ExtrusionEntityCollection original, thinned, repeat;
    uniform.fill_surface_extrusion(&surface, params, original.entities);
    biased.fill_surface_extrusion(&surface, params, thinned.entities);
    REQUIRE_FALSE(thinned.entities.empty());
    REQUIRE(thinned.total_volume() > 0.);
    REQUIRE(thinned.total_volume() < original.total_volume());
    biased.fill_surface_extrusion(&surface, params, repeat.entities);
    REQUIRE_THAT(repeat.total_volume(), Catch::Matchers::WithinAbs(thinned.total_volume(), EPSILON));
}

TEST_CASE("Voronoi inner hull removes interior paths and prints a continuous skin", "[FillVoronoi]")
{
    const ExPolygons outlines {ExPolygon{Point::new_scale(-25., -25.), Point::new_scale(25., -25.),
                                        Point::new_scale(25., 25.), Point::new_scale(-25., 25.)}};
    auto table = std::make_shared<Voronoi::CubicWallDistance>(std::vector<Voronoi::WallSlice>{{0., 50., &outlines}}, 0.);
    size_t samples = 0;
    auto cloud = std::make_shared<Voronoi::WallDistancePointCloud>([table, &samples](const auto &p) {
        samples += p.size();
        return table->sample_fields(p);
    });
    FillParams params;
    params.density = 0.25f;
    params.flow = Flow(0.45f, 0.2f, 0.4f);
    params.using_internal_flow = true;
    FillVoronoi baseline, hollow;
    for (auto *filler : {&baseline, &hollow}) {
        filler->set_point_cloud(cloud);
        filler->set_bounding_box(outlines.front().contour.bounding_box());
        filler->spacing = params.flow.spacing();
        filler->angle = 0.f;
        filler->z = 25.;
        filler->layer_id = 125;
    }
    const Surface surface(stInternal, outlines.front());
    const auto original = baseline.fill_surface(&surface, params);
    REQUIRE(hollow.fill_surface(&surface, params) == original);
    hollow.set_hull_threshold(0.5);
    const auto paths = hollow.fill_surface(&surface, params);
    const size_t before_repeat = samples;
    REQUIRE(hollow.fill_surface(&surface, params) == paths);
    REQUIRE(samples == before_repeat);
    const auto hull = Voronoi::hull_section(*cloud, hollow.bounding_box, 24.9, 0.2, params.flow.width(), 0.5);
    REQUIRE_FALSE(paths.empty());
    REQUIRE(intersection_pl(paths, offset_ex(hull.cavity, -float(scale_(0.05)))).empty());
    bool closed = false;
    for (const auto &path : paths) if (path.first_point() == path.last_point()) closed = true;
    REQUIRE(closed);
    ExtrusionEntityCollection output;
    hollow.fill_surface_extrusion(&surface, params, output.entities);
    const auto deposited = output.polygons_covered_by_width();
    double total = 0., missing = 0.;
    for (const auto &p : hull.skin) total += p.area();
    for (const auto &p : diff_ex(hull.skin, deposited)) missing += p.area();
    REQUIRE(missing / total < 0.02);
    const auto replacement = std::make_shared<Voronoi::PoissonPointCloud>();
    hollow.set_point_cloud(replacement);
    baseline.set_point_cloud(replacement);
    baseline.set_hull_threshold(0.5);
    REQUIRE(hollow.fill_surface(&surface, params) == baseline.fill_surface(&surface, params));
    REQUIRE_FALSE(hollow.fill_surface(&surface, params) == paths);
    const auto flattened = output.flatten();
    for (const auto *entity : flattened.entities) {
        const auto *path = dynamic_cast<const ExtrusionPath *>(entity);
        REQUIRE(path != nullptr);
        REQUIRE_THAT(path->mm3_per_mm, Catch::Matchers::WithinAbs(params.flow.mm3_per_mm(), 1e-9));
    }
}

TEST_CASE("Voronoi hull floors and ceilings receive native solid fill", "[FillVoronoi]")
{
    auto cloud = std::make_shared<Voronoi::WallDistancePointCloud>([](const auto &points) {
        std::vector<Voronoi::WallDistanceSample> samples;
        for (const auto &p : points) {
            const double d = -std::min({20. - std::abs(p.x()), 20. - std::abs(p.y()), p.z(), 20. - p.z()});
            samples.push_back({d, d});
        }
        return samples;
    });
    const ExPolygon outline{Point::new_scale(-20., -20.), Point::new_scale(20., -20.),
                            Point::new_scale(20., 20.), Point::new_scale(-20., 20.)};
    FillParams params;
    params.density = 0.25f;
    params.flow = Flow(0.45f, 0.2f, 0.4f);
    params.using_internal_flow = true;
    for (double z : {6.3, 13.9}) {
        FillVoronoi filler;
        filler.set_point_cloud(cloud);
        filler.set_hull_threshold(0.5);
        filler.set_bounding_box(outline.contour.bounding_box());
        filler.spacing = params.flow.spacing();
        filler.angle = 0.f;
        filler.z = z;
        filler.layer_id = size_t(z / 0.2);
        const Surface surface(stInternal, outline);
        ExtrusionEntityCollection output;
        filler.fill_surface_extrusion(&surface, params, output.entities);
        const ExPolygon center{Point::new_scale(-5., -5.), Point::new_scale(5., -5.),
                                Point::new_scale(5., 5.), Point::new_scale(-5., 5.)};
        const auto missing = diff_ex(ExPolygons{center}, output.polygons_covered_by_width());
        REQUIRE(missing.empty());
    }
}

TEST_CASE("Independent Voronoi hull loops can route through an intervening wall", "[FillVoronoi]")
{
    const auto box = [](double x) {
        return ExPolygon{Point::new_scale(x, 0.), Point::new_scale(x + 10., 0.),
                         Point::new_scale(x + 10., 10.), Point::new_scale(x, 10.)};
    };
    Voronoi::HullSection hull{{box(0.), box(30.)}, {}};
    const ExPolygon region{Point::new_scale(-5., -5.), Point::new_scale(45., -5.),
                           Point::new_scale(45., 15.), Point::new_scale(-5., 15.)};
    FillParams params;
    params.flow = Flow(0.45f, 0.2f, 0.4f);
    FillVoronoi filler;
    filler.set_bounding_box(region.contour.bounding_box());
    auto groups = filler.fill_hull(hull, region, params);
    REQUIRE(groups.size() == 2);
    const Point start = groups[0].paths.front().first_point();
    const Point next = groups[1].paths.front().first_point();
    const Polyline wall(start + Point::new_scale(0., 1.), next);
    groups.push_back({Polylines{wall}, true});
    const auto ordered = Voronoi::order_groups(std::move(groups), start);
    REQUIRE(ordered[0].paths.front().first_point() == start);
    REQUIRE(ordered[1].paths.front().points == wall.points);
    REQUIRE(ordered[2].paths.front().first_point() == next);
}

namespace {
class FiniteWallSites : public Voronoi::PoissonPointCloud {
public:
    std::vector<Vec3d> points_in(const BoundingBoxf3 &box, double) const override
    {
        std::vector<Vec3d> result;
        for (double x : {-5., 5.}) {
            const Vec3d site(x, 0., 0.);
            if ((site.array() >= box.min.array()).all() && (site.array() < box.max.array()).all()) result.push_back(site);
        }
        return result;
    }
    std::optional<BoundingBoxf3> extent() const override
    {
        BoundingBoxf3 box;
        box.merge(Vec3d(-5., 0., 0.));
        box.merge(Vec3d(5., 0., 0.));
        return box;
    }
};
}

TEST_CASE("Voronoi uses native multiline outlines instead of cell loops", "[FillVoronoi][Multiline]")
{
    const int count = GENERATE(1, 2, 3, 4);
    FillParams params;
    params.multiline = count;
    params.density = 0.2f;
    params.flow = Flow(0.45f, 0.2f, 0.4f);
    params.using_internal_flow = true;
    params.anchor_length_max = 0.f;
    const ExPolygon region {Point::new_scale(-10., -10.), Point::new_scale(10., -10.),
                            Point::new_scale(10., 10.), Point::new_scale(-10., 10.)};
    FillVoronoi filler;
    filler.set_point_cloud(std::make_shared<FiniteWallSites>());
    filler.set_bounding_box(region.contour.bounding_box());
    filler.spacing = params.flow.spacing();
    filler.z = 0.;
    filler.angle = 0.f;
    const Surface surface(stInternal, region);
    auto expected = FillVoronoi();
    expected.set_point_cloud(std::make_shared<FiniteWallSites>());
    expected.set_bounding_box(region.contour.bounding_box());
    expected.spacing = filler.spacing;
    expected.z = 0.;
    expected.angle = 0.f;
    FillParams single = params;
    single.multiline = 1;
    Polylines native = expected.fill_surface(&surface, single);
    multiline_fill(native, params, float(params.flow.spacing()));
    const auto clip = offset_ex(region, -float(scale_(0.5 * filler.spacing)));
    native = intersection_pl(std::move(native), clip);
    const auto actual = filler.fill_surface(&surface, params);
    Polygons actual_area = offset(actual, float(scale_(0.5 * params.flow.width())), jtRound, 0.01, etOpenRound);
    Polygons native_area = offset(native, float(scale_(0.5 * params.flow.width())), jtRound, 0.01, etOpenRound);
    REQUIRE(diff(actual_area, native_area).empty());
    REQUIRE(diff(native_area, actual_area).empty());
}

TEST_CASE("Voronoi honours native boundary anchor lengths", "[FillVoronoi][Anchors]")
{
    const ExPolygon region {Point::new_scale(-10., -10.), Point::new_scale(10., -10.),
                            Point::new_scale(10., 10.), Point::new_scale(-10., 10.)};
    FillParams params;
    params.density = 0.2f;
    params.flow = Flow(0.45f, 0.2f, 0.4f);
    FillVoronoi filler;
    filler.set_point_cloud(std::make_shared<FiniteWallSites>());
    filler.set_bounding_box(region.contour.bounding_box());
    filler.spacing = params.flow.spacing();
    filler.z = 0.;
    filler.angle = 0.f;
    const Surface surface(stInternal, region);
    params.anchor_length_max = 0.f;
    const auto unanchored = filler.fill_surface(&surface, params);
    params.anchor_length = params.anchor_length_max = 2.f;
    const auto anchored = filler.fill_surface(&surface, params);
    double original_length = 0., anchored_length = 0.;
    for (const auto &path : unanchored) original_length += path.length();
    for (const auto &path : anchored) anchored_length += path.length();
    REQUIRE(anchored_length > original_length);
    REQUIRE(diff_pl(anchored, offset_ex(region, float(SCALED_EPSILON))).empty());
}

TEST_CASE("Floor rows use existing wall anchors without reversing their sweep", "[FillVoronoi][Floors][Anchors]")
{
    const bool both_ends = GENERATE(false, true);
    FillParams params;
    params.flow = Flow(0.45f, 0.2f, 0.4f);
    const ExPolygon floor {Point::new_scale(-2.5, -10.), Point::new_scale(0., -10.),
                           Point::new_scale(0., 10.), Point::new_scale(-2.5, 10.)};
    ExPolygons support {
        ExPolygon{Point::new_scale(-3., -10.8), Point::new_scale(1., -10.8), Point::new_scale(1., -9.99), Point::new_scale(-3., -9.99)}};
    if (both_ends) support.push_back(ExPolygon{Point::new_scale(-3., 9.99), Point::new_scale(1., 9.99),
                                              Point::new_scale(1., 10.8), Point::new_scale(-3., 10.8)});
    FillVoronoi filler;
    filler.set_point_cloud(std::make_shared<FiniteWallSites>());
    filler.set_bounding_box(floor.contour.bounding_box());
    filler.spacing = params.flow.spacing();
    filler.z = 0.;
    filler.angle = 0.f;
    params.density = 0.2f;
    params.anchor_length_max = 0.f;
    const Surface surface(stInternal, floor);
    filler.fill_surface(&surface, params);
    filler.set_floor_support(support);
    const auto paths = filler.fill_band({floor, Vec2d(1., 0.), 0.}, params);
    REQUIRE_FALSE(paths.empty());
    bool has_rows = false;
    for (const auto &path : paths) {
        REQUIRE(path.first_point().x() >= path.last_point().x());
        for (size_t i = 1; i < path.points.size(); ++i)
            if (std::abs(path.points[i].y() - path.points[i - 1].y()) > scale_(10.)) {
                for (const Point &p : {path.points[i - 1], path.points[i]}) {
                    REQUIRE(expolygons_contain(support, p) == (both_ends || p.y() < 0));
                    if (!both_ends && p.y() > 0) REQUIRE(p.y() <= scale_(10.));
                }
                has_rows = true;
            }
    }
    REQUIRE(has_rows);
}

TEST_CASE("Native bridge coverage requires real support at both ends", "[FillVoronoi][Floors][Bridge]")
{
    const bool both_ends = GENERATE(false, true);
    const ExPolygon floor {Point::new_scale(-2.5, -10.), Point::new_scale(0., -10.),
                           Point::new_scale(0., 10.), Point::new_scale(-2.5, 10.)};
    ExPolygons support {ExPolygon{Point::new_scale(-3., -11.), Point::new_scale(1., -11.),
                                 Point::new_scale(1., -9.5), Point::new_scale(-3., -9.5)}};
    if (both_ends) support.push_back(ExPolygon{Point::new_scale(-3., 9.5), Point::new_scale(1., 9.5),
                                              Point::new_scale(1., 11.), Point::new_scale(-3., 11.)});
    BridgeDetector detector(floor, support, scale_(0.4));
    const auto coverage = detector.coverage(M_PI / 2.);
    double area = 0.;
    for (const auto &polygon : coverage) area += polygon.area();
    if (both_ends) REQUIRE(area > 0.8 * floor.area());
    else REQUIRE_THAT(area, Catch::Matchers::WithinAbs(0., 1.));
}
