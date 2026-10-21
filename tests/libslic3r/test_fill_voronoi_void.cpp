#include <catch2/catch_all.hpp>
#include "libslic3r/Fill/VoronoiPointCloud.hpp"
#include "libslic3r/InfillDensityField.hpp"
#include "catch2/catch_test_macros.hpp"
#include "libslic3r/VoronoiCloudMethod.hpp"
#include "catch2/catch_message.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Fill/VoronoiVoid.hpp"
#include "libslic3r/Fill/VoronoiLattice.hpp"
#include "libslic3r/Fill/VoronoiWallDistance.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Point.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

using namespace Slic3r;
using namespace Slic3r::Voronoi;

namespace {
class FixedCloud : public PointCloudProvider {
public:
    FixedCloud()
    {
        for (double x : {-30., 30.})
            for (double y : {-30., 30.})
                for (double z : {-30., 30.}) sites.emplace_back(x, y, z);
        sites.emplace_back(0., 0., -12.);
        sites.emplace_back(0., 0., 12.);
        sites.emplace_back(0., 0., 0.);
    }
    double spacing(double, double, double) const override { return 5.; }
    std::optional<BoundingBoxf3> extent() const override
        { return BoundingBoxf3(Vec3d::Constant(-30.), Vec3d::Constant(30.)); }
    std::vector<Vec3d> points_in(const BoundingBoxf3 &box, double) const override
    {
        std::vector<Vec3d> result;
        for (const auto &p : sites)
            if ((p.array() >= box.min.array()).all() && (p.array() < box.max.array()).all()) result.push_back(p);
        return result;
    }
private:
    std::vector<Vec3d> sites;
};

const BoundingBoxf3 model(Vec3d::Constant(-15.), Vec3d::Constant(15.));
std::vector<InfillDensitySample> samples(const std::vector<Vec3d> &points)
{
    std::vector<InfillDensitySample> result;
    for (const auto &p : points) {
        const double distance = p.cwiseAbs().maxCoeff() - 15.;
        result.push_back({distance, distance < -5. ? 0.1 : 1.});
    }
    return result;
}
std::vector<Vec3d> all(const PointCloudProvider &cloud)
{
    auto box = *cloud.extent();
    if (!box.defined) return {};
    box.min.array() -= 1e-5;
    box.max.array() += 1e-5;
    return cloud.points_in(box, 5.);
}
}

TEST_CASE("Void filtering removes low density sites and retains exterior neighbours", "[FillVoronoiVoid]")
{
    auto source = std::make_shared<LatticePointCloud>(VoronoiCloudMethod::Cubic);
    const VoidPointCloud cloud(source, model, 5., samples, 0.5, {35., 10., 0});
    REQUIRE(cloud.result().removed > 0);
    REQUIRE(cloud.result().inserted == 0);
    bool exterior = false;
    for (const auto &p : all(cloud)) {
        const auto value = samples({p})[0];
        if (value.distance >= 0.) exterior = true;
        else REQUIRE(value.density >= 0.5);
    }
    REQUIRE(exterior);
}

TEST_CASE("Void refinement reduces large shallow face cost with bounded point insertion", "[FillVoronoiVoid]")
{
    auto source = std::make_shared<FixedCloud>();
    const VoidPointCloud cloud(source, model, 5., samples, 0.5, {35., 10., 16});
    const auto &result = cloud.result();
    INFO("cost " << result.initial_cost << " -> " << result.final_cost << ", added " << result.inserted << ", remaining " << result.remaining_faces);
    REQUIRE(result.initial_cost > 0.);
    REQUIRE(result.inserted > 0);
    REQUIRE(result.inserted <= 16);
    REQUIRE(result.final_cost < result.initial_cost);
    REQUIRE(result.remaining_faces == 0);
    // Refinement may deliberately reintroduce support points below the cutoff.
    bool supported_void = false;
    for (const auto &p : all(cloud))
        if (samples({p})[0].distance < -5.) supported_void = true;
    REQUIRE(supported_void);
}

TEST_CASE("Void point queries reproduce their partitioned constellation", "[FillVoronoiVoid]")
{
    auto source = std::make_shared<PoissonPointCloud>();
    const VoidPointCloud cloud(source, model, 7., samples, 0.5, {35., 10., 2});
    const auto original = all(cloud);
    auto low = *cloud.extent();
    low.min.array() -= 1e-5;
    low.max.array() += 1e-5;
    auto high = low;
    low.max.z() = high.min.z() = 0.;
    auto split = cloud.points_in(low, 7.);
    const auto upper = cloud.points_in(high, 7.);
    split.insert(split.end(), upper.begin(), upper.end());
    REQUIRE(split.size() == original.size());
    for (size_t i = 0; i < original.size(); ++i) REQUIRE(split[i] == original[i]);
}

TEST_CASE("Void preparation propagates cancellation", "[FillVoronoiVoid]")
{
    auto source = std::make_shared<PoissonPointCloud>();
    REQUIRE_THROWS_AS(VoidPointCloud(source, model, 5., samples, 0.5, {},
        [](double fraction) { if (fraction >= 0.1) throw std::runtime_error("cancelled"); }), std::runtime_error);
}

TEST_CASE("Void cutoff equality uses the material field encoding and caches the prepared cloud", "[FillVoronoiVoid]")
{
    const ExPolygons contours {ExPolygon(Points {Point::new_scale(-15., -15.), Point::new_scale(15., -15.),
        Point::new_scale(15., 15.), Point::new_scale(-15., 15.)})};
    const auto slices = [&] { return std::vector<WallSlice>{{-15., 15., &contours}}; };
    WallDistanceCloudCache cache;
    const auto get = [&](double cutoff) {
        return cache.get(slices, {{AdaptiveInfillMode::Disabled, 0., 0.15, 0.15, AdaptiveInfillGradient::Linear}, 5., VoronoiCloudMethod::Cubic, 0, 0., false, -1, cutoff}, {}, {});
    };
    const auto uniform = get(0.);
    REQUIRE_FALSE(uniform->extent().has_value());
    const auto cut = get(0.15);
    const auto *prepared = dynamic_cast<const VoidPointCloud *>(cut.get());
    REQUIRE(prepared != nullptr);
    REQUIRE(prepared->result().removed == 0);
    REQUIRE(get(0.15) == cut);
    REQUIRE(get(0.) == uniform);
}

TEST_CASE("An entirely cut-off finite cloud terminates without expanding", "[FillVoronoiVoid]")
{
    const auto zero = [](const std::vector<Vec3d> &points) {
        return std::vector<InfillDensitySample>(points.size(), {-1., 0.});
    };
    auto provider = std::make_shared<VoidPointCloud>(std::make_shared<FixedCloud>(), model, 5., zero, 0.5);
    ActivePointCloud cloud(provider);
    const BoundingBox bounds(Point::new_scale(-1., -1.), Point::new_scale(1., 1.));
    REQUIRE_FALSE(cloud.update(bounds, 0., 0., 5., [](const auto &) {
        return std::numeric_limits<double>::infinity();
    }));
    REQUIRE(cloud.sites().empty());
}

TEST_CASE("Void filtering retains densities equal to the cutoff after float storage", "[FillVoronoiVoid][Regression]")
{
    auto source = std::make_shared<FixedCloud>();
    const auto constant = [](const std::vector<Vec3d> &points) {
        return std::vector<InfillDensitySample>(points.size(), {-1., double(float(0.7))});
    };
    const VoidPointCloud cloud(source, model, 5., constant, 0.7, {35., 10., 0});
    REQUIRE(cloud.result().removed == 0);
}
