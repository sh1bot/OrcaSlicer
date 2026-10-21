#include <vector>
#include "catch2/catch_test_macros.hpp"
#include "catch2/generators/catch_generators.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Flow.hpp"
#include "libslic3r/ExPolygon.hpp"
#include <catch2/catch_all.hpp>
#include "libslic3r/Fill/FillBase.hpp"
#include "libslic3r/Fill/AdaptiveInfillField.hpp"
#include "libslic3r/InfillDensityField.hpp"
#include "libslic3r/Surface.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include <memory>

using namespace Slic3r;

TEST_CASE("TPMS fillers consume the shared material field for reinforced and empty regions", "[AdaptiveInfillField][InfillDensityField]")
{
    const auto pattern = GENERATE(ipGyroid, ipTpmsD, ipTpmsFK);
    const double target = GENERATE(0., 0.5);
    const bool body_on_grid = GENERATE(false, true);
    const ExPolygons model {ExPolygon(Points {Point::new_scale(-20., -20.), Point::new_scale(20., -20.),
        Point::new_scale(20., 20.), Point::new_scale(-20., 20.)})};
    const ExPolygon center(Points {Point::new_scale(-10., -10.), Point::new_scale(10., -10.),
        Point::new_scale(10., 10.), Point::new_scale(-10., 10.)});
    DensityModifier modifier {3, target, 3., {{0., 40., model}}};
    if (target > 0.) modifier.interior_density = target;
    const InfillDensityField material({{0., 40., &model}}, 0., {}, AdaptiveInfillMode::Disabled, 0.1, {modifier}, 3);
    const ExPolygons empty;
    const std::vector<AdaptiveInfillField::Slice> slices {{0., 40., body_on_grid ? &model : &empty}};
    const AdaptiveInfillField radial(slices, get_extents(model), AdaptiveInfillMode::SmoothBlend, [] {});
    REQUIRE(radial.empty() == !body_on_grid);
    FillParams params;
    params.density = 0.1f;
    params.layer_height = 0.2;
    params.flow = Flow(0.45f, 0.2f, 0.4f);
    params.adaptive_infill_mode = AdaptiveInfillMode::SmoothBlend;
    params.adaptive_infill_interior_density = 0.05f;
    params.anchor_length = 0.;
    params.anchor_length_max = 0.;
    std::unique_ptr<Fill> fill(Fill::new_from_type(pattern));
    fill->set_bounding_box(get_extents(model));
    fill->z = 20.;
    fill->angle = 0.3f;
    fill->spacing = params.flow.spacing();
    fill->adaptive_infill_field = &radial;
    fill->material_density_field = &material;
    const Surface surface(stInternal, center);
    const auto result = fill->fill_surface(&surface, params);
    if (target == 0.) REQUIRE(result.empty());
    else {
        REQUIRE_FALSE(result.empty());
        fill->material_density_field = nullptr;
        double dense_length = 0., ordinary_length = 0.;
        for (const auto &line : result) dense_length += line.length();
        for (const auto &line : fill->fill_surface(&surface, params)) ordinary_length += line.length();
        REQUIRE(dense_length > 2. * ordinary_length);
    }
}
