#include "catch2/generators/catch_generators.hpp"
#include "libslic3r/InfillDensityModifier.hpp"
#include <limits>
// Regression test for the "option in def + UI but missing from preset key list"
// crash class.
//
// The print preset's DynamicPrintConfig is seeded with only the keys returned by
// Preset::print_options() (PresetBundle.cpp). A field added to PrintRegionConfig
// or PrintObjectConfig and registered via print_config_def plus a TabPrint
// optgroup, but left out of print_options(), still gets its control built; on tab
// activation reload_config -> get_config_value dispatches to opt_bool/opt_int on a
// DynamicPrintConfig with no entry for the key, and the accessor null-derefs the
// result of option<T>(key).
//
// The invariant asserted here is the inverse: every key declared on
// PrintRegionConfig and PrintObjectConfig appears in Preset::print_options() or
// Preset::filament_options(), the two preset key lists that seed a print preset's
// DynamicConfig.

#include <catch2/catch_all.hpp>
#include "libslic3r/VoronoiCloudMethod.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"
#include "libslic3r/libslic3r.h"
#include "libslic3r/Preset.hpp"
#include "libslic3r/PrintConfig.hpp"

#include "libslic3r/Config.hpp"
#include <set>
#include <string>

using namespace Slic3r;

namespace {

// Deprecated keys renamed in handle_legacy() (ironing_direction ->
// ironing_angle, wall_infill_order -> wall_sequence); neither is in a
// preset list. Register new options in a preset list, not here.
const std::set<std::string> kDeprecatedRegionFields = {
    "ironing_direction",
    "wall_infill_order",
};

void check_keys_are_in_a_preset(const t_config_option_keys& keys, const std::string& class_name)
{
    REQUIRE_FALSE(keys.empty());
    const auto& print_options    = Preset::print_options();
    const auto& filament_options = Preset::filament_options();
    const std::set<std::string> in_print(print_options.begin(), print_options.end());
    const std::set<std::string> in_filament(filament_options.begin(), filament_options.end());
    for (const std::string& key : keys) {
        DYNAMIC_SECTION(class_name << "::" << key)
        {
            INFO("'" << key << "' on " << class_name
                     << " is missing from "
                        "Preset::print_options()/filament_options(); add it to "
                        "s_Preset_print_options (or s_Preset_filament_options) in Preset.cpp.");
            const bool registered = in_print.count(key) || in_filament.count(key) || kDeprecatedRegionFields.count(key);
            REQUIRE(registered);
        }
    }
}

} // namespace

// Bodies are laid out like the rest of the test suite rather than collapsed
// onto the brace line.
// clang-format off
TEST_CASE("Every PrintRegionConfig field is registered in a preset key list", "[Preset][Config]")
{
    check_keys_are_in_a_preset(PrintRegionConfig::defaults().keys(), "PrintRegionConfig");
}

TEST_CASE("Every PrintObjectConfig field is registered in a preset key list", "[Preset][Config]")
{
    check_keys_are_in_a_preset(PrintObjectConfig::defaults().keys(), "PrintObjectConfig");
}
// clang-format on

TEST_CASE("Print presets initialize and retain shared Voronoi density settings and smoothing", "[Preset][Config]")
{
    DynamicPrintConfig config;
    config.apply_only(FullPrintConfig::defaults(), Preset::print_options());
    REQUIRE_THAT(config.opt_float("infill_modifier_blend"), Catch::Matchers::WithinAbs(-1., EPSILON));
    config.set_deserialize_strict("infill_modifier_blend", "4");
    REQUIRE_THAT(config.opt_float("infill_modifier_blend"), Catch::Matchers::WithinAbs(4., EPSILON));
    REQUIRE_THAT(config.opt_float("adaptive_infill_interior_density"), Catch::Matchers::WithinAbs(5., EPSILON));
    REQUIRE(config.opt_enum<AdaptiveInfillGradient>("adaptive_infill_gradient") == AdaptiveInfillGradient::Linear);
    REQUIRE_THAT(config.opt_float("infill_smoothing_sigma"), Catch::Matchers::WithinAbs(0., EPSILON));
    REQUIRE(config.opt_enum<VoronoiCloudMethod>("voronoi_cloud_method") == VoronoiCloudMethod::Random);
    REQUIRE(config.opt_int("voronoi_relaxation_iterations") == 0);
    config.set_deserialize_strict("voronoi_relaxation_iterations", "5");
    REQUIRE(config.opt_int("voronoi_relaxation_iterations") == 5);
    config.set_deserialize_strict("voronoi_cloud_method", "fcc");
    REQUIRE(config.opt_enum<VoronoiCloudMethod>("voronoi_cloud_method") == VoronoiCloudMethod::Fcc);
    config.set_deserialize_strict("voronoi_cloud_method", "cubic");
    REQUIRE(config.opt_enum<VoronoiCloudMethod>("voronoi_cloud_method") == VoronoiCloudMethod::Cubic);
    REQUIRE(config.option("voronoi_cloud_method")->serialize() == "cubic");
    REQUIRE_FALSE(config.opt_bool("voronoi_trace_shell"));
    config.set_deserialize_strict("voronoi_trace_shell", "1");
    REQUIRE(config.opt_bool("voronoi_trace_shell"));
    REQUIRE_THAT(config.option<ConfigOptionPercent>("voronoi_hull_threshold")->value, Catch::Matchers::WithinAbs(0., EPSILON));
    REQUIRE_FALSE(config.opt_bool("voronoi_snap_to_wall"));
    config.set_deserialize_strict("voronoi_snap_to_wall", "1");
    REQUIRE(config.opt_bool("voronoi_snap_to_wall"));
    config.set_deserialize_strict("adaptive_infill_interior_density", "15%");
    config.set_deserialize_strict("infill_smoothing_sigma", "12");
    REQUIRE_THAT(config.opt_float("adaptive_infill_interior_density"), Catch::Matchers::WithinAbs(15., EPSILON));
    REQUIRE_THAT(config.opt_float("infill_smoothing_sigma"), Catch::Matchers::WithinAbs(12., EPSILON));
}

TEST_CASE("Fused modifiers keep density controls and preserve inactive overrides", "[Preset][FillVoronoi]")
{
    const double blend = GENERATE(-1., 0., 4.);
    DynamicPrintConfig config;
    config.set_deserialize_strict({{"sparse_infill_density", "50%"},
        {"infill_smoothing_sigma", "4"}, {"wall_loops", "3"}, {"fill_multiline", "2"}, {"sparse_infill_pattern", "cubic"}});
    config.set_key_value("infill_modifier_blend", new ConfigOptionFloat(blend));
    const auto effective = Slic3r::effective_density_modifier_config(config);
    REQUIRE(Slic3r::fused_density_modifier(effective) == (blend >= 0.));
    REQUIRE_THAT(effective.opt_float("sparse_infill_density"), Catch::Matchers::WithinAbs(50., EPSILON));
    REQUIRE_THAT(effective.opt_float("infill_smoothing_sigma"), Catch::Matchers::WithinAbs(4., EPSILON));
    for (const std::string key : {"wall_loops", "fill_multiline", "sparse_infill_pattern"}) {
        REQUIRE(effective.has(key) == (blend < 0.));
        REQUIRE(config.has(key));
    }
    REQUIRE(config.opt_int("wall_loops") == 3);
    REQUIRE(config.opt_int("fill_multiline") == 2);
    config.set_key_value("infill_modifier_blend", new ConfigOptionFloat(-1.));
    const auto restored = Slic3r::effective_density_modifier_config(config);
    REQUIRE(restored.opt_int("wall_loops") == 3);
    REQUIRE(restored.opt_int("fill_multiline") == 2);
    REQUIRE(restored.opt_enum<InfillPattern>("sparse_infill_pattern") == ipCubic);
    REQUIRE(Slic3r::valid_density_modifier_blend(config));
    config.set_key_value("infill_modifier_blend", new ConfigOptionFloat(std::numeric_limits<double>::infinity()));
    REQUIRE_FALSE(Slic3r::valid_density_modifier_blend(config));
}

TEST_CASE("Modifier fusion eligibility reads percentage density from mutable parent settings", "[PresetOptions][Regression]")
{
    const auto pattern = GENERATE(ipVoronoi, ipGyroid, ipTpmsD, ipTpmsFK, ipCubic);
    const auto mode = GENERATE(AdaptiveInfillMode::Disabled, AdaptiveInfillMode::SmoothBlend);
    const double density = GENERATE(0., 15., 100.);
    DynamicConfig parent;
    parent.set_key_value("sparse_infill_pattern", new ConfigOptionEnumGeneric(nullptr, int(pattern)));
    parent.set_key_value("adaptive_infill_mode", new ConfigOptionEnumGeneric(nullptr, int(mode)));
    parent.set_key_value("sparse_infill_density", new ConfigOptionPercent(density));
    const bool expected = density > 0. && (pattern == ipVoronoi ||
        (is_adaptive_infill_pattern(pattern) && mode == AdaptiveInfillMode::SmoothBlend && density < 100.));
    REQUIRE(supports_density_field(parent) == expected);
}

TEST_CASE("Earlier development density keys load into the shared infill settings", "[Preset][Config]")
{
    DynamicPrintConfig config;
    config.set_deserialize_strict({{"voronoi_smoothing_sigma", "3"}, {"voronoi_modifier_blend", "2"}});
    REQUIRE_THAT(config.opt_float("infill_smoothing_sigma"), Catch::Matchers::WithinAbs(3., EPSILON));
    REQUIRE_THAT(config.opt_float("infill_modifier_blend"), Catch::Matchers::WithinAbs(2., EPSILON));
}
