#pragma once

#include "Config.hpp"
#include "PrintConfig.hpp"
#include <cmath>
#include <string>

namespace Slic3r {

inline bool supports_density_field(InfillPattern pattern, AdaptiveInfillMode mode, double density)
{
    return density > 0. && (pattern == ipVoronoi ||
        (is_adaptive_infill_pattern(pattern) && mode == AdaptiveInfillMode::SmoothBlend && density < 100.));
}

inline bool supports_density_field(const DynamicConfig &config)
{
    return supports_density_field(config.opt_enum<InfillPattern>("sparse_infill_pattern"),
        config.opt_enum<AdaptiveInfillMode>("adaptive_infill_mode"),
        config.option<ConfigOptionPercent>("sparse_infill_density")->value);
}

inline bool fused_density_modifier(const DynamicPrintConfig &config)
{
    const auto *blend = config.option<ConfigOptionFloat>("infill_modifier_blend");
    return blend && blend->value >= 0.;
}

inline bool density_modifier_option(const std::string &key)
{
    return key == "sparse_infill_density" || key == "infill_modifier_blend" ||
           key == "infill_smoothing_sigma" ||
           key == "adaptive_infill_mode" || key == "adaptive_infill_interior_density" || key == "adaptive_infill_gradient";
}

inline DynamicPrintConfig effective_density_modifier_config(const DynamicPrintConfig &config)
{
    if (!fused_density_modifier(config)) return config;
    DynamicPrintConfig effective;
    for (const auto &key : config.keys())
        if (density_modifier_option(key)) effective.set_key_value(key, config.option(key)->clone());
    return effective;
}

inline bool valid_density_modifier_blend(const DynamicPrintConfig &config)
{
    const auto *blend = config.option<ConfigOptionFloat>("infill_modifier_blend");
    return !blend || std::isfinite(blend->value);
}

} // namespace Slic3r
