#pragma once

#include <algorithm>

namespace Slic3r {

inline double modifier_blend(double inside_distance, double width)
{
    if (inside_distance < 0.) return 0.;
    if (width == 0.) return 1.;
    const double t = std::clamp(inside_distance / width, 0., 1.);
    return t * t * (3. - 2. * t);
}

} // namespace Slic3r
