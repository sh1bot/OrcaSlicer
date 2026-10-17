#include <catch2/catch_all.hpp>

#include "libslic3r/ShortestPath.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/Polyline.hpp"
#include <cstddef>
#include <random>
#include <vector>

using namespace Slic3r;

TEST_CASE("Path ordering visits every path and respects direction restrictions", "[ShortestPath]")
{
    const size_t count = GENERATE(size_t(0), size_t(1), size_t(512));
    const int direction = GENERATE(0, 1, 2);
    const bool fixed_start = GENERATE(false, true);
    std::mt19937 random(83719);
    std::uniform_int_distribution<int> coordinate(0, 100000000);
    std::vector<ExtrusionPath> paths;
    for (size_t i = 0; i < count; ++i) {
        const Point first(coordinate(random), coordinate(random));
        const Point last = i % 7 == 0 ? first : Point(coordinate(random), coordinate(random));
        paths.emplace_back(erInternalInfill, 0.08, 0.45f, 0.2f);
        paths.back().polyline = Polyline3(Polyline(first, last));
        if (direction == 2 || (direction == 1 && i % 3 == 0))
            paths.back().set_reverse();
    }
    std::vector<ExtrusionEntity *> entities;
    for (auto &path : paths)
        entities.push_back(&path);
    const Point origin(0, 0);
    const auto order = chain_extrusion_entities(entities, fixed_start ? &origin : nullptr);
    REQUIRE(order.size() == count);
    std::vector<bool> visited(count, false);
    for (const auto &[index, reversed] : order) {
        REQUIRE(index < count);
        REQUIRE_FALSE(visited[index]);
        visited[index] = true;
        if (!paths[index].can_reverse())
            REQUIRE_FALSE(reversed);
    }
}
