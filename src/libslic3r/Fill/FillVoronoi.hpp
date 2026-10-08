#pragma once

#include "FillBase.hpp"
#include "../BoundingBox.hpp"
#include "VoronoiInfill.hpp"
#include "VoronoiRouting.hpp"
#include "VoronoiHull.hpp"
#include <array>
#include <memory>
#include <vector>
#include <utility>

namespace Slic3r {

class FillAlignedRectilinear;

class FillVoronoi : public Fill
{
public:
    FillVoronoi() : m_provider(std::make_shared<Voronoi::PoissonPointCloud>()), m_infill(std::make_shared<VoronoiInfill>(m_provider)) {}
    FillVoronoi(const FillVoronoi &other)
        : Fill(other), m_surface_boundary(other.m_surface_boundary), m_floor_support(other.m_floor_support), m_hull(other.m_hull), m_hull_bounds(other.m_hull_bounds),
          m_hull_key(other.m_hull_key), m_provider(other.m_provider), m_infill(std::make_shared<VoronoiInfill>(*other.m_infill)), m_hull_threshold(other.m_hull_threshold) {}
    FillVoronoi &operator=(const FillVoronoi &other)
        { if (this != &other) *this = FillVoronoi(other); return *this; }
    FillVoronoi &operator=(FillVoronoi &&) = default;
    void set_point_cloud(std::shared_ptr<const Voronoi::PointCloudProvider> provider, VoronoiInfillCache *cache = nullptr)
    {
        m_provider = provider;
        m_hull_key = {};
        if (cache) {
            auto &engine = cache->clouds[provider];
            if (!engine) engine = std::make_shared<VoronoiInfill>(provider);
            m_infill = engine;
        } else {
            m_infill = std::make_shared<VoronoiInfill>(std::move(provider));
        }
    }
    void set_hull_threshold(double threshold) { m_hull_threshold = threshold; }
    void set_floor_support(ExPolygons support) { m_floor_support = std::move(support); }
    Polylines fill_surface(const Surface *surface, const FillParams &params) override;
    Fill *clone() const override { return new FillVoronoi(*this); }
    bool no_sort() const override { return true; }
    bool is_self_crossing() override { return false; }
    Polylines fill_band(const Voronoi::Band &band, const FillParams &params) const;
    std::vector<Voronoi::PathGroup> fill_hull(const Voronoi::HullSection &hull, const ExPolygon &region, const FillParams &params) const;
    std::vector<Voronoi::PathGroup> fill_paths(Voronoi::Layer layer, const FillParams &params) const;

protected:
    float _layer_angle(size_t) const override { return 0.f; }
    void _fill_surface_single(const FillParams &params, unsigned int thickness_layers,
                              const std::pair<float, Point> &direction, ExPolygon expolygon,
                              Polylines &polylines_out) override;

private:
    Polylines solid_fill(FillAlignedRectilinear &filler, const ExPolygon &polygon,
                        FillParams params, bool monotonic = false) const;
    ExPolygon m_surface_boundary;
    ExPolygons m_floor_support;
    Voronoi::HullSection m_hull;
    BoundingBox m_hull_bounds;
    std::array<double, 4> m_hull_key {};
    std::shared_ptr<const Voronoi::PointCloudProvider> m_provider;
    std::shared_ptr<VoronoiInfill> m_infill;
    double m_hull_threshold { 0. };
};

} // namespace Slic3r
